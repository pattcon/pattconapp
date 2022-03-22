from Bio import SeqIO
import re


class trieNo(object):

    def __init__(self, chtr):
        self.chr = chtr
        self.seq = []
        self.endOfWord = False
        self.seqCount = 0
        self.aminoCount = 1
        self.contentPerc = 60
        self.conservPerc = 60
        self.wordsize = 0
        self.fileName = ''
        self.fileFormat = 'fasta'
        self.specNumb = 0
        self.children = [None] * 31
        self.ranking = []
        self.isLeaf = False

    


    #Read the sequences from Fasta file
    def readSequences(self):
        listSequences = []
        for i in SeqIO.parse(self.fileName, self.fileFormat):
            listSequences.append(list(i.seq))
        return listSequences


    #Read the total of the residues from the species (number of columns)
    def readMaxSizeSequences(self, listSequencesInit):
        seq = self.readSequences()
        tam = len(seq[0])
        self.numTotalResiduos = tam
        return seq




    #Extract the name of the species from FASTA and build a list with this
    def readSpecies(self):
        listSpecies = []

        for i in SeqIO.parse(self.fileName, self.fileFormat):
            listTemp = []
            listTemp.append(re.split(r'_', i.id))

            if len(listTemp[0][0]) < 3:
                it = listTemp[0][0] + ' ' + listTemp[0][1]
            else:
                it = listTemp[0][0]

            listSpecies.append(it)

        return listSpecies


    #Count the total of the species present in the FASTA file
    def countSpecies(self):
        num = len(self.readSpecies())
        self.specNumb = num
        self.ranking = [0] * num
        return num



    #Calculate the amount of species in a column with no gap's to be consider as a conserved region
    def calcContPercent(self):
        return int(round(self.specNumb * self.contentPerc / 100))

    def calcConservPerc(self):
        return int(round(self.specNumb * self.conservPerc / 100))




    #Eliminate the columns with only gap's from sequências and verify if the column has minimun content
    #to be considered as an conserved region.
    def conservLocals(self, contPercent):
        listSequencesInit = self.readSequences()
        listConservLocals = []

        for j in range(len(listSequencesInit[0])):
            cont = 0
            colunas = []
            for k in range(len(listSequencesInit)):
                if (listSequencesInit[k][j] != '-'):
                    colunas.append(listSequencesInit[k][j])


            if len(colunas) >= self.calcContPercent():
                listConservLocals.append(j)

            colunas.clear()

        return listConservLocals




    #Grouping consecutive conserved regions (columns), return indexes from original list
    def groupConsLocals(self, listConservedLocals):
        listTemp = []
        listConservedFinal = []
        m = 0
        incr = 0
        size = len(listConservedLocals)

        while incr < size:
            chave = incr
            while listConservedLocals[incr] - listConservedLocals[chave] == m:
                listTemp.append(listConservedLocals[incr])
                m = m + 1
                incr = incr + 1
                if incr == size:
                    break
            listConservedFinal.append(listTemp.copy())
            listTemp.clear()
            m = 0
        return listConservedFinal





    #Create a line list from column list of data
    def genConserved(self, listConservedAminos, listSequencesInit):
        listTemp = []
        listAdd = []
        listLocal = []

        for it in listConservedAminos:
            min = it[0]
            size = len(it) - 1
            max = it[size]
            for p in range(len(listSequencesInit)):
                listTemp.append(listSequencesInit[p][min:max + 1])

            listAdd.append(listTemp.copy())

            listTemp.clear()

        #Concatenate simbols of each line, creating words to be used in the trie structure.
        lTemp = []
        listStr = []
        separator = ''
        for it in listAdd:
            for j in it:
                strAm = separator.join(j)
                lTemp.append(strAm)
            listStr.append(lTemp.copy())
            lTemp.clear()

        return listStr











    #Trie structure, modified add procedure.
    def addNode(self, key):

        curr = self
        exist = True
        gapsCount = 0
        empty = False

        for i in range(len(key)):
            if key[i] == '-':
                gapsCount = gapsCount + 1
        if gapsCount == len(key):
            empty = True
            return 0, True


        for i in range(len(key)):
            index = ord(key[i]) - ord('A')

            if index == -20:
                index = 29

            if curr.children[index] is None:
                curr.children[index] = trieNo(index)
                exist = False

            curr = curr.children[index]

            if exist == True:
                curr.seqCount = curr.seqCount + 1
        curr.isLeaf = True

        return curr.seqCount, empty








    # Trie structure, modified add procedure to generate the modifications find in each sequences
    # in compairson with the motif.
    def addNodeGenModif(self, key):

        listOut = []
        curr = self

        empty = False
        idxNew = 0
        idxSubst = -1

        #eliminate sequences with only GAP's
        cont = 0
        for k in range(len(key)):
            if key[k] == "-":
                cont = cont+1
        if cont == len(key)-1:
            empty = True



        for i in range(len(key)):
            exist = True

            if idxSubst == -1:
                index = ord(key[i]) - ord('A')
            else:
                index = idxSubst

            if index == -20:
                index = 29

            if curr.children[index] is None:
                curr.children[index] = trieNo(index)
                exist = False
                listOut.append("[")
                listOut.append(key[i])
                listOut.append("]")
                idxNew = index

                for j in range(31):
                    if curr.children[j] is not None and j != idxNew :
                        idxSubst = j





            if idxSubst == -1:
                curr = curr.children[index]
            else:
                curr = curr.children[idxSubst]


            if exist == True:
                curr.seqCount = curr.seqCount + 1
                listOut.append(key[i])

        curr.isLeaf = True

        return curr.seqCount, empty, listOut











    # Generate the most found sequences in each group
    def genMotifs(self, listSequences, listConserved):

        consLev = self.calcConservPerc()

        listMotifsFinal = []
        listCountTemp = []
        listLocalsFinal = []
        listConsAminFinal = []
        listCountFinal = []
        listPercentFinal = []
        listMotifsNumbers = []
        listOrig = []
        listItTemp = []

        listCountAcum = []
        listItAcum=[]

        idx = 0

        for group in listSequences:
            maxVal = 0
            seqPred = ''
            node = trieNo('+')
            currVal = 0

            for it in group:
                insertTree = node.addNode(it)
                currVal = insertTree[0]+1
                if insertTree[1] == False:

                    listCountTemp.append(currVal)
                    listItTemp.append(it)
                else:
                    currVal = 0
                    listCountTemp.append(currVal)
                    listItTemp.append(it)

            #print("Count occurrences ::::: ", listItTemp , ":::::", listCountTemp)

            valComp = 0
            valMaj = 0
            item = ""
            listaTempUniq = list(sorted(set(listItTemp), key=listItTemp.index))
            listTempUse = []



            for j in listaTempUniq:
                for k in range(len(listCountTemp)):
                    if j == listItTemp[k]:
                        listTempUse.append(listCountTemp[k])
                valMaj = max(listTempUse)
                #print("Compiled of each sequence occurrence :: ", j , " :: ", valMaj)
                listCountAcum.append(valMaj)
                listItAcum.append(j)
                listTempUse.clear()

                valMaj = 0
            listaTempUniq.clear()


            for k in range(len(listCountAcum)):

                if listCountAcum[k] >0 and listCountAcum[k] >= consLev:
                    listMotifsFinal.append(listItAcum[k])
                    listLocalsFinal.append(listConserved[idx])
                    listConsAminFinal.append(group)
                    listCountFinal.append(listCountAcum[k])
                    listPercentFinal.append(round(listCountAcum[k] / self.specNumb * 100))
                    listOrig.append(group)


            for p in range(len(listMotifsFinal)):
                listMotifsNumbers.append(p)

            idx = idx + 1

            listItTemp.clear()
            listCountTemp.clear()
            listItAcum.clear()
            listCountAcum.clear()

            #print("LISTA DE MOTIVOS FINAL :: ", listMotifsFinal)

        return listMotifsFinal, listLocalsFinal, listConsAminFinal, listCountFinal, listPercentFinal, listMotifsNumbers, listOrig














    #Insertion of data in a trie strucuture, generating the motifs and the regions with alterations
    def processSequences(self, listMotifs, listConservedAmin, listConservLocals, listPercOccurr):
        listStrTemp = []
        listStrFinal = []
        conserved = 1
        listTypTemp = []
        listTypeFinal = []
        listConsLocalsFinal = []
        listAltCountTemp = []
        listAltCountFinal = []
        listMotifsFinal = []
        listPercentOccurrFinal = listPercOccurr




        for i in range(len(listConservedAmin)):
            listConsLocalsFinal.append(listConservLocals[i])
            listMotifsFinal.append(listMotifs[i])
            

            for j in range(len(listConservedAmin[i])):
                word = ""
                conserved = 1
                alterations = 0
                for k in range(len(listConservedAmin[i][j])):
                    if listConservedAmin[i][j][k] == listMotifs[i][k]:
                        word = word + listConservedAmin[i][j][k]

                    else:

                        word1 = "["
                        word2 = listConservedAmin[i][j][k]
                        word3 = "]"
                        word = word + word1+word2+word3
                        conserved = 0
                        alterations = alterations + 1


                listAltCountTemp.append(alterations)

                listStrTemp.append(word)
                listTypTemp.append(conserved)


            listTypeFinal.append(listTypTemp.copy())
            listAltCountFinal.append(listAltCountTemp.copy())
            listStrFinal.append(listStrTemp.copy())
            listStrTemp.clear()
            listAltCountTemp.clear()
            listTypTemp.clear()




        return [listStrFinal, listConsLocalsFinal, listMotifsFinal, listAltCountFinal, listTypeFinal, listPercentOccurrFinal]





    def motifBySize(self, finalList, indexes, minSize, txConser, txContent):
        minConser = txConser / 100 * self.countSpecies()
        minContent = txConser / 100 * self.countSpecies()

        sizeSeqVar = 0
        strSeq = []
        sizeSeq = []
        groupSeq = []
        groupId = -1
        mersList = []
        motifsFinal = []
        sizesFinal = []
        suportFinal = []
        localsTemp = []
        localFinal = []

        listStrTemp = []
        listStrFinal = []
        listLocalFinal = []
        localOutTemp = []
        localOut = []



        #print("Indexes ::", finalList)

        for group in finalList:


            groupId = groupId + 1
            groupList = []




            for item in group:




                size = len(item)
                currVal = 0

                incr = -1

                #Mers generation
                for count in range(size):

                    incr = incr + 1

                    for leng in range(size, count, -1):
                        sizeSeqVar = leng-count

                        #Minimum size defined by user
                        if sizeSeqVar>=minSize:
                            strSeq.append(item[count:leng])
                            #List of size of each MER
                            sizeSeq.append(len(item[count:leng]))
                            groupSeq.append(groupId)
                            localsTemp.append(indexes[groupId][incr]+1)

                groupList.extend(strSeq.copy())
                localFinal.extend(localsTemp.copy())




                listStrFinal.append(strSeq.copy())
                listLocalFinal.append(localsTemp.copy())







                strSeq.clear()
                localsTemp.clear()



                listTree = []
                listTreeTemp = []
                maxTam = max(sizeSeq)

                listLocalsTemp = []
                listLocalsFF = []
                #Generate the list of sequence by sizes to be used in the tree
                for j in range(maxTam, minSize-1, -1):
                    for i in range(len(groupList)):
                        if sizeSeq[i] == j:
                            listTreeTemp.append(groupList[i])
                            #listLocalsFF.append(localsTemp[i])


                    listTree.append(listTreeTemp.copy())
                    listTreeTemp.clear()


            for i in range(len(listTree)):
                node = trieNo("+")
                contMer = 0
                majMotifSize = 0
                seqMaj = ""
                sizeMaj = 0
                lclTemp = 0

                for j in range(len(listTree[i])):
                    value = listTree[i][j]
                    contMer = node.addNode(value)[0]+1

                    if contMer>majMotifSize:
                        majMotifSize = contMer
                        seqMaj = value
                        sizeMaj = sizeSeq[i]





                #print("Lista Final ::: ", len(listStrFinal), "lista locals :: ", len(listLocalFinal))



                if majMotifSize>=minConser:

                    motifsFinal.append(seqMaj)
                    sizesFinal.append(sizeMaj)
                    suportFinal.append(majMotifSize)

                    for i in range(len(listStrFinal)):
                        for j in range(len(listStrFinal[i])):
                            if listStrFinal[i][j] == seqMaj:
                                localOutTemp.append(listLocalFinal[i][j])
                                listStrTemp.append(seqMaj)
                                #print("Motivo ::", listStrFinal[i][j], "Local ", listLocalFinal[i][j])

            listStrFinal.clear()
            listLocalFinal.clear()

            #localOutTemp = sorted(set(localOutTemp), key=localOutTemp.index)

            for i in range(len(motifsFinal)):
                for j in range(len(localOutTemp)):
                    if motifsFinal[i] == listStrTemp[j]:
                        localOut.append(localOutTemp[j])
                        break
                    break


        #loc = set(localOut)
        print("Motifs ::: ", motifsFinal)
        print("Size ::: ", sizesFinal)
        print("Suporte ::: ", suportFinal)
        print("Locals ::: ", localOutTemp)


        return suportFinal, motifsFinal,localOutTemp, sizesFinal







    def searchInGroup(self, finalList, group,groupNumber, sizeItem):

        listCounterMers = []
        listaTempMers = []
        listTempLocals = []

        listLocals = []
        listLocalsFinal = []
        listMers = []



        no = trieNo("+")




        for item in group:
            incr = 0
            max = len(item)-sizeItem+1
            for x in range(max):
                mer = item[x: x + sizeItem]
                listaTempMers.append(mer)
                listTempLocals.append(finalList[groupNumber][incr])
                incr = incr+1


                #listaTempMers = sorted(set(listaTempMers), key=listaTempMers.index)

        x=-1
        for mer in listaTempMers:
            listCounterMers.append(no.addNode(mer)[0] + 1)
            listMers.append(mer)
            localMer = listTempLocals[x+1]
            listLocals.append(listTempLocals[x+1])
            x = x + 1


        maxOccurMer = 0
        listMaxMers = []



        for i in range(len(listCounterMers)):
            if listCounterMers[i] > maxOccurMer:
                maxOccurMer = listCounterMers[i]

        for i in range(len(listCounterMers)):
            if listCounterMers[i] == maxOccurMer:
                listMaxMers.append(listMers[i])
                listLocalsFinal.append(listLocals[i])



        return maxOccurMer, listMaxMers, listLocalsFinal


