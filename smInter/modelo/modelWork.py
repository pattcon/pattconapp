from Bio import SeqIO
import re
import math


class trieNo(object):

    def __init__(self, chtr):
        self.chr = chtr
        self.seq = []
        self.endOfWord = False
        self.seqCount = 1
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

        self.listMotifs = []
        self.support = []
        self.localMotifs = []
        self.sizesMotifs = []
        self.linesMotifs = []
        self.groupMotifs = []
        self.beginMotifs =[]
        self.endMotifs = []


    # Read the sequences from Fasta file
    def readSequences(self):
        listSequences = []
        for i in SeqIO.parse(self.fileName, self.fileFormat):
            listSequences.append(list(i.seq))
        return listSequences

    # Read the total of the residues from the species (number of columns)
    def readMaxSizeSequences(self, listSequencesInit):
        seq = self.readSequences()
        tam = len(seq[0])
        self.numTotalResiduos = tam
        return seq

    # Extract the name of the species from FASTA and build a list with this
    def readSpecies(self):
        listSpecies = []

        for i in SeqIO.parse(self.fileName, self.fileFormat):
            listTemp = []
            listTemp.append(re.split(r'_', i.id))

            it = []
            it.append(listTemp[0][0])

            # if len(listTemp[0][0]) < 3:
            # it = listTemp[0][0] + ' ' + listTemp[0][1]
            # else:
            # it = listTemp[0][0]

            listSpecies.append(it)

        return listSpecies

    # Count the total of the species present in the FASTA file
    def countSpecies(self):
        num = len(self.readSpecies())
        self.specNumb = num
        self.ranking = [0] * num
        return num

    # Calculate the amount of species in a column with no gap's to be consider as a conserved region
    def calcContPercent(self):
        return int(round(self.specNumb * self.contentPerc / 100))

    def calcConservPerc(self):
        return int(round(self.specNumb * self.conservPerc / 100))

    # Eliminate the columns with only gap's from sequências and verify if the column has minimun content
    # to be considered as an conserved region.
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

    # Grouping consecutive conserved regions (columns), return indexes from original list
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

    # Create a line list from column list of data
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

        # Concatenate simbols of each line, creating words to be used in the trie structure.
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

    # Trie structure, modified add procedure.
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

        # eliminate sequences with only GAP's
        cont = 0
        for k in range(len(key)):
            if key[k] == "-":
                cont = cont + 1
        if cont == len(key) - 1:
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
                    if curr.children[j] is not None and j != idxNew:
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
        listItAcum = []

        idx = 0

        for group in listSequences:
            maxVal = 0
            seqPred = ''
            node = trieNo('+')
            currVal = 0

            for it in group:
                insertTree = node.addNode(it)
                currVal = insertTree[0] + 1
                if insertTree[1] == False:

                    listCountTemp.append(currVal)
                    listItTemp.append(it)
                else:
                    currVal = 0
                    listCountTemp.append(currVal)
                    listItTemp.append(it)

            # print("Count occurrences ::::: ", listItTemp , ":::::", listCountTemp)

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
                # print("Compiled of each sequence occurrence :: ", j , " :: ", valMaj)
                listCountAcum.append(valMaj)
                listItAcum.append(j)
                listTempUse.clear()

                valMaj = 0
            listaTempUniq.clear()

            for k in range(len(listCountAcum)):

                if listCountAcum[k] > 0 and listCountAcum[k] >= consLev:
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

            # print("LISTA DE MOTIVOS FINAL :: ", listMotifsFinal)

        return listMotifsFinal, listLocalsFinal, listConsAminFinal, listCountFinal, listPercentFinal, listMotifsNumbers, listOrig

    # Insertion of data in a trie strucuture, generating the motifs and the regions with alterations
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
                        word = word + word1 + word2 + word3
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

        return [listStrFinal, listConsLocalsFinal, listMotifsFinal, listAltCountFinal, listTypeFinal,
                listPercentOccurrFinal]

    def motifBySize(self, listIndxsCons, listStrsCons, minSize, txConser):


        numSp = self.countSpecies()

        lisLocalInit = []
        listStrTemp = []

        listStr, listOcc, listLine, listGroup, listFinalPos, listSize, listLocal = [], [], [], [], [], [], []


        for grp in range(len(listStrsCons)):




            for beg in range(len(listStrsCons[grp][0])):


                localSeq = listIndxsCons[grp][beg]

                for end in range(len(listStrsCons[grp][0]),beg+minSize-1,-1):




                    #largest subsequence identified by tree
                    lgstSeqFound = 0

                    for its in range(len(listStrsCons[grp])):

                        if beg<end:
                            seq = listStrsCons[grp][its][beg:end]
                            lisLocalInit.append(localSeq)
                            listStrTemp.append(seq)


                    node = trieNo("+")



                    #insert subsequences in the trie
                    lTempStr, lTempOcc, lTempLine, lTempLocal = [],[],[],[]
                    for iTree in range(len(listStrTemp)):
                        insTree = node.addNode(listStrTemp[iTree])
                        vl = listStrTemp[iTree]
                        lTempStr.append(vl)
                        cons = insTree[0]/numSp*100
                        lTempOcc.append(cons)
                        lTempLine.append(iTree)
                        lTempLocal.append(lisLocalInit[iTree])

                    listStrTemp.clear()
                    lisLocalInit.clear()


                    #get the percentage occ. of the most common sequence in the group
                    lgsSub = 0
                    for lgs in range(len(lTempOcc)):
                        lenSubs = lTempOcc[lgs]
                        if lenSubs >lgsSub:
                            lgsSub = lenSubs

                    #get the most common sequence
                    spUnit = 100/numSp

                    for tMax in range(len(lTempOcc)):
                        if lgsSub>spUnit and lTempOcc[tMax] == lgsSub:
                            vl = lTempStr[tMax]
                            listStr.append(vl)

                            listOcc.append(lTempOcc[tMax])
                            iniLocal = lTempLocal[tMax]
                            finalLocal = lTempLocal[tMax]+len(vl)
                            size = finalLocal-iniLocal
                            listLocal.append(iniLocal)
                            listFinalPos.append(finalLocal)
                            listGroup.append(grp)
                            listSize.append(size)
                    if len(listStr)>0:
                        self.listMotifs.append(listStr.copy())
                        self.localMotifs.append(listLocal.copy())
                        self.sizesMotifs.append(listSize.copy())
                        self.support.append(listOcc.copy())

                        self.groupMotifs.append(listGroup.copy())
                        self.beginMotifs.append(listLocal.copy())
                        self.endMotifs.append(listFinalPos.copy())

                    listGroup.clear()
                    listLocal.clear()
                    listFinalPos.clear()
                    listOcc.clear()
                    listLocal.clear()
                    listSize.clear()



                    #get the line where the common sequence occurs
                    listTmpLines = []
                    for i in range(len(listStr)):
                        for j in range(len(lTempStr)):
                            if listStr[i] == lTempStr[j]:
                                listTmpLines.append(j)
                    listLine.append(listTmpLines.copy())
                    listTmpLines.clear()
                    self.linesMotifs.append(listLine.copy())
                    listLine.clear()
                    listStr.clear()



    def processBySize(self):
        motifs = self.listMotifs
        supports = self.support
        sizes = self.sizesMotifs
        lines = self.linesMotifs
        group = self.groupMotifs
        begin = self.beginMotifs
        end = self.endMotifs


        listTempMotifs = []
        listTempSupports = []
        listTempSizes = []
        listTempLines = []
        listTempGroup = []
        listTempBegin = []
        listTempEnd = []

        listMotifs = []
        listSupports = []
        listSizes = []
        listLines = []
        listGroup = []
        listBegin = []
        listEnd = []



        listMotifs = []
        lgstMot = 0

        for i in range(len(motifs)):
            for j in range(len(motifs[i])):
                listTempMotifs.append(motifs[i][j])

                listTempSupports.append(supports[i][j])
                listTempSizes.append(sizes[i][j])
                listTempLines.append(lines[i])
                listTempGroup.append(group[i][j])
                listTempBegin.append(begin[i][j])
                listTempEnd.append(end[i][j])


        beginList = []
        begLst = list(set(listTempBegin))


        listLargMot = []
        listBeginL=[]
        for bg in begLst:
            maj = 0
            idx = 0
            for it in range(len(listTempMotifs)):
                if listTempBegin[it]==bg and len(listTempMotifs[it])>maj:
                    maj = len(listTempMotifs[it])
                    idx = it
                    beg = bg
            listLargMot.append(listTempMotifs[idx])



        for i in range(len(listLargMot)):
            for j in range(len(listTempMotifs)):
                if listLargMot[i] == listTempMotifs[j]:
                    listMotifs.append(listTempMotifs[j])
                    listSupports.append(listTempSupports[j])
                    listSizes.append(listTempSizes[j])
                    listBegin.append(listTempBegin[j])
                    listEnd.append(listTempEnd[j])




        return listMotifs, listSupports, listSizes, listBegin, listEnd



    def insertTreeBySize(self, listMers):

        valueMer = 0
        print(listMers)
        listLocals = []
        listLines = []
        listMotifs = []

        for i in range(len(listMers)):
            for j in range(len(listMers[i])):
                node = trieNo("+")
                valueStr = 0
                print(listMers[i][j])
                print("-------------------------------------")
                for k in range(len(listMers[i][j])):

                    for l in range(len(listMers[i][j][k])):
                        valueStr = node.addNode(listMers[i][j][k][l])[0]

                        print(listMers[i][j][k][l], " value ::", valueStr)
                    # print("-------------------------------------")

                # valueMer = valueStr
                # print(valueMer)

    def searchInGroup(self, finalList, group, groupNumber, sizeItem):

        listCounterMers = []
        listaTempMers = []
        listTempLocals = []

        listLocals = []
        listLocalsFinal = []
        listMers = []

        no = trieNo("+")

        for item in group:
            incr = 0
            max = len(item) - sizeItem + 1
            for x in range(max):
                mer = item[x: x + sizeItem]
                listaTempMers.append(mer)
                listTempLocals.append(finalList[groupNumber][incr])
                incr = incr + 1

                # listaTempMers = sorted(set(listaTempMers), key=listaTempMers.index)

        x = -1
        for mer in listaTempMers:
            listCounterMers.append(no.addNode(mer)[0] + 1)
            listMers.append(mer)
            localMer = listTempLocals[x + 1]
            listLocals.append(listTempLocals[x + 1])
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


root = trieNo('+')
root.fileName = "_testeAligned1.fas"

especiesList = root.readSpecies()

seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 90

# root.conservPerc = 70


locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

root.motifBySize(indices, conservados, 2, 80)
root.processBySize()
'''


for i in range(len(root.listMotifs)):

    print('Motif: ',root.listMotifs[i])

    print('Support: ',root.support[i])
    print('Size: ',root.sizesMotifs[i])
    print('Lines: ',root.linesMotifs[i])
    print('Group: ',root.groupMotifs[i])
    print('Begin: ',root.beginMotifs[i])
    print('Final: ',root.endMotifs[i])
    print('------------------------------------------------------')


'''