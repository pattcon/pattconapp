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
        return int(self.specNumb * self.contentPerc / 100)

    def calcConservPerc(self):
        return int(self.specNumb * self.conservPerc / 100)



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











    # Generate the most finding sequences in each group
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

        idx = 0

        for group in listSequences:
            maxVal = 0
            seqPred = ''
            node = trieNo('+')
            currVal = 0

            for it in group:

                if node.addNode(it)[1] == False:
                    currVal = node.addNode(it)[0] + 1
                    listCountTemp.append(currVal)
                else:
                    currVal = 0
                    listCountTemp.append(currVal)

                if currVal > maxVal:
                    maxVal = currVal

            for k in range(len(group)):
                if maxVal >= consLev:
                    if listCountTemp[k] >= maxVal:
                        listMotifsFinal.append(group[k])
                        listLocalsFinal.append(listConserved[idx])
                        listConsAminFinal.append(group)
                        listCountFinal.append(listCountTemp[k])
                        #listPercentFinal.append(round(maxVal / self.specNumb * 100))
                        listOrig.append(group)

            for p in range(len(listMotifsFinal)):
                listMotifsNumbers.append(p)

            idx = idx + 1

            listCountTemp.clear()

        return listMotifsFinal, listLocalsFinal, listConsAminFinal, listCountFinal, listPercentFinal, listMotifsNumbers, listOrig














    #Insertion of data in a trie strucuture, generating the motifs and the regions with alterations
    def processSequences(self, listMotifs, listConservedAmin, listConservLocals):
        listStrTemp = []
        listStrFinal = []
        conserved = 1
        listTypTemp = []
        listTypeFinal = []
        listConsLocalsFinal = []
        listAltCountTemp = []
        listAltCountFinal = []
        listMotifsFinal = []




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




        return [listStrFinal, listConsLocalsFinal, listMotifsFinal, listAltCountFinal, listTypeFinal]


