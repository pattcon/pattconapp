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


    def motifBySize(self, listIndxsCons, listStrsCons, minSize, txConser):


        numSp = self.countSpecies()

        lisLocalInit = []
        listStrTemp = []

        listStr, listOcc, listLine, listGroup, listFinalPos, listSize, listLocal = [], [], [], [], [], [], []
        listAllMot, listAllLoc, listAllOcc, listMotFinal, listLocFinal, listOccFinal= [],[],[],[],[],[]
        listAllSeqs, listAllNumAlter, listAllStrForm, lsAllStrNorm=[],[],[],[]
        lstSeqs = []

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
                    lTempStr, lTempOcc, lTempLocal = [],[],[]
                    for iTree in range(len(listStrTemp)):
                        insTree = node.addNode(listStrTemp[iTree])
                        vl = listStrTemp[iTree]


                        cons = insTree[0] / numSp * 100

                        #if cons>=txConser:
                        lTempStr.append(vl)
                        lTempOcc.append(cons)
                        lTempLocal.append(lisLocalInit[iTree])



                    #get the occurrence's count of the most frequent subsequence
                    lgsSub=0
                    for iit in range(len(lTempOcc)):
                        if lTempOcc[iit]>=lgsSub:
                            lgsSub =  lTempOcc[iit]


                    #get the positions of the most common substring


                    for iis in range(len(lTempStr)):

                        if lgsSub>= txConser and lTempOcc[iis] == lgsSub:


                            lstMosStr, lstMosOcc, lstMosLoc = [], [], []
                            vl = lTempStr[iis]


                            loc = lisLocalInit[iis]


                            for iit in range(len(lTempStr)):
                                if vl == lTempStr[iit]:
                                    lstMosOcc.append(iit)
                                lstSeqs.append(lTempStr[iit])
                            lstMosStr.append(vl)
                            lstMosLoc.append(loc)

                            listAllMot.append(lstMosStr)
                            listAllLoc.append(lstMosLoc)
                            listAllOcc.append(lstMosOcc)
                            listAllSeqs.append(lstSeqs)



                            lsTemp, lsAltTmp, lsTemp2,lsTemp3, lsTempAlt2 =[], [],[],[], []
                            for ls in range(len(lstSeqs)):
                                alter = 0
                                strng = ''
                                strng2 = ''
                                for lt in range(len(lstSeqs[ls])):
                                    for mot in range(len(lstMosStr)):
                                        strng2 = strng2+lstSeqs[ls][lt]
                                        if lstMosStr[mot][lt] == lstSeqs[ls][lt]:
                                            strng = strng+lstMosStr[mot][lt]
                                        else:
                                            strng = strng + '['+lstSeqs[ls][lt]+']'
                                            alter = alter+1
                                lsTemp.append(strng)

                                lsAltTmp.append(alter)
                                lsTemp3.append(strng2)

                            lsAllStrNorm.append(lsTemp3.copy())
                            listAllStrForm.append(lsTemp.copy())
                            listAllNumAlter.append(lsAltTmp.copy())
                            lsAltTmp.clear()
                            lsTemp.clear()
                            lsTemp3.clear()

                    lstSeqs.clear()
                    listStrTemp.clear()
                    lisLocalInit.clear()

                    lTempStr.clear()
                    lTempOcc.clear()



        lstTempMot, lstTempLoc, lstTempOcc, lstTempSupp = [], [], [], []
        lstAllStrForm, lstAllNumAlter, lstAllStrNorm = [],[],[]


        for i in range(len(listAllMot)):
            for j in range(len(listAllMot[i])):
                lstTempMot.append(listAllMot[i][j])
                lstTempLoc.append(listAllLoc[i][j])
                lstTempOcc.append(listAllOcc[i])
                lstTempSupp.append(len(lstTempOcc[i])/self.specNumb*100)
                lstAllStrForm.append(listAllStrForm[i])
                lstAllNumAlter.append(listAllNumAlter[i])
                lstAllStrNorm.append(lsAllStrNorm[i])


        lstLocAlt, lstLocNoAlt, lstSpecAlt, lstSpecNoAlt = [], [], [], []
        lTmpSpecAlt, lTmpSpecNoAlt, lsTmpAlt, lsTmpNoAlt = [],[],[],[]
        lTempStrAlt, lTempStrNoAlt, lstStrAlt, lstStrNoAlt = [],[],[],[]


        #identify alterations locals in the sequences
        for ab in range(len(lstTempMot)):
            for ac in range(len(lstAllNumAlter[ab])):

                if lstAllNumAlter[ab][ac] > 0:
                    lTmpSpecAlt.append(ac)
                    lTmpSpecNoAlt.append(0)
                    lsTmpAlt.append(lstAllNumAlter[ab][ac])
                    lsTmpNoAlt.append(0)
                    lTempStrAlt.append(listAllStrForm[ab][ac])
                else:

                    lTmpSpecAlt.append(0)
                    lTmpSpecNoAlt.append(ac)
                    lsTmpAlt.append(0)
                    lsTmpNoAlt.append(lstAllNumAlter[ab][ac])
                    lTempStrNoAlt.append(listAllStrForm[ab][ac])

            lstLocAlt.append(lsTmpAlt.copy())
            lstLocNoAlt.append(lsTmpNoAlt.copy())
            lstSpecAlt.append(lTmpSpecAlt.copy())
            lstSpecNoAlt.append(lTmpSpecNoAlt.copy())
            lstStrAlt.append(lTempStrAlt.copy())
            lstStrNoAlt.append(lTempStrNoAlt.copy())

            lTmpSpecAlt.clear()
            lTmpSpecNoAlt.clear()
            lsTmpAlt.clear()
            lsTmpNoAlt.clear()
            lTempStrAlt.clear()
            lTempStrNoAlt.clear()



        listSpAltF, listSpNoAltF = [], []
        lsTmpAlt, lsTmpNoAlt = [],[]

        #get the name of species
        for spA in range(len(lstSpecAlt)):
            for spB in range(len(lstSpecAlt[spA])):
                if lstSpecAlt[spA][spB]>0:
                    lsTmpAlt.append(self.readSpecies()[spB])
                else:
                    lsTmpNoAlt.append(self.readSpecies()[spB])
            listSpAltF.append(lsTmpAlt.copy())
            listSpNoAltF.append(lsTmpNoAlt.copy())
            lsTmpAlt.clear()
            lsTmpNoAlt.clear()






        idx = 0
        curr = 0


        #eliminate the smaller motifs that is in the other largest motif
        lstLocFinal, lstMotFinal, lstOccFinal, lstSuppFinal, lstMotifEnd =[], [], [], [], []
        listSpAltFinal, listSpNoAltFinal, listAllMotFinal =[],[],[]
        listStrAltFinal, listStrNoAltFinal, listAllStrNorm, listNumAltFinal = [],[], [],[]


        while idx<len(lstTempMot):
            lstLocFinal.append(lstTempLoc[idx])
            lstMotFinal.append(lstTempMot[idx])
            lstOccFinal.append(lstTempOcc[idx])
            lstSuppFinal.append(lstTempSupp[idx])
            lstMotifEnd.append(lstTempLoc[idx]+len(lstTempMot[idx]))

            listSpAltFinal.append(listSpAltF[idx])
            listSpNoAltFinal.append(listSpNoAltF[idx])
            listStrAltFinal.append(lstStrAlt[idx])
            listStrNoAltFinal.append(lstStrNoAlt[idx])
            listAllStrForm.append(lstAllStrForm[idx])
            listAllStrNorm.append(lstAllStrNorm[idx])
            listNumAltFinal.append(lstLocAlt[idx])


            while lstTempMot[idx] in lstTempMot[curr] and lstTempLoc[idx] <= lstTempLoc[curr]+len(lstTempMot[curr]):


                idx = idx + 1
                size = len(lstTempMot)
                if idx <size:


                    continue
                else:
                    break
            curr = idx

        return lstMotFinal, lstLocFinal, lstOccFinal, lstSuppFinal, listSpAltFinal, listSpNoAltFinal, listStrAltFinal, listStrNoAltFinal, listAllStrNorm,listNumAltFinal
    
'''
root = trieNo('+')
root.fileName = "testsequence.fas"

especiesList = root.readSpecies()

seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 90

# root.conservPerc = 70


locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

root.motifBySize(indices, conservados, 4, 50)

'''