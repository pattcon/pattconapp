
from .model import trieNo


class conserv_alter():
    nomeArq = '001.fas'
    percCont = 0
    percCons = 0
    lstSpAlter = []
    lstSpNoAlter = []



    def executeBySize(self, motsize, txCons, txCont):



        root = trieNo('+')
        root.fileName = self.nomeArq
        root.conservPerc = txCons
        root.contentPerc = txCont

        especiesList = root.readSpecies()
        seqAminos = root.readSequences()

        qtSpecies = root.countSpecies()
        porceCon = root.calcConservPerc()
        locaisConservados = root.conservLocals(porceCon)
        indices = root.groupConsLocals(locaisConservados)

        conservados = root.genConserved(indices, seqAminos)

        listMotifs = root.motifBySize(indices,conservados,motsize, txCons)[0]
        listLocals = root.motifBySize(indices,conservados,motsize, txCons)[1]
        listSpec = root.motifBySize(indices,conservados,motsize, txCons)[2]
        listSupport = root.motifBySize(indices,conservados,motsize, txCons)[3]
        listEndLocal = []
        listAlt = root.motifBySize(indices, conservados, motsize, txCons)[4]
        listNoAlt = root.motifBySize(indices, conservados, motsize, txCons)[5]

        listStrAlt = root.motifBySize(indices, conservados, motsize, txCons)[6]
        listStrNoAlt = root.motifBySize(indices, conservados, motsize, txCons)[7]
        listStrNorm = root.motifBySize(indices, conservados, motsize, txCons)[8]
        listNumAlt = root.motifBySize(indices, conservados, motsize, txCons)[9]



        listSpecN = []
        lsTemp =[]


        for m in range(len(listSpec)):
            for n in range(len(listSpec[m])):
                idx = listSpec[m][n]
                sp = especiesList[idx]
                lsTemp.append(sp)
            listSpecN.append(lsTemp.copy())
            lsTemp.clear()

        listSpecNames = []
        for p in range(len(listSpecN)):
            listSpecNames.append(listSpecN[p])

        listSizes = []
        listLocalsFinal = []
        for k in range(len(listMotifs)):
            listSizes.append(len(listMotifs[k]))
            listEndLocal.append(listLocals[k]+len(listMotifs[k]))
            listLocalsFinal.append(listLocals[k]+1)








        return listMotifs, listLocalsFinal, listSpec, listSupport, listSpecNames, listSizes, listEndLocal, listAlt, listNoAlt, listStrAlt, listStrNoAlt, especiesList, listStrNorm, listNumAlt

