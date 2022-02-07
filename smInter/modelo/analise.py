
from .model import trieNo


class conserv_alter():
    nomeArq = '001.fas'
    percCont = 0
    percCons = 0

    def executar(self):
        root = trieNo('+')
        root.fileName = self.nomeArq

        root.conservPerc = self.percCons
        root.contentPerc = self.percCont

        especiesList = root.readSpecies()
        seqAminos = root.readSequences()

        qtdeEspecies = root.countSpecies()
        porceCon = root.calcConservPerc()
        locaisConservados = root.conservLocals(porceCon)
        indices = root.groupConsLocals(locaisConservados)



        conservados = root.genConserved(indices, seqAminos)
        listResult = root.genMotifs(conservados, indices)
        #predominantes = root.geraPred(conservados)

        listConsOrig = listResult[6]
        lista = root.processSequences(listResult[0], listResult[2], listResult[1],listResult[4])
        motivos =lista[2]
        locais = lista[1]

        listaStrFim = lista[0]
        listaContaAlt = lista[3]
        listaTipo = lista[4]

        #numerMotivo = listResult[5]
        numerMotivo = []

        qtde = len(motivos)
        for i in range(0, qtde):
            numerMotivo.append(i)



        inicioLocal = []
        fimLocais = []

        for k in range(len(locais)):
            size = len(locais[k])
            inicioLocal.append(locais[k][0])
            fimLocais.append(locais[k][size-1])

        inFimJuntos = []
        for l in range(len(inicioLocal)):
            if l>1:
                inFimJuntos.append([numerMotivo[l-1], inicioLocal[l], fimLocais[l]])

        listaEsp = especiesList

        listaAltStr = []
        listaAltEsp = []
        listaAltQtde = []

        listaConsStr = []
        listaConsEsp = []


        porcentOcorr = listResult[4]

        print("Motivos :: ", motivos)
        print("Numeros :: ", numerMotivo)


        return numerMotivo, locais, motivos, inicioLocal, fimLocais, inFimJuntos,listaStrFim, listaContaAlt, listaTipo, listaEsp,listaAltStr,listaAltEsp,listaAltQtde,listaConsStr,listaConsEsp, listConsOrig, porcentOcorr, porceCon



    def executeBySize(self, motsize, txCons, txCont):
        root = trieNo('+')
        root.fileName = self.nomeArq
        root.conservPerc = txCons
        root.contentPerc = txCont

        especiesList = root.readSpecies()
        seqAminos = root.readSequences()

        qtdeEspecies = root.countSpecies()
        porceCon = root.calcConservPerc()
        locaisConservados = root.conservLocals(porceCon)
        indices = root.groupConsLocals(locaisConservados)

        conservados = root.genConserved(indices, seqAminos)

        listMotifsRes = root.motifBySize(conservados, indices,motsize, txCons, txCont)
        #listLocalsRes = root.motifBySize(conservados, indices, motsize, txCons)[1]
        #listSuportsRes = root.motifBySize(conservados, indices, motsize, txCons)[2]


        """
        tamMax = 0

        for p in range(motsize, len(listMotifsRes)):
            for q in range(len(listMotifsRes[p])):
                tamCurr = len(listMotifsRes[p][q])
                if tamCurr>tamMax:
                    tamMax = tamCurr
        """



        return listMotifsRes

