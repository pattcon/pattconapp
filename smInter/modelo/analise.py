
from .model import trieNo


class conserv_alter():
    nomeArq = '001.fas'
    por = 60

    def executar(self):
        root = trieNo('+')
        root.nomeDoArquivo = self.nomeArq
        root.porcenCon = self.por

        especiesList = root.leEspecies()
        seqAminos = root.leSequencias()

        qtdeEspecies = root.contaEspecies()
        porceCon = root.porcConservacao(qtdeEspecies)
        locaisConservados = root.locaisConserv(porceCon)
        indices = root.agrupaSitesCons(locaisConservados)

        conservados = root.geraConserv(indices, seqAminos)
        listResult = root.geraFim(conservados, indices)
        #predominantes = root.geraPred(conservados)

        listConsOrig = listResult[6]
        lista = root.geraListaFim(listResult[0], listResult[2], listResult[1])
        motivos =lista[2]
        locais = lista[1]

        listaStrFim = lista[0]
        listaContaAlt = lista[3]
        listaTipo = lista[4]

        #numerMotivo = listResult[5]
        numerMotivo = []

        qtde = len(motivos)
        for i in range(1, qtde):
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



        return numerMotivo, locais, motivos, inicioLocal, fimLocais, inFimJuntos,listaStrFim, listaContaAlt, listaTipo, listaEsp,listaAltStr,listaAltEsp,listaAltQtde,listaConsStr,listaConsEsp, listConsOrig, porcentOcorr



    def executeBySize(self, motsize, txCons, txCont):
        root = trieNo('+')
        root.nomeDoArquivo = self.nomeArq
        root.porcenCon = self.por

        especiesList = root.leEspecies()
        seqAminos = root.leSequencias()

        qtdeEspecies = root.contaEspecies()
        porceCon = root.porcConservacao(qtdeEspecies)
        locaisConservados = root.locaisConserv(porceCon)
        indices = root.agrupaSitesCons(locaisConservados)

        conservados = root.geraConserv(indices, seqAminos)

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

