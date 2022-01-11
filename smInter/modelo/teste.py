from Bio import SeqIO
import re


class trieNo(object):

    def __init__(self, chtr):
        self.chr = chtr
        self.seq = []
        self.fimpalavra = False
        self.contaSeq = 0
        self.contaAmino = 1
        self.porcenCon = 60
        self.tipo = 1
        self.ranking = []
        self.contaDiferentes = 0
        self.wordsize = 0
        self.nomeDoArquivo = ''
        self.formatoArquivo = 'fasta'
        self.numTotalResiduos = 0
        self.numEspecies = 0

        self.children = [None] * 31

    # lê a sequencia de aminos no arquivo Fasta
    def leSequencias(self):

        listaSequencias = []

        for i in SeqIO.parse(self.nomeDoArquivo, self.formatoArquivo):
            listaSequencias.append(list(i.seq))

        return listaSequencias

    def leQtdeResiduos(self, listaSeq):
        seq = self.leSequencias()
        tam = len(seq[0])
        self.numTotalResiduos = tam
        return seq

    # formata e adiciona os nomes das espécies à uma lista
    def leEspecies(self):
        listaEspecies = []

        for i in SeqIO.parse(self.nomeDoArquivo, self.formatoArquivo):
            listaTemp = []
            listaTemp.append(re.split(r'_', i.id))

            if len(listaTemp[0][0]) < 3:
                item = listaTemp[0][0] + ' ' + listaTemp[0][1]
            else:
                item = listaTemp[0][0]

            listaEspecies.append(item)

        return listaEspecies

    def contaEspecies(self):
        num = self.numEspecies = len(self.leEspecies())
        self.ranking = [0] * num
        return num

    # porcentagem a se considerar como site conservado - 60%
    def porcConservacao(self, numeroEspecies):
        return int(numeroEspecies * self.porcenCon / 100)

    # elimina os GAP's das sequências e verifica os aminoácidos remanescentes se atendem à quantidade mínima exigida para serem
    # considerados como uma região conservada.
    def locaisConserv(self, porConserv):
        listaSequencias = self.leSequencias()

        # for i in SeqIO.parse(self.nomeDoArquivo, self.formatoArquivo):
        # listaSequencias.append(i.seq)

        listaEspecies = self.leEspecies()

        locaisConservados = []
        sitesConservados = []
        especiesConserv = []

        for j in range(len(listaSequencias[0])):
            cont = 0
            colunas = []
            for k in range(len(listaSequencias)):
                if (listaSequencias[k][j] != '-'):
                    colunas.append(listaSequencias[k][j])

            if len(colunas) >= self.porcConservacao(self.numEspecies):
                locaisConservados.append(j)

            colunas.clear()

        return locaisConservados

    # agrupa os sites conservados consecutivos
    def agrupaSitesCons(self, listCon):
        listaTemp = []
        listaFim = []
        m = 0
        incr = 0
        tam = len(listCon)

        while incr < tam:
            chave = incr
            while listCon[incr] - listCon[chave] == m:
                listaTemp.append(listCon[incr])
                m = m + 1
                incr = incr + 1
                if incr == tam:
                    break
            listaFim.append(listaTemp.copy())
            listaTemp.clear()
            m = 0
        return listaFim

        # OK - indices

    # pega os aminoacidos conservados por linha dentro da taxa de conservação estabelecida:
    def geraConserv(self, listaCon, listaSeq):
        listTemp = []
        listaAdd = []
        listaLocal = []

        for item in listaCon:
            min = item[0]
            tam = len(item) - 1
            max = item[tam]
            for p in range(len(listaSeq)):
                listTemp.append(listaSeq[p][min:max + 1])

            listaAdd.append(listTemp.copy())

            listTemp.clear()
        ##
        # converte cada linha em uma palavra pra ser usada na Trie
        lTemp = []
        listaStr = []
        separator = ''
        for item in listaAdd:
            for j in item:
                strAm = separator.join(j)
                lTemp.append(strAm)
            listaStr.append(lTemp.copy())
            lTemp.clear()

        return listaStr

    # ok - conservados - linhas

    # gera as sequências mais encontradas em cada grupo de sites
    def motifBySize(self, listaStr, listaFim, motSize):
        listQtMers = []
        listaMotivosFim = []
        listaTempQuant = []
        listaLocaisFim = []
        listaAmConservFim = []
        listaQuantFim = []
        listaPorcentFim = []
        listMotifsNumbers = []
        listOrig = []

        idx = 0

        for grupo in listaStr:
            seqPred = ''
            no = trieNo('+')
            maior = 0
            valor = 0

            motif = ""

            for item in grupo:
                tam = len(item)
                print("Sequence ::" , item)


                listAdd = []
                for x in range(len(item)):
                    part = item[x:x + motSize]

                    if len(part) == motSize:
                        listAdd.append(part)

                    myList = set(listAdd)
                    myList = list(myList)

                print("Mers List :::: ", myList)



                for mer in myList:
                    print("MER :: ", mer)
                    valor = no.adicionaNo(mer) + 1
                    listaTempQuant.append(valor)
                    listQtMers.append(mer)
                    semMotivo = 0

                    if valor > maior:
                        maior = valor
                        motif = mer


                if maior > 1:
                    listaMotivosFim.append(motif)
                    listaLocaisFim.append(listaFim[idx])
                    listaAmConservFim.append(grupo)
                    listaQuantFim.append(listaTempQuant[0])
                    listaPorcentFim.append(round(maior / self.numEspecies * 100))
                    listOrig.append(grupo)
            print("Motif ::: ", motif)

            for p in range(len(listaMotivosFim)):
                listMotifsNumbers.append(p)

            idx = idx + 1

            listaTempQuant.clear()

        return listaMotivosFim, listaLocaisFim, listaAmConservFim, listaQuantFim, listaPorcentFim, listMotifsNumbers, listOrig







    def adicionaNo(self, key):

        curr = self
        exist = True
        contaGaps = 0

        for i in range(len(key)):
            if key[i] == '-':
                contaGaps = contaGaps + 1

        if contaGaps - 1 != len(key):

            for i in range(len(key)):
                index = ord(key[i]) - ord('A')

                if index == -20:
                    index = 29

                if curr.children[index] is None:
                    curr.children[index] = trieNo(index)
                    exist = False

                curr = curr.children[index]

                if exist == True:
                    curr.contaSeq = curr.contaSeq + 1

        return curr.contaSeq




    def geraListaFim(self, listaMotivos, listaAminos, listaLocCons):
        listaStrTemp = []
        listaStrFim = []
        conservada = 1
        listaTipoTemp = []
        listaTipoFim = []
        listaLocaisConsFim = []
        listaContaAltTemp = []
        listaContaAltFim = []
        listaMotivoFim = []

        for i in range(len(listaAminos)):
            listaLocaisConsFim.append(listaLocCons[i])
            listaMotivoFim.append(listaMotivos[i])

            for j in range(len(listaAminos[i])):
                word = ""
                # print("Motivo em teste :::", listaMotivos[i], "Sequencia em teste ::", listaAminos[i][j])
                conservada = 1

                # listaLocaisConsFim.append(listaLocCons.copy())

                alteracoes = 0
                for k in range(len(listaAminos[i][j])):
                    if listaAminos[i][j][k] == listaMotivos[i][k]:
                        word = word + listaAminos[i][j][k]

                    else:
                        word = word + "[" + listaAminos[i][j][k] + "]"
                        conservada = 0
                        alteracoes = alteracoes + 1

                listaContaAltTemp.append(alteracoes)

                listaStrTemp.append(word)
                listaTipoTemp.append(conservada)

            listaTipoFim.append(listaTipoTemp.copy())
            listaContaAltFim.append(listaContaAltTemp.copy())
            listaStrFim.append(listaStrTemp.copy())
            listaStrTemp.clear()
            listaContaAltTemp.clear()
            listaTipoTemp.clear()

        return [listaStrFim, listaLocaisConsFim, listaMotivoFim, listaContaAltFim, listaTipoFim]


root = trieNo('+')
root.nomeDoArquivo = "001.fas"
root.porcenCon = 60

especiesList = root.leEspecies()
seqAminos = root.leSequencias()

qtdeEspecies = root.contaEspecies()
porceCon = root.porcConservacao(qtdeEspecies)
locaisConservados = root.locaisConserv(porceCon)
indices = root.agrupaSitesCons(locaisConservados)

conservados = root.geraConserv(indices, seqAminos)

listResult = root.motifBySize(conservados, indices, 2)

#lista = root.geraListaFim(listResult[0], listResult[2], listResult[1])

# print("Motivos: ", listResult[4], "Total: ", len(listResult[4]))
