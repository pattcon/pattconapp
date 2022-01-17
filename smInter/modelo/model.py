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
        self.numTotalResiduos= 0
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

        #for i in SeqIO.parse(self.nomeDoArquivo, self.formatoArquivo):
            #listaSequencias.append(i.seq)

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
    def geraFim(self, listaStr, listaFim):

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
            maior = 0
            seqPred = ''
            no = trieNo('+')
            valor = 0


            for item in grupo:
                valor = no.adicionaNo(item)+1
                listaTempQuant.append(valor)
                semMotivo = 0


                if valor > maior:
                    maior = valor

            for k in range(len(grupo)):
                if maior >1:
                    if listaTempQuant[k] >= maior:
                        listaMotivosFim.append(grupo[k])
                        listaLocaisFim.append(listaFim[idx])
                        listaAmConservFim.append(grupo)
                        listaQuantFim.append(listaTempQuant[k])
                        listaPorcentFim.append(round(maior / self.numEspecies * 100))
                        listOrig.append(grupo)
            """
                else:
                    semMotivo = 1

            if semMotivo == 1:
                for l in range(len(item)):
                    seqPred = seqPred + '-'
                listaMotivosFim.append(seqPred)
                listaLocaisFim.append(listaFim[idx])
                listaAmConservFim.append(grupo)
                listaQuantFim.append(0)
                listaPorcentFim.append(0)
                listOrig.append(grupo)
                seqPred = ''
                semMotivo = 0

            """


            for p in range(len(listaMotivosFim)):
                listMotifsNumbers.append(p)


            idx = idx + 1

            listaTempQuant.clear()


        return listaMotivosFim, listaLocaisFim, listaAmConservFim, listaQuantFim, listaPorcentFim, listMotifsNumbers, listOrig

    def motifBySize(self, listaStr, finalList, minSize, txConser, txContent):
        minConser = txConser / 100 * self.contaEspecies()
        minContent = txConser / 100 * self.contaEspecies()


        listTempMotif = []
        listQtMers = []
        listaMotivosFim = []

        listaLocaisFim = []

        listaPorcentFim = []
        listMotifsNumbers = []

        listLocals = []
        listLocalsMotifs = []
        listLocalsMotifsFim = []
        idxMotif = 0

        listFinalLocals = []
        listCFinalMotifs = []
        listaFinalCont = []
        groupNumber = 0
        result = []

        for grupo in listaStr:




            groupSize = len(grupo[0])

            for gSz in range(0, groupSize-minSize+1):
                searchMotif = self.searchInGroup(finalList, grupo, groupNumber, groupSize-gSz)
                if searchMotif[0] >= minConser:
                    #value = "Motif :", searchMotif[1], "- Local :: ", searchMotif[2],
                    value = "Motif: {} - Local: {}".format(searchMotif[1], searchMotif[2])
                    result.append(value)
                    result.append("\n")
                    print(result)
                    break

            groupNumber = groupNumber+1

            occurrList = searchMotif[0]
            motifsList = searchMotif[1]
            localsList  = searchMotif[2]






        """
                    listLocals.append(listSpeciesTemp.copy())

                    listTempLocals = []
                    listCounterMers = []
                    listaTempMerOccur = []

                    incr = 0



                    for mer in myList:

                        listCounterMers.append(no.adicionaNo(mer) + 1)
                        listaTempMerOccur.append(mer)
                        listTempLocals.append(listaFim[grupo][incr])
                        incr = incr + 1

                listaLocaisFim.extend(listTempLocals)

                maxOccur = 0
                for r in range(len(listCounterMers)):
                    if listCounterMers[r]>minContent:

                        listFinalLocals.append(listTempLocals[r])
                        listCFinalMotifs.append(listaTempMerOccur[r])
                        listaFinalCont.append(listCounterMers[r])

                        maxOccur = r
                listaTempMerOccur.clear()
                listTempLocals.clear()
                listCounterMers.clear()

                print("motivos :: ", listCFinalMotifs)
                print("-------------------------------------------------------------------------")
                #print("contagem :: ", listaFinalCont)
                #print("local :: ", listFinalLocals)








                


                
                listMersTemp = sorted(set(listaTempQuant), key=listaTempQuant.index)
                
                

                
                pos = 0
                listCountMers = []
                listLocalsMotifsTemp = []
                for p in range(len(listMersTemp)):
                    for q in range(len(listaTempQuant)):
                        if listMersTemp[p] == listaTempQuant[q]:
                            pos = q
                    listCountMers.append(listCounterMers[pos])
                    listLocalsMotifsTemp.append(listaLocaisFim[pos])
                listaLocaisFim.clear()

                listSuportTemp = []

                for r in range(len(listCountMers)):
                    if listCountMers[r] >= minConser:
                        listTempMotif.append(listMersTemp[r])
                        listSuportTemp.append(round(listCountMers[r] / self.numEspecies * 100))
                        listLocalsMotifs.append(listLocalsMotifsTemp[r])

                if len(listTempMotif) > 0:
                    listaPorcentFim.append(listSuportTemp.copy())
                    listLocalsMotifsFim.append(listLocalsMotifs.copy())
                    listaMotivosFim.append(listTempMotif.copy())

                listLocalsMotifs.clear()
                listLocalsMotifsTemp.clear()
                listaTempQuant.clear()
                listCountMers.clear()
                listMersTemp.clear()
                listSuportTemp.clear()
                listCounterMers.clear()

                listTempMotif.clear()

            listQtMers.clear()
            listaTempQuant.clear()

            for p in range(len(listaMotivosFim)):
                listMotifsNumbers.append(p)

            idx = idx + 1
            """
        return motifsList, occurrList, localsList, result



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
            listCounterMers.append(no.adicionaNo(mer) + 1)
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





    def adicionaNo(self, key):

        curr = self
        exist = True
        contaGaps = 0

        for i in range(len(key)):
            if key[i] == '-':
                contaGaps = contaGaps + 1

        if contaGaps-1 != len(key):

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


    def busca2(self, key):

        curr = self
        saida = []
        numAlteracoes = 0

        for i in range(len(key)):

            if key[i] == '-':
                index = ord('^') - ord('A')
            else:
                index = ord(key[i]) - ord('A')

            if curr.children[index] is None:
                numAlteracoes = numAlteracoes + 1
                saida.append('{}{}{}'.format('[', key[i], ']'))

                for j in range(31):
                    if curr.children[j] is not None:
                        index = j
                        break

            else:
                saida.append(key[i])
            curr = curr.children[index]


        return [numAlteracoes, saida, i]



    def busca(self, listaMotivos, listaSequencias):


        temp = []
        numAlteracoes = 0
        saidaFim = []

        for j in range(len(listaMotivos)):
            curr = self
            motivo = listaMotivos[j]
            curr.adicionaNo(motivo)




            for k in range(len(listaSequencias)):

                for l in range(len(listaSequencias[k])):
                    key = listaSequencias[k][l]
                    endWord = ""

                    for i in range(len(key)):

                        if key[i] == '-':
                            index = ord('^') - ord('A')
                        else:
                            index = ord(key[i]) - ord('A')


                        if curr is not None:


                            if curr.children[index] is None:
                                numAlteracoes = numAlteracoes + 1
                                #temp.append('{}{}{}'.format('[', key[i], ']'))
                                st = '['+ key[i]+']'
                                endWord = endWord+st

                                for j in range(31):
                                    if curr.children[j] is not None:
                                        index = j
                                        break

                            else:

                                #temp.append(key[i])
                                endWord = endWord+key[i]

                            curr = curr.children[index]
                    temp.append(endWord)



                    saidaFim.append(temp.copy())
                    temp.clear()





        return [numAlteracoes, saidaFim, i]



    # gera ranking de espécies mais conservadas
    def insereAlte(self, indice):
        val = self.ranking[indice]
        val = val + 1
        self.ranking[indice] = val

    def geraRankOrd(self, listaEsp):

        vetEspOrd = [i for i in listaEsp]
        vetEspDes = [i for i in listaEsp]

        vetOrd = [int(i) for i in self.ranking]
        vetDes = [int(i) for i in self.ranking]

        for i in range(len(vetOrd)):

            min = i

            for j in range(i + 1, len(vetOrd)):
                if vetOrd[min] > vetOrd[j]:
                    min = j

            temp = vetOrd[i]
            vetOrd[i] = vetOrd[min]
            vetOrd[min] = temp

            temp2 = vetEspOrd[i]
            vetEspOrd[i] = vetEspOrd[min]
            vetEspOrd[min] = temp2

        return [vetDes, vetOrd, vetEspDes, vetEspOrd]


    def compara(self, listPred, listAminos, listEsp, listAgrup):

        listaAlterada = []
        listaConservada = []
        listaAltTemp = []
        listaConsTemp = []
        espConsTemp = []
        espAltTemp = []
        espConsFinal = []
        espAltFinal = []
        seqPred = []
        listNumAlterac = []
        numAltera = 0

        no = trieNo('+')
        no.adicionaNo(listPred)
        seqPred.append(listPred)

        for j in range(len(listAminos)):

            item = listAminos[j]
            busca = no.busca(item)
            if busca[0] > 0:
                listaAltTemp.append(busca[1])
                espAltTemp.append(listEsp[j])
                listNumAlterac.append(busca[0])
                #print("Especie com alterações:::", listEsp[j])
                #print(" Qtde:  ", busca[0])



            else:
                listaConsTemp.append(busca[1])
                espConsTemp.append(listEsp[j])

        # print("Espcies conservadas ::", espConsTemp)
        # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        listaAlterada.append(listaAltTemp.copy())
        listaConservada.append(listaConsTemp.copy())

        espConsFinal.append(espConsTemp.copy())
        espAltFinal.append(espAltTemp.copy())

        espConsTemp.clear()
        espAltTemp.clear()
        listaAltTemp.clear()
        listaConsTemp.clear()

        return [listaAlterada, listaConservada, espConsFinal, espAltFinal, listAgrup, seqPred, listNumAlterac]


    def compara2(self, listPred, listAminos, listEsp, listAgrup):

        listaAlterada = []
        listaConservada = []
        listaAltTemp = []
        listaConsTemp = []
        espConsTemp = []
        espAltTemp = []
        espConsFinal = []
        espAltFinal = []
        seqPred = []
        listNumAlterac = []
        numAltera = 0
        listaMotivos = []

        no = trieNo('+')
        no.adicionaNo(listPred)
        seqPred.append(listPred)


        for i in range(len(listPred)):
            for j in range(len(listAminos)):

                item = listAminos[j]
                busca = no.busca(item)
                if busca[0] > 0:
                    listaAltTemp.append(busca[1])
                    espAltTemp.append(listEsp[j])
                    listNumAlterac.append(busca[0])
                    #print("Especie com alterações:::", listEsp[j])
                    #print(" Qtde:  ", busca[0])



                else:
                    listaConsTemp.append(busca[1])
                    espConsTemp.append(listEsp[j])

                listaMotivos.append(listPred[i])

        # print("Espcies conservadas ::", espConsTemp)
        # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        listaAlterada.append(listaAltTemp.copy())
        listaConservada.append(listaConsTemp.copy())

        espConsFinal.append(espConsTemp.copy())
        espAltFinal.append(espAltTemp.copy())

        espConsTemp.clear()
        espAltTemp.clear()
        listaAltTemp.clear()
        listaConsTemp.clear()

        return [listaAlterada, listaConservada, espConsFinal, espAltFinal, listAgrup, seqPred, listNumAlterac]


    def ordena(self, listaAl, listaEsp, listNum):
        if len(listaAl) > 1:
            mid = len(listaAl) // 2

            lefthalf = listaAl[:mid]
            righthalf = listaAl[mid:]

            lefthalfEsp = listaEsp[:mid]
            righthalfEsp = listaEsp[mid:]

            lefthalfNum = listNum[:mid]
            righthalfNum = listNum[mid:]

            self.ordena(lefthalf, lefthalfEsp, lefthalfNum)
            self.ordena(righthalf, righthalfEsp, righthalfNum)

            i = 0
            j = 0
            k = 0
            while i < len(lefthalf) and j < len(righthalf):
                if lefthalf[i] < righthalf[j]:
                    listaAl[k] = lefthalf[i]
                    listaEsp[k] = lefthalfEsp[i]
                    listNum[k] = lefthalfNum[i]

                    i = i + 1
                else:
                    listaAl[k] = righthalf[j]
                    listaEsp[k] = righthalfEsp[j]
                    listNum[k] = righthalfNum[j]
                    j = j + 1
                k = k + 1

            while i < len(lefthalf):
                listaAl[k] = lefthalf[i]
                listaEsp[k] = lefthalfEsp[i]
                listNum[k] = lefthalfNum[i]
                i = i + 1
                k = k + 1

            while j < len(righthalf):
                listaAl[k] = righthalf[j]
                listaEsp[k] = righthalfEsp[j]
                listNum[k] = righthalfNum[j]
                j = j + 1
                k = k + 1

        return listaAl, listaEsp, listNum


    def rankingBS(self, listaEsp, listAlter):


        for j in range(len(listaEsp)-1,0,-1):
            for s in range(j):

                if listAlter[s] > listAlter[s+1]:
                    temp = listAlter[s]
                    listAlter[s] = listAlter[s + 1]
                    listAlter[s + 1] = temp

                    temp2 = listaEsp[s]
                    listaEsp[s] = listaEsp[s + 1]
                    listaEsp[s + 1] = temp2


        return listaEsp, listAlter


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

                #listaLocaisConsFim.append(listaLocCons.copy())

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


"""
root = trieNo('+')
root.nomeDoArquivo = "sbtest.fas"

especiesList = root.leEspecies()
seqAminos = root.leSequencias()

qtdeEspecies = root.contaEspecies()
porceCon = root.porcConservacao(qtdeEspecies)
locaisConservados = root.locaisConserv(porceCon)
indices = root.agrupaSitesCons(locaisConservados)

conservados = root.geraConserv(indices, seqAminos)

root.motifBySize(conservados, indices, 2, 50, 50)

"""
