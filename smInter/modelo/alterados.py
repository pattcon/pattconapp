
from model import trieNo
import re

class conserv_alter:

    def executar(self, arquivo):
        root = trieNo('+')
        root.nomeDoArquivo = arquivo

        especiesList = root.leEspecies()
        seqAminos = root.leSequencias()

        qtdeEspecies = root.contaEspecies()
        porceCon = root.porcConservacao(qtdeEspecies)
        locaisConservados = root.locaisConserv(porceCon)
        indices = root.agrupaSitesCons(locaisConservados)

        conservados = root.geraConserv(indices, seqAminos)
        predominantes = root.geraFim(conservados, indices)



    # gera arquivo de saida

        arquivo = "generatedfiles/saidaAlt.txt"
        with open(arquivo, "w") as arq:


            text09 = '\nSequências Conservadas e Alteradas\n'
            print('\nSequências Conservadas e Alteradas\n')


            arq.write(text09)

            for i in range(len(indices)) :
                listResult = root.compara(predominantes[i], conservados[i], especiesList, indices)






                strT = str(listResult[5])
                seqT = re.search("_", strT)

                if seqT != None:
                    seqAd = "Não Existe Uma Sequência Consenso"
                else:
                    seqAd = listResult[5]
                print('Sequencia: ', seqAd)
                print('Locais: ', indices[i])
                print('Tamanho da Sequência: ', len(indices[i]))





                espConserv = []
                espAlter = []
                numAlt = []


                for j in range(len(listResult[2])):

                    espConserv.append(listResult[2][j])



                    listaEspAlt = listResult[3]
                    listaNumAlt = listResult[6]



                    #Ordenação da lista
                    for m in range(len(listaNumAlt)):
                        for n in range(len(listaNumAlt)-1):
                            if listaNumAlt[n] > listaNumAlt[n + 1]:
                                temp = listaNumAlt[n]
                                listaNumAlt[n] = listaNumAlt[n+1]
                                listaNumAlt[n + 1] = temp

                                temp = listaEspAlt[0][n]
                                listaEspAlt[0][n] = listaEspAlt[0][n + 1]
                                listaEspAlt[0][n + 1] = temp



                espConserv.clear()

                espAlter.clear()

                numAlt.clear()









execu = conserv_alter()
execu.executar("teste.fas")







