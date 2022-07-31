from django.http import FileResponse
from django.shortcuts import render
from .modelo.analise import conserv_alter
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from django.core.files.storage import FileSystemStorage

import os, urllib, io, base64


def index(request):
    return render(request, 'index.html')


def analise(request):

    searchType = request.POST.get("typesearch")
    tsearch = ""
    context = {}

    tsearch = "Minimun Size"
    fastaloc = ''
    percCons = 0
    percCont = 0
    minSize = 0

    if request.method == 'POST' and request.FILES['fileseq']:
        myfile = request.FILES['fileseq']
        fs = FileSystemStorage()

        percCont = int(request.POST.get("contpercent"))
        percCons = int(request.POST.get("conspercent"))
        minSize = int(request.POST.get("tammin"))
        filename = fs.save(myfile.name, myfile)
        fastaloc = str("smamF/fasta/input/" + myfile.name)

    new = conserv_alter()
    new.nomeArq = fastaloc
    new.percCons = percCons
    new.percCont = percCont

    lstMotifs = new.executeBySize(minSize, percCons, percCont)[0]
    lstLocals = new.executeBySize(minSize, percCons, percCont)[1]
    lstSpec = new.executeBySize(minSize, percCons, percCont)[2]
    lstSupprt = new.executeBySize(minSize, percCons, percCont)[3]
    lstSpecN = new.executeBySize(minSize, percCons, percCont)[4]
    lstSizes = new.executeBySize(minSize, percCons, percCont)[5]
    lstEndLocal = new.executeBySize(minSize, percCons, percCont)[6]

    lstSpAlter = new.executeBySize(minSize, percCons, percCont)[7]
    lstSpNoAlter = new.executeBySize(minSize, percCons, percCont)[8]

    lstStrAlter = new.executeBySize(minSize, percCons, percCont)[9]
    lstStrNoAlter = new.executeBySize(minSize, percCons, percCont)[10]

    lstSpList = new.executeBySize(minSize, percCons, percCont)[11]
    lstStrNorm = new.executeBySize(minSize, percCons, percCont)[12]
    lstNumAlt = new.executeBySize(minSize, percCons, percCont)[13]

    request.session['lstLocals'] = lstLocals
    request.session['lstEndLocal'] = lstEndLocal

    request.session['lstSpAlt'] = lstSpAlter
    request.session['lstSpNoAlt'] = lstSpNoAlter

    request.session['lstStrAlt'] = lstStrAlter
    request.session['lstStrNoAlt'] = lstStrNoAlter

    request.session['lstMotifs'] = lstMotifs
    request.session['lstSpec'] = lstSpList
    request.session['lstStrNorm'] = lstStrNorm
    request.session['lstNumAlt'] = lstNumAlt
    request.session['lstSizes'] = lstSizes

    lstMotifNumber, lsNumOrig = [],[]
    for m in range(len(lstMotifs)):
        lstMotifNumber.append(m + 1)
        lsNumOrig.append(m)



        lstQty = []
        for t in range(len(lstSpec)):
            lstQty.append(10)

    request.session['listNumMot'] = lsNumOrig

    lstSpecNames = []
    for n in range(len(lstSpecN)):
        vl = lstSpecN[n]
        lstSpecNames.append(vl)

    lstSpecNumber = []
    for n in range(len(lstSpecNames)):
        lstSpecNumber.append(str(n + 1))

    lstZip = zip(lstMotifs, lstLocals, lstSpec, lstSupprt, lstMotifNumber, lstSpecNames, lstSizes, lstEndLocal)

    context = {"tsearch": tsearch, "minsize": minSize, "motifs": lstMotifs,
               "locals": lstLocals, "spec": lstSpec, "supports": lstSupprt, "ziplist": lstZip
               }

    return render(request, "resultbysize.html", context=context)


def motbyspecies(request):
    lstSpAlter = list(request.session.get('lstSpAlt'))
    lstSpNoAlter = list(request.session.get('lstSpNoAlt'))

    lstStrAlter = list(request.session.get('lstStrAlt'))
    lstStrNoAlter = list(request.session.get('lstStrNoAlt'))

    lstMotifs = list(request.session.get('lstMotifs'))
    lstSpec = list(request.session.get('lstSpec'))





    lstSpAlterFinal, lstSpNoAlterFinal, lstMotifsFinal = [],[], []
    lstStrAlterFinal, lstStrNoAlterFinal = [], []



    lsTemp, lsTemp2, lstNumMotifs, lstIdxMot = [],[], [], []

    for t in range(len(lstMotifs)):
        lstMotifsFinal.append(lstMotifs[t])
        lstNumMotifs.append(t+1)
        lstIdxMot.append(t)





    for j in range(len(lstSpNoAlter)):
        for k in range(len(lstStrNoAlter[j])):
            lsTemp.append(lstSpNoAlter[j][k])
        lstSpNoAlterFinal.append(lsTemp.copy())
        lsTemp.clear()

    for l in range(len(lstStrNoAlter)):
        for m in range(len(lstStrNoAlter[l])):
            lsTemp.append(lstStrNoAlter[l][m])
        lstStrNoAlterFinal.append(lsTemp.copy())
        lsTemp.clear()



    #species with alterations

    for j in range(len(lstSpAlter)):
        for k in range(len(lstSpAlter[j])):
            lsTemp.append(lstSpAlter[j][k])
        lstSpAlterFinal.append(lsTemp.copy())
        lsTemp.clear()

    for p in range(len(lstStrAlter)):
        for q in range(len(lstStrAlter[p])):
            lsTemp.append(lstStrAlter[p][q])
        lstStrAlterFinal.append(lsTemp.copy())
        lsTemp.clear()




    listJoinNoAlt,listJoinAlt, tmp, tmp2=[],[],[],[]
    for d in range(len(lstSpNoAlterFinal)):
        for e in range(len(lstSpNoAlterFinal[d])):
                tmp.append(str(lstSpNoAlterFinal[d][e])+" - "+str(lstStrNoAlterFinal[d][e]))
        listJoinNoAlt.append(tmp.copy())
        tmp.clear()

    for m in range(len(lstSpAlterFinal)):
        for n in range(len(lstSpAlterFinal[m])):
                tmp2.append(str(lstSpAlterFinal[m][n])+" - "+str(lstStrAlterFinal[m][n]))
        listJoinAlt.append(tmp2.copy())
        tmp2.clear()



    listStr = zip(lstIdxMot, lstNumMotifs, lstMotifsFinal, listJoinNoAlt, listJoinAlt)
    listMotsNumb = zip(lstNumMotifs, lstMotifsFinal)
    listMM = zip(lstIdxMot, lstNumMotifs, lstMotifsFinal, listJoinNoAlt, listJoinAlt)






    context = {"lstSpNoAlterFinal": lstSpNoAlterFinal, "lstStrNoAlterFinal": lstStrNoAlterFinal, "lstMotifs": lstMotifsFinal, "listStr": listStr,
               "listSpec": lstSpec, "listMotsNumb": listMotsNumb, "lstMotifsFinal": lstMotifsFinal, "listMM":listMM
               }

    return render(request, "motbysp.html", context=context)


def motiflogo(request):
    nmber = 0

    """
    dir = "smInter/static/logoseq/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))
    """

    if request.POST.get("motlogonum"):
        nmber = request.POST.get("motlogonum")
        idx = int(nmber)

    strMot = list(request.session.get('lstStrNorm'))
    motifs = list(request.session.get('lstMotifs'))

    motifSelect = motifs[idx]


    lstMotStr = strMot[idx]

    fastalogo = "smInter/static/logoseq/fastalogo.fas"

    file = open(fastalogo, "w")
    lisNumbers = []
    for p in range(len(lstMotStr)):
        nums = ">Motif Number " + str(p+1)
        lisNumbers.append(p)
        file.write(nums)
        file.write("\n")
        mot = lstMotStr[p]
        file.write(str(mot))
        file.write("\n")
    file.close()



    crp_sites_df = pd.read_csv("smInter/static/logoseq/fastalogo.fas", comment='>', names=['site'])
    crp_sites_list = crp_sites_df['site'].values

    crp_counts_df = logomaker.alignment_to_matrix(sequences=crp_sites_list, to_type='counts')
    crp_weight_df = logomaker.alignment_to_matrix(sequences=crp_sites_list, to_type='weight')
    crp_prob_df = logomaker.alignment_to_matrix(sequences=crp_sites_list, to_type='probability')
    crp_inf_df = logomaker.alignment_to_matrix(sequences=crp_sites_list, to_type='information')

    dir = "smInter/static/logoseq/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))

    logo = logomaker.Logo(crp_counts_df, color_scheme="skylign_protein")
    plt.savefig("smInter/static/logoseq/lgcount.png")

    logo = logomaker.Logo(crp_weight_df, color_scheme="skylign_protein")
    plt.savefig("smInter/static/logoseq/lgweight.png")

    logo = logomaker.Logo(crp_prob_df, color_scheme="skylign_protein")
    plt.savefig("smInter/static/logoseq/lgprob.png")

    logo = logomaker.Logo(crp_inf_df, color_scheme="skylign_protein")
    plt.savefig("smInter/static/logoseq/lginf.png")

    context = {"strMot": strMot, "number": nmber, "motifSelect": motifSelect, "lstMotStr": lstMotStr}

    return render(request, 'motlogo.html', context=context)


def fastaspecies(request):
    specie = ''

    dir = "smamF/fasta/output/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))

    listMotifs = request.session.get('lstMotifs')
    lstSpec = request.session.get('lstSpec')
    lstLocals = list(request.session.get('lstLocals'))
    lstEndLocal = list(request.session.get('lstEndLocal'))

    lBegin, lEnds = [],[]

    for p in range(len(lstLocals)):
        lBegin.append(lstLocals[p])
        lEnds.append(lstEndLocal[p])


    lstStrNorm = request.session.get('lstStrNorm')
    lstNumAlt = request.session.get('lstNumAlt')



    lSpec = []
    for q in range(len(lstSpec)):
        lSpec.append(lstSpec[q])

    #for q in range(len(lstStrNorm)):
        #print(lstStrNorm)

    lstSpecF, lstSeqF, lstBegF, lstEndF = [],[],[],[]
    cont = 0
    for v in range(len(listMotifs)):
        for s in range(len(lstStrNorm)):
            for t in range(len(lstStrNorm[s])):

                if listMotifs[v] == lstStrNorm[s][t]:
                    lstSpecF.append(lstSpec[t])
                    lstSeqF.append(lstStrNorm[s][t])
                    lstBegF.append(lstLocals[s])
                    lstEndF.append(lstEndLocal[s])



    lstSpecFinal, lstSeqFinal, lstBegFinal, lstEndFinal = [],[],[],[]
    lstSpecTemp, lstSeqTemp, lstBegTemp, lstEndTemp= [],[],[],[]


    for q in range(len(lstSpec)):
        cont = 0
        for r in range(len(lstSpecF)):
            if lstSpec[q] == lstSpecF[r]:
                lstSeqTemp.append(lstSeqF[r])
                lstSpecTemp.append(lstSpecF[r])
                lstBegTemp.append(lstBegF[r])
                lstEndTemp.append(lstEndF[r])
        lstSpecFinal.append(lstSpecTemp.copy())
        lstSeqFinal.append(lstSeqTemp.copy())
        lstBegFinal.append(lstBegTemp.copy())
        lstEndFinal.append(lstEndTemp.copy())

        lstSeqTemp.clear()
        lstSpecTemp.clear()
        lstBegTemp.clear()
        lstEndTemp.clear()


    lNumbers = []

    for j in range(len(listMotifs)):
        lNumbers.append(j+1)


    # Catch the specie name from HTML form
    if request.POST.get("sppsel"):
        specie = request.POST.get("sppsel")

    idx = 0
    for j in range(len(lstSpec)):
        if specie == lstSpec[j]:
            idx = j



    listres = []

    # Create new FASTA file
    fastaf = "smamF/fasta/output/fastasp.fas"


    file = open(fastaf, "w")

    # Populate the FASTA file and the list
    if specie != "":
        listres.append("Specie")
        listres.append(specie)
        listres.append("\n")
        for i in range(len(lstSeqFinal[idx])):

            mot = str(lstSeqFinal[idx][i])
            motins = ">" + specie + " Motif "+str(i) + "\n"
            file.write(motins)
            motins = mot
            file.write(motins)
            file.write("\n")



            listres.append("Motif: "+ mot)

            listres.append("Local: "+str(lstBegFinal[idx][i])+" to " + str(lstEndFinal[idx][i]))
            listres.append("")

        file.close()

    context = {'listsp': listres}
    return render(request, "fastaspec.html", context=context)


def fastageneral(request):
    # Delete other FASTA files
    dir = "smamF/fasta/outputgen/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))

    # Create new FASTA file
    fastafile = "smamF/fasta/outputgen/motifgen.fas"
    txtfile = "smamF/fasta/outputgen/motifgen.txt"

    filegen = open(fastafile, "w")
    filetxt = open(txtfile, "w")

    motifs = request.session.get('lstMotifs')

    # Create a list of motifs, without nulls
    listMotifs = []
    for i in range(len(motifs)):
        cont = 0
        for j in range(len(motifs[i])):
            if motifs[i][j] == "-":
                cont += 1
        if cont == 0:
            listMotifs.append(motifs[i])

    listres = []
    # Populate the FASTA file and the list

    for i in range(1,len(listMotifs)):
        mot = str(listMotifs[i])
        motins = ">Motif Number " + str(i)
        motins += "\n"
        motins += str(mot)
        motins += "\n"

        listres.append(">Motif Number " + str(i))
        listres.append(mot)
        listres.append("")
        listres.append("")
        filegen.write(motins)
        filetxt.write(motins)

    filegen.close()
    filetxt.close()


    del request.session['lstMotifs']
    context = {"fastagen": filegen, 'listgen': listres, "filetxt":filetxt}
    return render(request, "fastagen.html", context=context)


def downfastsp(request):

    return FileResponse(open("smamF/fasta/output/fastasp.fas", "rb"), as_attachment=True, content_type="text/plain")


def downfastgen(request):
    return FileResponse(open("smamF/fasta/outputgen/motifgen.fas", "rb"), as_attachment=True, content_type="text/plain")




def downfastgenTxt(request):
    return FileResponse(open("smamF/fasta/outputgen/motifgen.txt", "rb"), as_attachment=True, content_type="text/plain")


def novapag(request):
    return render(request, "resultbysize.html")


def report(request):


    listMotifs = list(request.session.get('lstMotifs'))
    listMotifBegin= list(request.session.get('lstLocals'))
    listMotifEnd = list(request.session.get('lstEndLocal'))

    listSpecies = list(request.session.get('lstSpec'))
    listStr = list(request.session.get('lstStrAlt'))
    listContAlter = list(request.session.get('lstNumAlt'))
    listSizes = list(request.session.get('lstSizes'))






    fileloc = "smamF/reports/report.txt"

    file = open(fileloc, "w")


    lsTmpStr, lsTmpAlt, lsTmpSpecie, lsTmpSpecieNoAlt , lsTmpNoAlt = [],[],[],[],[]
    lstStr, lstAlt, lstSpecie, lstSpecieNoAlt, lstNoAlt = [], [], [], [],[]
    for i in range(len(listMotifs)):

        for j in range(len(listStr[i])):

            lsTmpStr.append(listStr[i][j])

        lstStr.append(lsTmpStr.copy())

        lsTmpStr.clear()

        for k in range(len(listContAlter[i])):
            if listContAlter[i][k] > 0:
                lsTmpAlt.append(listContAlter[i][k])
                lsTmpSpecie.append(listSpecies[k])
            else:
                lsTmpNoAlt.append(listContAlter[i][k])
                lsTmpSpecieNoAlt.append(listSpecies[k])


        lstAlt.append(lsTmpAlt.copy())
        lsTmpAlt.clear()
        lstSpecie.append(lsTmpSpecie.copy())
        lsTmpSpecie.clear()

        lstNoAlt.append(lsTmpNoAlt.copy())
        lsTmpNoAlt.clear()
        lstSpecieNoAlt.append(lsTmpSpecieNoAlt.copy())
        lsTmpSpecieNoAlt.clear()


    for m in range(len(listMotifs)):
        file.write("::Motif " + str(m) + " - " + str(listMotifs[m]) + '\n')
        file.write("::Local : " + str(listMotifBegin[m]) + ' to ' + str(listMotifEnd[m]))
        file.write(" - Size: " + str(listSizes[m]) + '\n')
        file.write("::Species With Motifs : " + '\n')
        for s in range(len(lstSpecieNoAlt[m])):
            file.write(" - Specie: " + str(lstSpecieNoAlt[m][s]) + "\n")




        file.write("::Species With Alterations : " + '\n')
        for l in range(len(lstStr[m])):

                file.write(" - Specie: " + str(lstSpecie[m][l]) + "\n")


                file.write(" - Sequence: " + str(lstStr[m][l]) + "\n")
                file.write(" - Alterations: "+str(lstAlt[m][l])+"\n")
                file.write("\n")





        file.write("----------------------------------------------------------------------------")
        file.write("\n")

    file.close()

    return FileResponse(open("smamF/reports/report.txt", "rb"), as_attachment=True, content_type="text/plain")


def tutorial(request):
    return render(request, "tutorial.html")


def downsequencetest(request):
    return FileResponse(open("smamF/fasta/seqtest/testsequence.fas", "rb"), as_attachment=True,
                        content_type="text/plain")


def contact(request):
    return render(request, "contact.html")
