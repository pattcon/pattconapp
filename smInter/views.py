
from django.http import FileResponse
from django.shortcuts import render
from .modelo.analise import conserv_alter
import logomaker
import pandas as pd
import matplotlib.pyplot as plt


from django.core.files.storage import FileSystemStorage

import os, urllib, io, base64


def index(request):


    return render(request, 'index.html')


def analise(request):
    # Delete all files from input directory
    searchType = request.POST.get("typesearch")
    tsearch = ""
    context = {}

    if searchType == "largestmot":
        tsearch = "Largests Motifs"

        dir = "smamF/fasta/input/"
        for f in os.listdir(dir):
            os.remove(os.path.join(dir, f))



        # Get input data from template form
        fastaloc = ""
        porcentg = 0
        if request.method == 'POST' and request.FILES['fileseq']:
            myfile = request.FILES['fileseq']
            fs = FileSystemStorage()

            porcentg = int(request.POST.get("porcent"))
            filename = fs.save(myfile.name, myfile)
            fastaloc = str("smamF/fasta/input/" + myfile.name)



        novo = conserv_alter()
        novo.nomeArq = fastaloc
        novo.por = porcentg

        res = novo.executar()



        request.session["finalList"] = res


        rang = len(res[0])
        tam = []

        for i in range(rang):
            inicio = res[3][i]
            fim = res[4][i]
            tam.append(fim - inicio + 1)

        inicio = res[3]
        fim = res[4]
        numero = res[0]
        local = res[1]
        motivo = res[2]
        species = res[9]



        rang = len(res[0])

        motifsizes = []
        for i in range(rang):
            bg = res[3][i]
            ed = res[4][i]
            motifsizes.append(bg - ed + 1)



        porcOcorr = res[16]

        listazip = zip(numero, local, motivo,inicio,fim, tam)

        context = {

            "vl": "fim", "resultado": res, "numero": res[0], "locais": local, "motivo": motivo, "inicio": inicio,
            "fim": fim, "tam": tam, "inFimJuntos": res[5], "listazip": listazip, "spec": species,
            "porcOcorr": porcOcorr,
            "tsearch": tsearch
        }



    else:
        if searchType == "minsizemotif":
            tsearch = "Minimun Size"


            if request.method == 'POST' and request.FILES['fileseq']:
                myfile = request.FILES['fileseq']
                fs = FileSystemStorage()

                porcentg = int(request.POST.get("porcent"))
                txcont = int(request.POST.get("conserv"))
                minSize = int(request.POST.get("tammin"))
                filename = fs.save(myfile.name, myfile)
                fastaloc = str("smamF/fasta/input/" + myfile.name)

            new = conserv_alter()
            new.nomeArq = fastaloc
            new.por = porcentg

            occurr = new.executeBySize(minSize, porcentg, txcont)[0]
            motifs = new.executeBySize(minSize, porcentg, txcont)[1]
            locals = new.executeBySize(minSize, porcentg, txcont)[2]


            #request.session["finalListBySize"] = motifList

            listazip = zip(motifs, occurr, locals)

            context = {"tsearch":tsearch,"motifList":motifs, "minsize":minSize, "listazip":listazip}

            """
                        ","minsize":minSize, localsList":localsList,#"suportsList":suportsList"""



        return render(request, "resultbysize.html", context=context)




    return render(request, "result.html", context=context)








def motivores(request):

    if request.POST:
        numero = request.POST.get("valor","")
        numero = int(numero)
        indice = numero-1

    res = request.session.get('finalList')
    listOrigi = res[15]
    listaStr = res[6]

    listStr = listaStr[indice]
    listaAltStr = res[10]

    motivo = res[2][indice]

    listaContAlt = res[7]
    listaTipo = res[8]
    listaEsp = res[9]
    listaModif = listaAltStr
    listaComp = []

    for j in range(len(listStr)):
        listaComp.append("Espécie:")
        listaComp.append(listaEsp[j])
        listaComp.append("Sequência:")
        listaComp.append(listStr[j])
        listaComp.append("\n")

    fastalogo = "smamF/fasta/output/fastalogo.fas"

    request.session["fastalogo"] = fastalogo

    file = open(fastalogo, "w")
    listaTeste = []
    for p in range(len(listOrigi[indice])):
        nums = ">Motivo Numero "+str(p)
        file.write(nums)
        file.write("\n")
        mot = str(listOrigi[indice][p])
        file.write(mot)
        file.write("\n")
    file.close()



    crp_sites_df = pd.read_csv("smamF/fasta/output/fastalogo.fas", comment='>', names=['site'])
    crp_sites_list = crp_sites_df['site'].values

    crp_counts_df = logomaker.alignment_to_matrix(sequences=crp_sites_list, to_type='counts')
    logo = logomaker.Logo(crp_counts_df, color_scheme="skylign_protein")

    plt.savefig("smInter/static/logoseq/lgseq.png")

    img = "smInter/static/logoseq/lgseq.png"




    contexto = {"numero": numero, "listaStr":listStr, "motivo":motivo, "listaEsp":listaEsp, "listaTipo":listaTipo,
                "listaContaAlt":listaContAlt, "listaMod":listaModif, "listaComp":listaComp, "logo":img,
                }

    return render(request, 'motresumo.html', context=contexto )




def fastaspecies(request):
    specie = ''

    dir = "smamF/fasta/output/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))





    res = request.session.get('finalList')
    motivo = res[2]
    numeros = res[0]



    #Create a list of motifs, without nulls
    listMotifs = []
    listNumbers =[]
    for i in range(len(motivo)):
        cont = 0
        for j in range(len(motivo[i])):
            if motivo[i][j] == "-":
                cont += 1
        if cont == 0:
            listMotifs.append(motivo[i])
            listNumbers.append(i+1)

    #Catch the specie name from HTML form
    if request.POST.get("spp"):
        specie = request.POST.get("spp")

    listres = []

    # Create new FASTA file
    fastaf = "smamF/fasta/output/"+specie+".fas"

    request.session["fastasp"] = fastaf


    file = open(fastaf, "w")

    #Populate the FASTA file and the list
    if specie != "":
        for i in range(len(listMotifs)):
            mot = str(listMotifs[i])
            nums = str(listNumbers[i])
            motins = ">"+specie+" Motivo Numero: "+nums +"\n"
            file.write(motins)
            motins = mot
            file.write(motins)
            file.write("\n")

            listres.append(">"+str(specie))
            listres.append(mot)
            listres.append("")

        file.close()




    context = {"fastasp":file, 'listsp':listres}
    return render(request, "fastaspec.html", context=context)





def fastageneral(request):

    # Delete other FASTA files
    dir = "smamF/fasta/outputgen/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))

    # Create new FASTA file
    fastafile = "smamF/fasta/outputgen/motifgen.fas"



    filegen = open(fastafile, "w")

    res = request.session.get('finalList')
    motivo = res[2]

    #Create a list of motifs, without nulls
    listMotifs = []
    for i in range(len(motivo)):
        cont = 0
        for j in range(len(motivo[i])):
            if motivo[i][j] == "-":
                cont += 1
        if cont == 0:
            listMotifs.append(motivo[i])

    

    listres = []
    #Populate the FASTA file and the list
    
    for i in range(len(listMotifs)):
        mot = str(listMotifs[i])
        motins = ">Motif Number "+ str(i+1)
        motins += "\n"
        motins += str(mot)
        motins += "\n"

        listres.append(">Motif Number "+str(i+1))
        listres.append(mot)
        listres.append("")
        listres.append("")
        filegen.write(motins)
    filegen.close()





    context = {"fastagen":filegen, 'listgen':listres}
    return render(request, "fastagen.html", context=context)



def downfastsp(request):

    fastaspe = request.session.get('fastasp')

    return FileResponse(open(fastaspe, "rb"), as_attachment=True, content_type="text/plain")





def downfastgen(request):


    return FileResponse(open("smamF/fasta/outputgen/motifgen.fas", "rb"), as_attachment=True, content_type="text/plain")


def novapag(request):

    return render(request, "resultbysize.html")


def report(request):
    res = request.session.get('finalList')
    listContAlt = res[15]
    listStr = res[6]
    listSpecies = res[9]
    listAlter = res[10]
    listContAlter = res[7]
    listMotifs = res[2]
    listMotifsNumbers = res[0]
    listMotifBegin = res[3]
    listMotifEnd  = res[4]

    fileloc = "smamF/reports/report.txt"

    file = open(fileloc, "w")



    for i in range(len(listStr)):

        file.write("::Motivo " + str(i) + " - " + str(listMotifs[i]) + '\n')
        file.write("::Locals : " + str(listMotifBegin[i]) + ' a ' + str(listMotifEnd[i]) + '\n')
        file.write("::Species : " + '\n')

        for j in range(len(listStr[i])):
            file.write("Specie: "+ str(listSpecies[j]) + " - Sequência "+ str(listStr[i][j]) + " Alterações: "+ str(listContAlter[i][j]) )
            file.write("\n")

        file.write("----------------------------------------------------------------------------")
        file.write("\n")


    file.close()

    return FileResponse(open("smamF/reports/report.txt", "rb"), as_attachment=True, content_type="text/plain")
