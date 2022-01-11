
from django.http import FileResponse
from django.shortcuts import render
from .modelo.analise import conserv_alter



from django.core.files.storage import FileSystemStorage

import os


def index(request):


    return render(request, 'index.html')


def analise(request):
    # Delete all files from input directory


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





    listazip = zip(numero, local, motivo,inicio,fim, tam)



    context = {

        "vl": "fim", "resultado": res, "numero": res[0], "locais": local, "motivo": motivo, "inicio": inicio,
        "fim":fim, "tam":tam, "inFimJuntos": res[5],  "listazip":listazip, "spec":species
    }

    return render(request, "result.html", context=context)








def motivores(request):

    if request.POST:
        numero = request.POST.get("valor","")
        numero = int(numero)
        indice = numero-1

    res = request.session.get('finalList')
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




    contexto = {"numero": numero, "listaStr":listStr, "motivo":motivo, "listaEsp":listaEsp, "listaTipo":listaTipo,
                "listaContaAlt":listaContAlt, "listaMod":listaModif, "listaComp":listaComp
                }

    return render(request, 'motresumo.html', context=contexto )




def fastaspecies(request):
    specie = ''

    dir = "smamF/fasta/output/"
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))



    #Create new FASTA file
    fastaf = "smamF/fasta/output/motifspec.fas"

    file = open(fastaf, "w")

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

    #Catch the specie name from HTML form
    if request.POST.get("spp"):
        specie = request.POST.get("spp")

    listres = []


    #Populate the FASTA file and the list
    if specie != "":
        for i in range(len(listMotifs)):
            mot = str(listMotifs[i])
            motins = ">"+specie+"\n"
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


    return FileResponse(open("smamF/fasta/output/motifspec.fas", "rb"), as_attachment=True, content_type="text/plain")





def downfastgen(request):


    return FileResponse(open("smamF/fasta/outputgen/motifgen.fas", "rb"), as_attachment=True, content_type="text/plain")