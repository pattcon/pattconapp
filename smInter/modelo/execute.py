from Bio import SeqIO
import re


from model import trieNo

"""
root = trieNo('+')
root.fileName = "teste003.fas"

especiesList = root.readSpecies()
seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 20
root.conservPerc = 60

locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

proc = root.genMotifs(conservados, indices)

finalList = root.processSequences(proc[0], proc[2], proc[1])



print(finalList[0])
print("Alterações :::: ", finalList[3])

"""


root = trieNo('+')
root.fileName = "01Choanephora_curcubitarum.fas"

especiesList = root.readSpecies()
seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 20
root.conservPerc = 30


locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

proc = root.genMotifs(conservados, indices)

finalList = root.processSequences(proc[0], proc[2], proc[1], proc[4])

print(finalList[2])
print(finalList[1])