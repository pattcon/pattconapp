from Bio import SeqIO
import re
import numpy


from model import trieNo

root = trieNo('+')
root.fileName = "testeBySize.fas"

especiesList = root.readSpecies()

print("List species :", especiesList)
seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 20
root.conservPerc = 30


locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

root.motifBySize(conservados, indices, 3,80,30)









