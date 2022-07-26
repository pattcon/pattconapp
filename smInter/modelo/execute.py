from Bio import SeqIO
import re
import numpy
import time

from model import trieNo

start_time = time.time()

root = trieNo('+')
root.fileName = "seqUnaligned.fasta"

especiesList = root.readSpecies()

seqAminos = root.readSequences()

qtdeEspecies = root.countSpecies()
root.contentPerc = 50


root.conservPerc = 70


locaisConservados = root.conservLocals(root.contentPerc)
indices = root.groupConsLocals(locaisConservados)

conservados = root.genConserved(indices, seqAminos)

lista = root.motifBySize(indices, conservados, 4,70)[0]

print(lista)


end_time = time.time()
print(end_time-start_time)






