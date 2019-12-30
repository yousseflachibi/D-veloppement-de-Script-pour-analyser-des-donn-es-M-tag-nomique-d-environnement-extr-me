import debruijn as db

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#-----------------		Metthode 1		 ------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------



#--------------------------------
#--------- etap 1 ---------------
#--------------------------------

name_file = "1128.fna"
print("etap 1")
fh = open(name_file,'r')
        #line = fh.read()
line = fh.readline()
x = [[],[]]
meta = list()
sequence = ""

v = 0
while line:
    line = line.rstrip('\n')
    if '>' in line:
        x[0].append(line)
        v = 1
    else:
    	sequence = sequence + line
    	if(v==1):
        	x[1].append(sequence)
        	sequence = ""
        	v = 0
    	
    line = fh.readline()

#print("len x[1]")
#for i in range(len(x[1])):
#	print(x[1][i])

#--------------------------------
#--------- etap 2 ---------------
#--------------------------------
print("etap 2")
# read all protein Id from Meta => x[0] => list()
def get_Protein_name(meta):
	start = "[protein="
	end="]"
	name = (meta.split(start))[1].split(end)[0]
	return name

names_of_proteins = []
for i in range(len(x[0])):
	names_of_proteins.append(get_Protein_name(x[0][i]))

# distinct names of Proteins from => names_of_proteins => list()
def Remove(ducplicate):
	final_list = []
	for num in ducplicate:
		if num not in final_list:
			final_list.append(num)
	return final_list
names_of_proteins_distinct = Remove(names_of_proteins)

#--------------------------------
#--------- etap 3 ---------------
#--------------------------------
print("etap 3")
# given two sequences and an offset, count the number of matching bases
def score(sequence1,sequence2,offset):
    return sum([1 for position in range(max(0-offset,0), min([len(sequence2)-offset, len(sequence2), len(sequence1)-offset])) 
    	if sequence2[position] == sequence1[position+offset]])

# given two sequences, find the offset which gives the best score
def find_best_offset(sequence1,sequence2):
    return max([(score(sequence1,sequence2,offset),offset,sequence2,sequence1) for offset in range(1-len(sequence2),len(sequence1))])

# given a single sequence and a collection of others, find the other sequence with the best match score
def find_best_match(sequence,others):
    return max([find_best_offset(sequence,sequence2) for sequence2 in others if sequence2 != sequence])

# given two sequences and an offset, calculate the consensus
def consensus(score,offset,sequence1,sequence2):
    return sequence2[0:max(0,offset)] + sequence1 +  sequence2[len(sequence1)+offset:]

# given a sequence and collection of others, return the complete consensus using recursion
def assemble(sequence, others):
    return consensus(*find_best_match(sequence, others)) if len(others) == 1 else assemble(consensus(*find_best_match(sequence, others)), [ y for y in others if y != find_best_match(sequence, others)[2]])

# given a collection of sequences, call assemble() to start the recursion
def assemble_helper(dnas):
    return assemble(dnas[0],dnas[1:])

assembler_result = [[],[]]
print("etap 4")
print(len(names_of_proteins_distinct))
print(len(x[0]))
# get all sequence of one proteins names  =>  => list()

sequences = list()
for i in range(0,len(names_of_proteins_distinct)):
	print(i)
	all_index_Protein = [f for f, ff in enumerate(x[0]) if get_Protein_name(ff) == names_of_proteins_distinct[i]]
	for j in range(0,len(all_index_Protein)):
		sequences.append(x[1][all_index_Protein[j]])
	sequences = [k for k in sequences if k]
	file = open("demofile3.fasta", "w")
	for k in range(0,len(sequences)):
		names_proteins = str(">"+names_of_proteins_distinct[i]+"\n")
		seq_proteins = str(sequences[k]+"\n")
		file.write(names_proteins)
		file.write(seq_proteins)
	file.close()
	#assembler_result[0].append(names_of_proteins_distinct[i])
	#assembler_result[1].append(assemble_helper(sequences))
	# Main script
	fname = 'demofile3.fasta'
	reads = db.read_reads(fname)
	# print reads

	test = ['bcdefg', 'defghi', 'abcd']
	# g = construct_graph(test, 3)
	g = db.construct_graph(reads, 11)
	# print_graph(g)
	# for k in g.keys():
	#   print k, g[k]
	# g = construct_graph(reads)
	contig = db.output_contigs(g)
	#print (contig)
	'''
	print(i)
	for j in range(0,len(x[0])):
		name_P = get_Protein_name(x[0][j])
		if name_P == names_of_proteins_distinct[i]:
			sequences.append(x[1][j])
	sequences = [k for k in sequences if k]
	assembler_result[0].append(names_of_proteins_distinct[i])
	assembler_result[1].append(assemble_helper(sequences))
	'''
	#break
#print(assembler_result[0])
print("etap 5")



file = open("Results/Methode-1/Rtl.fasta", "w")
for i in range(0,len(names_of_proteins_distinct)):
	names_proteins = str(assembler_result[0][i]+"\n")
	seq_proteins = str(assembler_result[1][i]+"\n")
	file.write(names_proteins)
	file.write(seq_proteins)
file.close()



#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#-----------------		Metthode 2		 ------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

