import random
import gc
from edit_distance_DP import edDistDp

def calculate_metrics(i, j, seq1, seq2):

	edit_distance = edDistDp(seq1, seq2)
	output = "{} \t {} \t {} \t {} \n".format((i, j), len(seq1), len(seq2), edit_distance) # 7 outputs
	return output

def simulated_data():
	
	filename = 'data/random_10000bp_50seq.txt'
	reads = []
	with open(filename) as f:
		seq = f.readline()
		while seq:
			reads.append(seq)
			seq = f.readline()

	fo_name = 'random_simulated_true_editdistance.txt'
	fo = open(fo_name, 'w')
	fo.write("{} \t {} \t {} \t {} \n".format("indexID", "len(seq1)", "len(seq2)", "edit_distance"))
	for i in range(len(reads)):
		for j in range(i+1, len(reads)):
			output = calculate_metrics(i, j, reads[i], reads[j])
			fo.write(output)
			print("%d, %d completed" % (i, j))
	fo.close

def real_ecoli_data():

	filename = 'data/ecoli_realdata.txt' # 10000 lines
	reads =[]
	i = 0
	with open(filename) as f:
		seq = f.readline()
		i += 1
		while seq:
			reads.append(seq)
			seq = f.readline()
			if i == 1000:
				break
			i += 1

	fo_name = 'ecoli_realdata_true_editdistance.txt'
	fo = open(fo_name, 'w')
	fo.write("{} \t {} \t {} \t {} \n".format("indexID", "len(seq1)", "len(seq2)", "edit_distance"))
	for i in range(0, len(reads), 2):
		output = calculate_metrics(i, i+1, reads[i], reads[i+1])
		fo.write(output)
		print("%d, %d completed" % (i, i+1))
	fo.close()

if __name__ == '__main__':
	
	simulated_data()
	print("Finished running simulated data")
	real_ecoli_data()
	print("Finished running real ecoli data")

	