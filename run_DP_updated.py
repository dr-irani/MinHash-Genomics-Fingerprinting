import random
import gc
from edit_distance_DP import edDistDp
import re
import glob

def calculate_metrics(seq1, seq2):
	'''
	Format output data for synthetic data
	'''
	edit_distance = edDistDp(seq1, seq2)
	output = "{} \t {} \t {} \n".format(len(seq1), len(seq2), edit_distance) # 7 outputs
	return output

def calculate_metrics_ecoli(i, j, seq1, seq2):
	'''
	Format output data for ecoli data
	'''
	edit_distance = edDistDp(seq1, seq2)
	output = "{} \t {} \t {} \t {} \n".format((i, j), len(seq1), len(seq2), edit_distance) # 7 outputs
	return output

def simulated_data():
	'''
	Run Needlemen-Wunsch on synthetic data
	'''
	filename = 'data/synthdata_8001bp.txt'
	reads =[]
	with open(filename) as f:
		seq = f.readline()
		while seq:
			reads.append(seq)
			seq = f.readline()

	fo_name = 'output/synthdata_true_editdistance.txt'
	fo = open(fo_name, 'w')
	
	for i in range(0, len(reads), 2):
		output = calculate_metrics(i, i+1, reads[i], reads[i+1])
		fo.write(output)
		print("%d, %d completed" % (i, i+1))
	fo.close()

def real_ecoli_data():
	'''
	Run Needlemen-Wunsch on Ecoli data
	'''
	filename = 'data/ecoli_realdata.txt' # 10000 lines
	reads =[]
	i = 0
	with open(filename) as f:
		seq = f.readline()
		i += 1
		while seq:
			reads.append(seq)
			seq = f.readline()
			if i == 600:
				break
			i += 1

	fo_name = 'output/ecoli_true_editdistance.txt'
	fo = open(fo_name, 'w')
	fo.write("{} \t {} \t {} \t {} \n".format("indexID", "len(seq1)", "len(seq2)", "edit_distance"))
	for i in range(0, len(reads), 2):
		output = calculate_metrics_ecoli(i, i+1, reads[i], reads[i+1])
		fo.write(output)
		print("%d, %d completed" % (i, i+1))
	fo.close()

# Main.
if __name__ == '__main__':
	
	simulated_data()
	print("Finished running simulated data")
	real_ecoli_data()
	print("Finished running real ecoli data")

	