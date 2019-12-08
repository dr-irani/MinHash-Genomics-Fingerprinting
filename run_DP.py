import random
import gc
from edit_distance_DP import edDistDp
import re
import glob

def calculate_metrics(seq1, seq2):

	edit_distance = edDistDp(seq1, seq2)
	output = "{} \t {} \t {} \n".format(len(seq1), len(seq2), edit_distance) # 7 outputs
	return output

def simulated_data():
	
	filenames1 = ['synth_data/synth_{}.txt'.format(i) for i in range(1,11)]
	#filenames2 = ['data/synth_{}_edited.txt'.format(i) for i in range(1, 11)]

	for i in range(1, 11):
		fn1 = filenames1[i-1]
		filenames2 = glob.glob('synth_data/synth_1_edited_*.txt')

		with open(fn1) as f1:
			seq1 = f1.readline()
			for j in range(5):
				fn2 = filenames2[j]
				with open(fn2) as f2:
					seq2 = f2.readline()
				fo_name = 'synth_output/synth_{}_{}_editdistance.txt'.format(i, j)
				fo = open(fo_name, 'w')
				output = calculate_metrics(seq1, seq2)
				fo.write(output)
				print("%d completed" % i)
				fo.close

	# with open(filename) as f:
	# 	seq = f.readline()
	# 	while seq:
	# 		reads.append(seq)
	# 		seq = f.readline()

	# fo_name = 'random_simulated_true_editdistance.txt'
	# fo = open(fo_name, 'w')
	# fo.write("{} \t {} \t {} \t {} \n".format("indexID", "len(seq1)", "len(seq2)", "edit_distance"))
	# for i in range(len(reads)):
	# 	for j in range(i+1, len(reads)):
	# 		output = calculate_metrics(i, j, reads[i], reads[j])
	# 		fo.write(output)
	# 		print("%d, %d completed" % (i, j))
	# fo.close

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

# Main.
if __name__ == '__main__':
	
	simulated_data()
	print("Finished running simulated data")
	# real_ecoli_data()
	# print("Finished running real ecoli data")

	