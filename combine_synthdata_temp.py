import os
import sys
import glob
''' Quick code for combining dataset '''

filenames1 = ['data/synth_{}.txt'.format(i) for i in range(1,11)]
filenames2 = ['data/synth_{}_edited.txt'.format(i) for i in range(1, 11)]

filename = "synthdata_8001bp.txt"
final = open(filename, "w")

for i in range(1, 11):
	fn1 = filenames1[i-1]
	fn2 = filenames2[i-1]
	with open(fn1) as f:
		seq1 = f.readline().strip()
	with open(fn2) as f:		
		seq2 = f.readline().strip()
	final.write(seq1 + '\n')
	final.write(seq2 + '\n')

filenames1 = ['synth_data/synth_{}.txt'.format(i) for i in range(1,11)]

for i in range(1, 11):
	fn1 = filenames1[i-1]
	filenames2 = glob.glob('synth_data/synth_1_edited_*.txt')

	with open(fn1) as f1:
		seq1 = f1.readline().strip()
		for j in range(5):
			fn2 = filenames2[j]
			with open(fn2) as f2:
				seq2 = f2.readline().strip()
				final.write(seq1 + '\n')
				final.write(seq2 + '\n')

final.close()
