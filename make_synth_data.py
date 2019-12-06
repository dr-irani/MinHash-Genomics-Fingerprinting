from generate_data import *

if __name__ == '__main__':
	seq_len = 8000
	num_seq = 10
	num_edits = [i*100 for i in range(1,11)] # num_edits = [100-1000]
	# insertions, deletions, subs
	# Split 25-25-50 ratio for equal len
	edits = [[int(x/4), int(x/4), int(x/2)] for x in num_edits]
	# Store in 'data' directory
	filenames = ['data/synth_{}.txt'.format(i) for i in range(1,11)]

	for i,filename in enumerate(filenames):
		gen_random_seq(filename, seq_len)
		num_ins, num_del, num_swps = edits[i]
		perturb_seq(filename, num_ins, num_del, num_swps)