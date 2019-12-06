from generate_data import *
import random

if __name__ == '__main__':
	sigma = 10
	seq_len = 8000
	num_seq = 10
	num_edits = [i*100 for i in range(1,11)] # num_edits = [100-1000]
	seq_per_edit = 5
	# insertions, deletions, subs
	# Split 25-25-50 ratio for equal len
	edits = [[int(x/4), int(x/4), int(x/2)] for x in num_edits]
	# Store in 'data' directory
	filenames = ['synth_data/synth_{}.txt'.format(i) for i in range(1,11)]

	for i,filename in enumerate(filenames):
		gen_random_seq(filename, seq_len)
		num_ins, num_del, num_swps = edits[i]
		additional_swps = random.sample(range(-1*sigma, sigma),seq_per_edit)
		for additional_swp in additional_swps:
			new_swps = num_swps + additional_swp # add variance to swp
			perturb_seq(filename, num_ins, num_del, new_swps)