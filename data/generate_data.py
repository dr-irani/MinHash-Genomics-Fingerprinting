import sys
import random
import os
import argparse
import numpy as np
import gc # garbage collection; call with gc.collect()

bases = list('GTCA')

def random_seq(seq_len):
	return ''.join(np.random.choice(bases,seq_len))
def gen_random_seq(filename, seq_len):
	'''
	Generate and return random GTCA sequence of length seq_len.
	'''
	f = open(filename,'w')
	f.write(random_seq(seq_len)) 
	f.close()

def sample_indices(size, num_indices):
	return sorted(random.sample(range(size),num_indices))

def gen_repeating_seq(filename, pattern_len, pattern_freq, seq_len):
	'''
	Generate a random GTCA sequence containing 
	a repeating pattern of length pattern_len.
	Pattern occurs pattern_freq amount of times, 
	sequence is seq_len long.
	'''
	pattern = random_seq(pattern_len)
	print(pattern)

	# Get indices to insert repeating pattern
	remaining_bases = seq_len - (pattern_len * pattern_freq)
	indices = sample_indices(remaining_bases,pattern_freq)
	padded_bases = \
		[indices[0]] + \
		[t - s for s, t in zip(indices, indices[1:])] + \
		[remaining_bases - indices[-1]]
	
	# Insert random number of bases between patterns
	f = open(filename, 'w')
	f.write(random_seq(padded_bases[0]))
	for num_bases in padded_bases[1:]:
		f.write(pattern + random_seq(num_bases))
	
	# return seq


def perturb_seq(filename, num_ins=0, num_del=0, num_swp=0):
	'''
	Given sequence seq, perturb with deletions/insertions/swaps
	'''
	# Generate list of indices to insert edits (random pernutation)
	# Assign edit to each index, then sort
	f = open(filename, 'r')
	seq = f.read()
	f.close()
	seq_len = len(seq)
	num_edits = num_ins + num_del + num_swp
	edit_indices = random.sample(range(seq_len), num_edits)
	edit_labels = (['INS']*num_ins) + (['DEL']*num_del) + (['SWP']*num_swp)
	edits = sorted(zip(edit_indices, edit_labels))

	dst_file = "{}_edited_{}i_{}d_{}s.txt".format(filename, num_ins, num_del, num_swp)
	f = open(dst_file, 'w')
	# Construct edited sequence and write t0 file
	prev = 0
	for i, edit in edits:
		if edit == 'INS':
			f.write(seq[prev:i+1] + random.choice(bases))
		if edit == 'DEL':
			f.write(seq[prev:i])
		if edit == 'SWP':
			curr = seq[i]
			exclusion = [x for x in bases if x!=curr]
			f.write(seq[prev:i] + random.choice(exclusion))
		prev = i+1
	f.write(seq[prev:seq_len])
	f.close()
	# return edit_seq

def get_args():
	parser = argparse.ArgumentParser(description="Generate random data")

	parser.add_argument("--random-type", type=str, help="Completely random or with repeating substring (either \'normal\' or \'repeating\'",)
	parser.add_argument("--seq-len", type=int, help="Length of randomly generated sequence")
	parser.add_argument("--pattern-len", type=int, help="Length of repeating substring if random-type = repeating")
	parser.add_argument("--pattern-freq", type=int, help="Frequency of repeating substring if random-type = repeating")
	parser.add_argument("--filename", type=str, help="Filename for edited sequence")

	# edit mode
	parser.add_argument("--edit-mode", type=bool, help="Edit sequence if true, else do nothing")

	parser.add_argument("--num-inserts", type=int, help="Number of inserts")
	parser.add_argument("--num-deletes", type=int, help="Number of deletes")
	parser.add_argument("--num-swaps", type=int, help="Number of swap")

	# TODO: Add optional command-line arguments as necessary.
	args = parser.parse_args()
	return args

if __name__ == '__main__':
	'''
	python generate_data.py <length> <ins> <del> <swp> <true_filename> <edited_filename>
							1		 2     3	 4 		5 				6
	'''
	args = get_args()
	seq_len = args.seq_len
	if args.random_type == 'normal':
		gen_random_seq(args.filename, args.seq_len)
	if args.random_type == 'repeating':
		gen_repeating_seq(args.filename, args.pattern_len, args.pattern_freq, args.seq_len)

	if args.edit_mode == True:
		perturb_seq(args.filename, args.num_inserts, args.num_deletes, args.num_swaps)

