import sys
import random
import os
import argparse
import gc # garbage collection; call with gc.collect()

bases = list('GTCA')

def gen_random_seq(filename, seq_len):
	'''
	Generate and return random GTCA sequence of length seq_len.
	'''
	open(filename,'w').write(''.join(np.random.choice(bases,seq_len))) 

def sample_indices(size, num_indices):
	return sorted(random.sample(range(size),num_indices))

def gen_repeating_seq(filename, pattern_len, pattern_freq, seq_len):
	'''
	Generate a random GTCA sequence containing 
	a repeating pattern of length pattern_len.
	Pattern occurs pattern_freq amount of times, 
	sequence is seq_len long.
	'''
	pattern = gen_random_seq(pattern_len)

	# Get indices to insert repeating pattern
	remaining_bases = seq_len - (pattern_len * pattern_freq)
	indices = sample_indices(remaining_bases,pattern_freq)
	padded_bases = \
		[indices[0]] + \
		[t - s for s, t in zip(indices, indices[1:])] + \
		[remaining_bases - indices[-1]]
	
	# Insert random number of bases between patterns
	seq = gen_random_seq(padded_bases[0])
	for num_bases in padded_bases[1:]:
		seq += pattern + gen_random_seq(num_bases)
	
	return seq


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
	# Construct edited sequence 
	edit_seq = ''
	prev = 0
	for i, edit in edits:
		if edit == 'INS':
			edit_seq += seq[prev:i+1] + random.choice(bases)
		if edit == 'DEL':
			edit_seq += seq[prev:i]
		if edit == 'SWP':
			curr = seq[i]
			exclusion = [x for x in bases if x!=curr]
			edit_seq += seq[prev:i] + random.choice(exclusion)
		prev = i+1
	edit_seq += seq[prev:seq_len]

	return edit_seq

def get_args():
	parser = argparse.ArgumentParser(description="Generate random data")

	parser.add_argument("--random-type", type=str, help="Completely random or with repeating substring (either \'normal\' or \'repeating\'", default=True)
	parser.add_argument("--seq-len", type=int, help="Length of randomly generated sequence")
	parser.add_argument("--pattern-len", type=int, help="Length of repeating substring if random-type = repeating")
	parser.add_argument("--pattern-freq", type=int, help="Frequency of repeating substring if random-type = repeating")
	parser.add_argument("--filename", type=str, help="Filename for edited sequence")

	parser.add_argument("--edit-mode", type=bool, help="Edit sequence if true, else do nothing")
	# parser.add_argument("--edit-file", type=str, help="Destination filename for edited sequence")
	parser.add_argument("--num_inserts", type=int, help="Number of inserts")
	parser.add_argument("--num_deletes", type=int, help="Number of deletes")
	parser.add_argument("--num_swaps", type=int, help="Number of swap")

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
	if args.random_type.lower() == 'normal':
		gen_random_seq(args.filename, args.seq_len)
	elif args.random_type.lower() == 'repeating':
		gen_repeating_seq(args.filename, args.pattern_len, args.pattern_freq, args.seq_len)

	if args.edit_mode == True:
		perturb_seq(args.filename, args.num_inserts, args.num_deletes, args.num_swaps)
	# seq_len = int(sys.argv[1])
	# seq = gen_random_seq(seq_len)
	# seq = gen_repeating_seq(20,500,100000)
	# write_file = sys.argv[4]
	# if write_file:
	# 	with open(write_file, "w") as f:
	# 		f.write(seq)
	# seq = gen_random_seq(20)
	# edit = perturb_seq(seq, 2,5,3)
	# print(seq)
	# print(edit)
	# sys.stdout.write(seq)
