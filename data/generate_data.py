import sys
import random
import numpy as np
import gc # garbage collection; call with gc.collect()

bases = ['G','T','C','A']

def gen_random_seq(seq_len):
	'''
	Generate and return random GTCA sequence of length seq_len.
	'''
	return ''.join(np.random.choice(bases,seq_len))

def gen_repeating_seq(pattern_len, pattern_freq, seq_len):
	'''
	Generate a random GTCA sequence containing 
	a repeating pattern of length pattern_len.
	Pattern occurs pattern_freq amount of times, 
	sequence is seq_len long.
	'''
	pattern = gen_random_seq(pattern_len)

	# Get indices to insert repeating pattern
	remaining_bases = seq_len - (pattern_len * pattern_freq)
	indices = sorted(np.random.randint(0,remaining_bases,pattern_freq))
	padded_bases = \
		[indices[0]] + \
		[t - s for s, t in zip(indices, indices[1:])] + \
		[remaining_bases - indices[-1]]
	
	# Insert random number of bases between 
	seq = gen_random_seq(padded_bases[0])
	for num_bases in padded_bases[1:]:
		seq += pattern + gen_random_seq(num_bases)

	return seq


def perturb_seq(seq, num_deletions, num_insertions, num_swaps):
	pass


if __name__ == '__main__':
	seq_len = int(sys.argv[1])
	# seq = gen_random_seq(seq_len)
	seq = gen_repeating_seq(6,6,100)
	sys.stdout.write(seq)
