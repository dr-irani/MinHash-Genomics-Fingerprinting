import generate_data as gen
from skbio import alignment

if __name__ == '__main__':
	# Get sequence pair
	orig_seq = gen.gen_random_seq(1000)
	edit_seq = 'GTCA' + orig_seq[4:]

	matrix = [[1]*(i) + [0] + [1]*(3-i) for i in range(4)]
	bases = ['G','T','C','A']
	transitions = {}
	for i,u in enumerate(bases):
		transitions[u] = {}
		for j,v in enumerate(bases):
			transitions[u][v] = matrix[i][j]

	info = alignment.global_pairwise_align(orig_seq, edit_seq, 1,1,transitions)
	print(info)