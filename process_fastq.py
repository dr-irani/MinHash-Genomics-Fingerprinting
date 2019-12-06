import sys

def parse_fastq(f, fo):
	''' 
	13124 sequences in ecoli_p6_25x.filtered.fastq. 
	Will only take the first 10000 sequences. 
	'''
	reads = []
	i = 0
	max_l = 0
	min_l = 100000

	for line in iter(f.readline, r''):
		name = line[1:].rstrip()
		seq = f.readline().rstrip()
		f.readline()
		qual = f.readline().rstrip()
		reads.append(seq)
		fo.write("%s\n" % seq)

		max_l = max(max_l, len(seq))
		min_l = min(min_l, len(seq))
		i += 1
		
		if i == 10000:
			break

	print(max_l)
	print(min_l)

	return reads

if __name__ == '__main__':

	input_file = 'data/ecoli_p6_25x.filtered.fastq'
	output_file = 'data/ecoli_realdata.txt'
	fi = open(input_file, 'r')
	fo = open(output_file, 'w')
	reads = parse_fastq(fi, fo)
	print(len(reads))
	fi.close()
	fo.close()