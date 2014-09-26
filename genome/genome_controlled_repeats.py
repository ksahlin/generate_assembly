import random

def generate_genome(repeat_length, mid_seq_length, nr_repeats):
	"""
	Generates a genome that looks like ARBRC...RN
	where A,B,C,...,N are unique regions with length mid_seq_length and
	R anre a perfect repeat with length repeat_length.
	"""
	genome = ''
	bp = ['A','G','C','T']
	repeat = [random.choice(bp) for nucl in range(repeat_length)]
	for i in range(nr_repeats):
		mid_seq = [random.choice(bp) for nucl in range(mid_seq_length)]
		genome += mid_seq + repeat

	mid_seq = [random.choice(bp) for nucl in range(mid_seq_length)]
	genome += mid_seq

	return genome