import random

def generate_genome(genome_length, nr_repeats, mid_seq_min, mid_seq_max, repeat_density):
    """
    Generates a genome that looks like A_1,...A_n,R,B_1,..,B_n,R,C...R,...
    where A_1,..,A_n,B_1,...,C_1,..., are unique regions with length mid_seq_length and
    R are one out of nr_repeats repeats.
    """
    genome = []
    contigs = {}
    bp = ['A','G','C','T']
    repeats = []

    for i in range(nr_repeats):
    	repeat_length = 100*(i+1)
    	repeat = [random.choice(bp) for nucl in range(repeat_length)]
    	repeats.append(repeat)
        contigs[">{0}_repeat".format(i)] = "".join(r for r in repeat)

    pos = 0
    iterator = 1
    # print(len(repeats))
    while pos < genome_length:
    # for i in range(nr_repeats):
        if random.random() < repeat_density:
            repeat_copy = random.randint(0, nr_repeats-1)
            repeat = repeats[repeat_copy]
            genome.append(repeat)
            pos += len(repeat)
        else:
            ctg_len = int(random.uniform(mid_seq_min, mid_seq_max+1))
            ctg_seq = [random.choice(bp) for nucl in range(ctg_len)]
            genome.append(ctg_seq)
            contigs[">{0}_pos_{1}_unique".format(iterator, pos)] = "".join([nucl for nucl in ctg_seq])
            pos += ctg_len

        iterator += 1

    # mid_seq = [random.choice(bp) for nucl in range(mid_seq_length)]
    # genome += mid_seq
    genome = "".join([nucl for contig in genome for nucl in contig])
    return genome, contigs

