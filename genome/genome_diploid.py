import random


def insertion():
     return(''.join([random.choice('AGCT') for i in range(random.randint(1,10))]))


def deletion():
    return(abs(int(random.gauss(4,2))))


def mutation():
    return(random.choice('AGCT'))



def diploid(genome,insertion_rate, deletion_rate, mutation_rate):
    genome_copy = []
    i = 0
    while i < len(genome):
        if random.uniform(0,1) < mutation_rate:
            genome_copy.append(mutation())
            #print 'here'

        elif random.uniform(0,1) < deletion_rate:
            i+= deletion()
            #print 'here2'
        elif random.uniform(0,1) < insertion_rate:
            genome_copy.append(insertion()) 
            #print 'here3'
        else:
            genome_copy.append(genome[i])
        i+=1
    return ''.join([nucl for nucl in genome_copy])
