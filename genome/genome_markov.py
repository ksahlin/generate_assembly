# -*- coding: utf-8 -*-



import random
#import numpy
#import sys
#import check_GC_content
import argparse
import Genome2
from collections import Counter

# Markov Chain: Transition matrix
transition = {
    'Sub' : {'Sub':0.8, 'TDup':0.08, 'IDup':0.02, 'Ins':0.05,'Inv':0.05},
    'TDup' : {'Sub':0.8, 'TDup':0.06, 'IDup':0.03, 'Ins':0.055,'Inv':0.055},
    'IDup' : {'Sub':0.8, 'TDup':0.085, 'IDup':0.01, 'Ins':0.0525,'Inv':0.0525},
    'Ins' : {'Sub':0.8, 'TDup':0.085, 'IDup':0.025, 'Ins':0.03,'Inv':0.06},
    'Inv' : {'Sub':0.8, 'TDup':0.085, 'IDup':0.025, 'Ins':0.06,'Inv':0.03}}

##transition = {
##    'Sub' : {'Sub':0.1, 'TDup':0.4, 'IDup':0.2, 'Ins':0.1,'Inv':0.2},
##    'TDup' : {'Sub':0.1, 'TDup':0.25, 'IDup':0.25, 'Ins':0.15,'Inv':0.25},
##    'IDup' : {'Sub':0.1, 'TDup':0.45, 'IDup':0.1, 'Ins':0.1,'Inv':0.25},
##    'Ins' : {'Sub':0.1, 'TDup':0.425, 'IDup':0.225, 'Ins':0.05,'Inv':0.2},
##    'Inv' : {'Sub':0.1, 'TDup':0.45, 'IDup':0.25, 'Ins':0.1,'Inv':0.1}}





def Sub(genome):

# Substitutes one nt to another
# Arguments: Genome

    #print('sub')
    gen = len(genome)
    number = random.randrange(0,gen)
    substitutes = ['A','C','G','T']
    replace = random.choice(substitutes)
    genome[number] = replace
    return(genome)
    





def TDup(genome):

# Makes tandem duplications; repeats of a sequence next to it
# Arguments: Genome

    #print('T')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position

    stop = position + copy_length
    copy = [ genome[u] for u in range(position,stop) ]
    # for u in range (position,stop):
    #     copy.append(genome[u])
    #     u += 1
    genome[stop:stop] = copy
    return(genome)






def IDup(genome):
    
# Makes an interspread duplication; repeats of a sequence on another place in the genome
# Arguments: Genome

    #print('I')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position
    stop = position + copy_length
    copy = [ genome[u] for u in range(position,stop) ]
    # for u in range (position,stop):
    #     copy.append(genome[u])
    #     u += 1
    inter = random.randrange(0,gen)
    genome[inter:inter] = copy
    return(genome)






def Ins(genome):

# Inserts a sequence
# Arguments: Genome

    #print('Ins')
    gen = len(genome)
    #copy = []
    copy_length = abs(int(random.expovariate(1/10.0)))
    #gen2 = len(copy)
    nucleobases = ['A','C','G','T']
    copy = [ random.choice(nucleobases) for i in range(copy_length) ]
    # while gen2 <= copy_length-1:
    #     copy.append(random.choice(nucleobases))
    #     gen2 = len(copy)
    inter = random.randrange(0,gen)
    genome[inter:inter] = copy
    return(genome)
    





def Inv(genome):

# Inserts a inverted sequence in the genome
# Arguments: Genome

    #print('Inv')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position
    stop = position + copy_length
    # print genome[position:stop]
    # for u in range (position,stop):
    #     copy.insert(0,genome[u])
    #     u += 1
    # print copy
    copy = [genome[nucl] for nucl in range(stop -1, position -1, -1)]

    #print inv
    #assert inv == copy
    inter = random.randrange(0,gen)
    #print genome
    genome[inter:inter] = copy
    return genome




def probability_func(state):

# Markov chain - change current state depending on what the last current state was. Uses the transition matrix
# Argument: dictionary with states and probabilities

    count = 0
    choice = random.random()
    for state_key, prob in state.iteritems():
        count += prob
        if choice <= count:
            current_state = state_key
            return current_state



def new_state(genome, finish_length):

# Sends the genome to the different functions above depending on the current state
# Arguments: The genome
    
    current_state = 'Sub'
    genome_length = len(genome)
    genome = list(genome)
    counter = 0
    while genome_length < finish_length:
        counter += 1
        if counter % 100 == 0:
            print genome_length
        current_state = probability_func(transition[current_state])
        if current_state == 'Sub':
            genome = Sub(genome)
            genome_length = len(genome)
        elif current_state == 'TDup':
            genome = TDup(genome)
            genome_length = len(genome)
        elif current_state == 'IDup':
            genome = IDup(genome)
            genome_length = len(genome)
        elif current_state == 'Ins':
            genome = Ins(genome)
            genome_length = len(genome)
        elif current_state == 'Inv':
            genome = Inv(genome)
            genome_length = len(genome)
        else:
            print('something is wrong')
    #print('hej3')
    return genome


def count_kmers(genome,k):
    kmer_counter= Counter()
    for i in range(len(genome)-k+1):
        kmer_counter[genome[i:i+k]] += 1
        if i % 10000 == 0:
            print i

    print kmer_counter.most_common(10)



def main(args):

# Makes a random genome with 500 nt and then increases the length of the genome using the functions above.
# Arguments: length of the random genome, length of the finished genome)

    nucleobases = ['A','C','G','T']
    i = 0
    genome = ''
    while i < args.start_length:
        rand = random.choice(nucleobases)
        genome = genome + rand
        i += 1
    length_genome = len(genome)
    # print genome
    print ('length of genome:', length_genome)
    genome = new_state(genome, args.finish_length)


    genome_str = "".join(genome)
    #print('hej1')
    outfile = open(args.outfile,'w')
    print >>outfile, genome_str
    count_kmers(genome_str,2000)
    # print genome
    #print ('length of genome', len(genome))
    #GC_content = check_GC_content.check_GC_content(genome)
    #print GC_content

    #genome2_str = Genome2.diploid(genome, 0.0001, 0.0001, 0.0001)    
    #print('hej2')
    #outfile2 = open(args.outfile+'copy2.fa','w')
    #print >>outfile2, genome2_str




if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "Simulate a genome of desired length")
    parser.add_argument('start_length', type=int, help='The length of the first genome that is random. ')
    parser.add_argument('finish_length', type=int, help='The approximate size of the finished genome. ')
    parser.add_argument('outfile', type=str, help='outfile prefix. ')

    args = parser.parse_args()
    main(args)



