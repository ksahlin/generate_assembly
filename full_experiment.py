import sys
import os
import argparse
import subprocess

from simulate import simulate_v2 as sim

def main():
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Generate a folder with, genome (diploid or haploid), reads (with\
        or without contamination), contigs and mapped-sorted bam files.")
    
    # Genome sequence options
    parser.add_argument('--diploid', action = "store_true", help='Generate a diploid genome. (default is haploid) ')
    parser.add_argument('-g', dest='genome',type=float, nargs=5, help='Genome stats as: bp composition, length e.g. 0.25 0.25 0.25 0.25 1000000. ')
  
    # Genomic libraries
    parser.add_argument('nr_pe', type=int, help='Number of PE libraries. ')
    parser.add_argument('nr_mp', type=int, help='Number of MP libraries. ')

    parser.add_argument('mean_pe', type=int, nargs='+', help='PE mean lib insert size. ')
    parser.add_argument('std_dev_pe', type=int, nargs='+', help='PE lib Standard deviation of insert size ')

    parser.add_argument('mean_mp', type=int, nargs='+', help='MP mean lib insert size. ')
    parser.add_argument('std_dev_mp', type=int, nargs='+', help='MP lib Standard deviation of insert size ')

    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    parser.add_argument('coverage', type=int, nargs='+', help='Coverage for read library. ')

    parser.add_argument('-c_rate',dest='c_rate', type=float, default=False, nargs ='+',  help='If contamine is specified, \
        the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
        this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')
    parser.add_argument('-c_mean',dest='c_mean', type=int, default=False, nargs ='+',  help='If contamine is specified, \
        the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
        this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')
    parser.add_argument('-c_stddev',dest='c_stddev', type=int, default=False, nargs ='+',  help='If contamine is specified, \
        the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
        this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')

 
    # contig generation

    parser.add_argument('-ctgs', dest='contigs', default=False, type=float, nargs=5, help='Generate contigs as min_size, max_size, min_gap, max_gap large_ctg_weight. ')
    parser.add_argument('--minia', dest='minia', default=False, action="store_true", help='Generate contigs as min_size, max_size, min_gap, max_gap large_ctg_weight. ')


  
    parser.add_argument('outfolder', type=str, help='Path to output location (folder). ')


 
    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    genome_file = open(os.path.join(args.outfolder,'haplotype1.fa'), 'w')
    read1_file = open(os.path.join(args.outfolder,'PE1.fa'), 'w')
    read2_file = open(os.path.join(args.outfolder,'PE2.fa'), 'w')
    ctg_file = open(os.path.join(args.outfolder,'contigs.fa'), 'w')
    map_file = os.path.join(args.outfolder,'mapped')



    if args.diploid:
        pass
    else:
        genome = sim.Genome(args.genome[:-1],int(args.genome[-1]),'genome_copy_1')
        fasta_string = genome.genome_fasta_format()
        genome_file.write(fasta_string)

    if args.contamine:
        DNA_lib = sim.DNAseq(args.mean,args.std_dev,args.read_length, args.contamine[0],args.contamine[1],args.contamine[2])
        DNA_lib.simulate_pe_reads(genome,args.coverage)
    else:
        DNA_lib = sim.DNAseq(args.mean,args.std_dev,args.coverage)

    read_iterator = DNA_lib.fasta_format()
    for read1 in read_iterator:
        read2 = read_iterator.next()
        print >>read1_file, read1, 
        print >>read2_file, read2, 


    ctg_generator = sim.generate_contigs_two_sizes(genome.sequence,int(args.contigs[0]),int(args.contigs[1]), int(args.contigs[2]), int(args.contigs[3]), args.contigs[4])
    for ctg in ctg_generator:
        print >> ctg_file, ctg,

    #with open( pe1_output, "w" ) as pe1_file:
    subprocess.check_call( [ "align", read1_file.name, read2_file.name, ctg_file.name, map_file , '-sort' ], 
        stdout = open( os.path.join(args.outfolder,'mapped.1'),'w'), stderr = open( os.path.join(args.outfolder,'mapped.0'),'w') )


if __name__=='__main__':
    
    main()






