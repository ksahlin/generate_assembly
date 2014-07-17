import sys
import os
import argparse
import subprocess

from genome import genome_markov
from reads import generate_library
from mapping import align
from assembly import minia

#from simulate import simulate_v2 as sim



class ConfigParams(object):
    """docstring for ConfigParams"""
    def __init__(self):
        super(ConfigParams, self).__init__()
        self.libs = []

    def read_cfg(self,config_file):
        for line in open(config_file,'r'):
            if line[:1] == "1":
                self.output_path = line.split()[1] 
            
            elif line[:1] == "2":
                self.genome_length = int(line.split()[1])
            
            elif line[:1] == "3":
                input_lib = line.split()
                if  input_lib[1] == "pe":
                    read_length, coverage, mean, sd = map(lambda x: int(x), input_lib[2:])
                    lib = generate_library.DNAseq(read_length, coverage, mean, sd)
                elif input_lib[1] == "mp":
                    read_length, coverage, mean, sd, c_ratio, c_mean,c_sd = map(lambda x: int(x), input_lib[2:])
                    c_ratio = float(c_ratio)/100.0
                    lib = generate_library.DNAseq(read_length, coverage, mean, sd, contamination_rate = c_ratio, contamine_mean= c_mean, contamine_std_dev= c_sd)
                self.libs.append(lib)

            elif line[:1] == "4":
                assembler = line.split()[1]   

            elif line[:1] == "5":
                self.aligner, kmer_size, abundance = tuple(line.split()[1:])
                self.kmer_size, self.abundance = kmer_size, abundance
     




def main():
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Generate a scaffold experiment.")
    parser.add_argument('config_file', type=str, help='Path to config file. ')
    # # Genome sequence options
    # parser.add_argument('--diploid', action = "store_true", help='Generate a diploid genome. (default is haploid) ')
    # parser.add_argument('-g', dest='genome',type=float, nargs=5, help='Genome stats as: bp composition, length e.g. 0.25 0.25 0.25 0.25 1000000. ')
  
    # # Genomic libraries
    # parser.add_argument('nr_pe', type=int, help='Number of PE libraries. ')
    # parser.add_argument('nr_mp', type=int, help='Number of MP libraries. ')

    # parser.add_argument('mean_pe', type=int, nargs='+', help='PE mean lib insert size. ')
    # parser.add_argument('std_dev_pe', type=int, nargs='+', help='PE lib Standard deviation of insert size ')

    # parser.add_argument('mean_mp', type=int, nargs='+', help='MP mean lib insert size. ')
    # parser.add_argument('std_dev_mp', type=int, nargs='+', help='MP lib Standard deviation of insert size ')

    # parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    # parser.add_argument('coverage', type=int, nargs='+', help='Coverage for read library. ')

    # parser.add_argument('-c_rate',dest='c_rate', type=float, default=False, nargs ='+',  help='If contamine is specified, \
    #     the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
    #     this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')
    # parser.add_argument('-c_mean',dest='c_mean', type=int, default=False, nargs ='+',  help='If contamine is specified, \
    #     the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
    #     this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')
    # parser.add_argument('-c_stddev',dest='c_stddev', type=int, default=False, nargs ='+',  help='If contamine is specified, \
    #     the library will contian a ration of reads orientated in RC order with a fragment distribution of mu and sigma.\
    #     this argument should be specified as -c 0.1 500 50 if contamine level of 10%  from a N(500,50) distribution. ')

 
    # # contig generation

    # parser.add_argument('-ctgs', dest='contigs', default=False, type=float, nargs=5, help='Generate contigs as min_size, max_size, min_gap, max_gap large_ctg_weight. ')
    # parser.add_argument('--minia', dest='minia', default=False, action="store_true", help='Generate contigs as min_size, max_size, min_gap, max_gap large_ctg_weight. ')


  
    # parser.add_argument('outfolder', type=str, help='Path to output location (folder). ')


 
    args = parser.parse_args()

    cfg = ConfigParams()
    cfg.read_cfg(args.config_file)

    if not os.path.exists(cfg.output_path):
        os.makedirs(cfg.output_path)

    #genome
    genome = genome_markov.Genome([0.25]*4, 200, '>genome_{0}bp'.format(cfg.genome_length))
    genome.sequence = genome_markov.generate_genome(genome.sequence, cfg.genome_length)
    genome_file = open(os.path.join(cfg.output_path,'haplotype1.fa'), 'w')
    genome_file.write(genome.genome_fasta_format())

    kmer_count_file = open(os.path.join(cfg.output_path,'kmer_count.txt'), 'w')
    print >> kmer_count_file, 'k=31' , genome_markov.count_kmers(genome.sequence, 31)
    print >> kmer_count_file, 'k=41', genome_markov.count_kmers(genome.sequence, 41)
    print >> kmer_count_file, 'k=51', genome_markov.count_kmers(genome.sequence, 51)
    print >> kmer_count_file, 'k=61', genome_markov.count_kmers(genome.sequence, 61)

    #reads

    for lib in cfg.libs:
        if lib.contamination_rate:
            lib.simulate_mp_reads(genome)
            read1_file = open(os.path.join(cfg.output_path,'MP.{0}.1.fa'.format(lib.lib_mean)), 'w')
            read2_file = open(os.path.join(cfg.output_path,'MP.{0}.2.fa'.format(lib.lib_mean)), 'w')
            read_iterator = lib.fasta_format()
            for read1 in read_iterator:
                read2 = read_iterator.next()
                print >>read1_file, read1, 
                print >>read2_file, read2, 
        else:
            lib.simulate_pe_reads(genome)
            read1_file = open(os.path.join(cfg.output_path,'PE.{0}.1.fa'.format(lib.lib_mean)), 'w')
            read2_file = open(os.path.join(cfg.output_path,'PE.{0}.2.fa'.format(lib.lib_mean)), 'w')
            read_iterator = lib.fasta_format()
            for read1 in read_iterator:
                read2 = read_iterator.next()
                print >>read1_file, read1, 
                print >>read2_file, read2, 


    # assembly
    subprocess.check_call( [ "minia", read1_file.name, cfg.kmer_size, cfg.abundance, '2000000', 
        os.path.join(cfg.output_path,'minia')], 
        stdout = open( os.path.join(cfg.output_path,'minia.stdout'),'w'), 
        stderr = open( os.path.join(cfg.output_path,'minia.stderr'),'w') )




    # ctg_file = open(os.path.join(args.outfolder,'contigs.fa'), 'w')
    # map_file = os.path.join(args.outfolder,'mapped')


    # DNA_lib = sim.DNAseq(args.mean,args.std_dev,args.coverage)




    # ctg_generator = sim.generate_contigs_two_sizes(genome.sequence,int(args.contigs[0]),int(args.contigs[1]), int(args.contigs[2]), int(args.contigs[3]), args.contigs[4])
    # for ctg in ctg_generator:
    #     print >> ctg_file, ctg,

    # #with open( pe1_output, "w" ) as pe1_file:
    # subprocess.check_call( [ "align", read1_file.name, read2_file.name, ctg_file.name, map_file , '-sort' ], 
    #     stdout = open( os.path.join(args.outfolder,'mapped.1'),'w'), stderr = open( os.path.join(args.outfolder,'mapped.0'),'w') )


if __name__=='__main__':
    
    main()






