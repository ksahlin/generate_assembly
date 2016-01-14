import sys
import os
import argparse
import subprocess

from genome import genome_controlled_repeats
from reads import generate_library
from mapping import align
from assembly import minia
from svsim.reads import simulator

def main():
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Generate a scaffold experiment.")
    parser.add_argument('output_folder', type=str, help='Path to output folder. ')
    # parser.add_argument('nr_pe', type=int, help='Number of PE libraries. ')
    # parser.add_argument('nr_mp', type=int, help='Number of MP libraries. ')

    # Genome sequence options
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

    # cfg = ConfigParams()
    # cfg.read_cfg(args.config_file)

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    genome_file = open(os.path.join(args.output_folder,'genome.fa'), "r+")
    contig_file = open(os.path.join(args.output_folder,'contigs.fa'), 'w')
    reads_prefix = os.path.join(args.output_folder,'reads')

    #genome
    genome_length, nr_repeats, mid_seq_min, mid_seq_max, repeat_density = 500000, 7, 500, 6000, 0.1
    genome, contigs = genome_controlled_repeats.generate_genome(genome_length, nr_repeats, mid_seq_min, mid_seq_max, repeat_density)
    genome_file.write(">genome_{0}\n{1}".format(len(genome), genome))

    for acc, seq in contigs.iteritems():
        contig_file.write(">{0}\n{1}\n".format(acc, seq))
    genome_file.seek(0)
    lib = generate_library.LogNSimulator(read_length=100, coverage=30, mean=7, std=1.08, distribution='lognormal')
    lib.simulate(genome_file, reads_prefix)

    # lib.simulate_pe_reads(genome_file)
    # read1_file = open(os.path.join(cfg.output_path,'PE.{0}.1.fa'.format(lib.lib_mean)), 'w')
    # read2_file = open(os.path.join(cfg.output_path,'PE.{0}.2.fa'.format(lib.lib_mean)), 'w')
    # read_lib_paths.append((read1_file,read2_file))
    # read_iterator = lib.fasta_format()
    # for read1 in read_iterator:
    #     read2 = read_iterator.next()
    #     print >>read1_file, read1, 
    #     print >>read2_file, read2, 




    #bwa mem
    # for i,(read1_file,read2_file) in enumerate(read_lib_paths):
        #subprocess.check_call( [ "python","mapping/align.py"])
    subprocess.check_call( [ "align", reads_prefix+"_pe1.fa",  reads_prefix+"_pe2.fa", 
        os.path.join(args.output_folder,'genome.fa'), os.path.join(args.output_folder,'mapped_genome'), 
        '-sort', '--mem' ], 
        stdout = open( os.path.join(args.output_folder, 'bwa_mem.genome.stdout'),'w'), 
        stderr = open( os.path.join(args.output_folder, 'bwa_mem.genome.stderr'),'w') )
    subprocess.check_call( [ "align", reads_prefix+"_pe1.fa",  reads_prefix+"_pe2.fa", 
        os.path.join(args.output_folder,'genome.fa'), os.path.join(args.output_folder,'mapped_contigs'), 
        '-sort', '--mem' ], 
        stdout = open( os.path.join(args.output_folder, 'bwa_mem.contigs.stdout'),'w'), 
        stderr = open( os.path.join(args.output_folder, 'bwa_mem.contigs.stderr'),'w') )






    # ctg_generator = sim.generate_contigs_two_sizes(genome.sequence,int(args.contigs[0]),int(args.contigs[1]), int(args.contigs[2]), int(args.contigs[3]), args.contigs[4])
    # for ctg in ctg_generator:
    #     print >> ctg_file, ctg,

    # #with open( pe1_output, "w" ) as pe1_file:
    # subprocess.check_call( [ "align", read1_file.name, read2_file.name, ctg_file.name, map_file , '-sort' ], 
    #     stdout = open( os.path.join(args.outfolder,'mapped.1'),'w'), stderr = open( os.path.join(args.outfolder,'mapped.0'),'w') )


if __name__=='__main__':
    
    main()






