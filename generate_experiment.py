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

    
    args = parser.parse_args()


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
    lib = generate_library.LogNSimulator(read_length=100, coverage=30, mean=7, std=0.7, distribution='lognormal')
    lib.simulate(genome_file, reads_prefix)




    subprocess.check_call( [ "align", reads_prefix+"_pe1.fa",  reads_prefix+"_pe2.fa", 
        os.path.join(args.output_folder,'genome.fa'), os.path.join(args.output_folder,'mapped_genome'), 
        '-sort', '--mem' ], 
        stdout = open( os.path.join(args.output_folder, 'bwa_mem.genome.stdout'),'w'), 
        stderr = open( os.path.join(args.output_folder, 'bwa_mem.genome.stderr'),'w') )
    subprocess.check_call( [ "align", reads_prefix+"_pe1.fa",  reads_prefix+"_pe2.fa", 
        os.path.join(args.output_folder,'contigs.fa'), os.path.join(args.output_folder,'mapped_contigs'), 
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






