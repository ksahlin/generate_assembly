#!/usr/bin/env python 

import argparse
import os
import subprocess
import sys
import tempfile
import shutil

import pysam
from time import time

##
# Converts a sam file to a bam file using pysam.
#
# @param sam_path Path of the .sam file.
# @param bam_path Path of the resulting .bam file.
#
def sam_to_bam(sam_path, bam_path):
    sam_file = pysam.Samfile( sam_path, "r" )
    bam_file = pysam.Samfile( bam_path, "wb", template = sam_file )

    for alignment in sam_file:
        bam_file.write( alignment )

##
# Maps the given paired end reads using bwa, and writes a
# sorted .bam file in the given output file.
#
# @param pe1_path Path of the first reads.
# @param pe2_path Path of the second reads.
# @param genome_path Path to the reference genome.
# @param output_path Path of the output file without extension ".bam".
#
def map_paired_reads(pe1_path, pe2_path, genome_path, output_path, args):
    work_dir = tempfile.mkdtemp( )
    genome_db = os.path.join( work_dir, "genome" )
    pe1_output = os.path.join( work_dir, "pe1.sai" )
    pe2_output = os.path.join( work_dir, "pe2.sai" )
    bwa_output = os.path.join( work_dir, "output.sam" )
    
    null = open( "/dev/null" ) #open("/tmp/bwa_out")#
    subprocess.check_call( [ "bwa", "index", "-p", genome_db, genome_path ], stderr = null )
    with open( pe1_output, "w" ) as pe1_file:
        subprocess.check_call( [ "bwa", "aln", genome_db, pe1_path ], stdout = pe1_file, stderr = null )
    
    with open( pe2_output, "w" ) as pe2_file:
        subprocess.check_call( [ "bwa", "aln", genome_db, pe2_path ], stdout = pe2_file, stderr = null )
    
    with open( bwa_output, "w" ) as bwa_file:
        subprocess.check_call( [ "bwa", "sampe",
                                "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                                genome_db,
                                pe1_output, pe2_output,
                                pe1_path, pe2_path ], stdout = bwa_file, stderr = null )
 


    if args.sam:
        shutil.move(bwa_output ,output_path+'.sam')
        #os.rename(bwa_output ,output_path+'.sam')
    else:
        sam_to_bam( bwa_output, bwa_output + ".bam" )
        if args.sort:
            # coordinate sort the file
            pysam.sort( bwa_output + ".bam", output_path )
            pysam.index(output_path+'.bam')
        else:
            shutil.move(bwa_output +".bam",output_path+'.bam')
            #os.rename(bwa_output + ".bam",output_path+'.bam')

    #else:
    #    os.rename(bwa_output + ".bam",output_path+'.bam')

def map_single_reads(pe_path, genome_path, output_path):
    work_dir = tempfile.mkdtemp( )
    genome_db = os.path.join( work_dir, "genome" )
    pe_output = os.path.join( work_dir, "pe.sai" )
    bwa_output = os.path.join( work_dir, "output.sam" )
    
    null = open( "/dev/null" )
    subprocess.check_call( [ "bwa", "index", "-p", genome_db, genome_path ], stderr = null )
    with open( pe_output, "w" ) as pe_file:
        subprocess.check_call( [ "bwa", "aln", genome_db, pe_path ], stdout = pe_file, stderr = null )
    
    with open( bwa_output, "w" ) as bwa_file:
        subprocess.check_call( [ "bwa", "samse",
                                "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                                genome_db,
                                pe_output,
                                pe_path ], stdout = bwa_file, stderr = null )

    shutil.move(bwa_output ,output_path+'.bam')
    #os.popen('mv '+bwa_output + ' '+ output_path)

    #sam_to_bam( bwa_output, bwa_output + ".bam" )
    #pysam.sort( bwa_output + ".bam", output_path )
    #pysam.index(output_path+'.bam')


def bwa_mem(pe1_path, pe2_path, genome_path, output_path, args):
    print 'Aligning with bwa mem'
    start = time()
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    stderr_file = open(output_path+'.bwa.1','w')

    #null = open("/dev/null")
    subprocess.check_call([ "bwa", "index", "-p", genome_db, genome_path ], stderr=stderr_file)
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call([ "bwa", "mem", "-t", args.threads,
                                genome_db, pe1_path, pe2_path ],
                              stdout=bwa_file,
                              stderr=stderr_file)

    elapsed = time() - start
    print 'Time elapsed for bwa mem: ', elapsed

    if args.sam:
        shutil.move(bwa_output ,output_path+'.sam')
        #os.rename(bwa_output ,output_path+'.sam')
    else:
        sam_to_bam( bwa_output, bwa_output + ".bam" )
        if args.sort:
            print 'here'
            # coordinate sort the file
            pysam.sort( bwa_output + ".bam", output_path )
            pysam.index(output_path+'.bam')
        else:
            shutil.move(bwa_output +".bam",output_path+'.bam')
            #os.rename(bwa_output + ".bam",output_path+'.bam')

    # sam_to_bam(bwa_output, bwa_output + ".bam")
    # pysam.sort(bwa_output + ".bam", output_path)
    # pysam.index(output_path + '.bam')

def bowtie(pe1_path, pe2_path, genome_path, output_path, args):
    print 'Aligning with bowtie'
    start = time()
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bowtie_output1 = os.path.join(work_dir, "bowtie1.sam")
    bowtie_output2 = os.path.join(work_dir, "bowtie2.sam")
    stderr_file = open(output_path+'.bowtie.1','w')

    #null = open("/dev/null")
    subprocess.check_call([ "bowtie-build", genome_path, genome_db ], stderr=stderr_file)
    with open(bowtie_output1, "w") as bowtie_file:
        subprocess.check_call([ "bowtie", "-f", "-S",
                                genome_db, pe1_path ],
                              stdout=bowtie_file,
                              stderr=stderr_file)

    with open(bowtie_output2, "w") as bowtie_file:
        subprocess.check_call([ "bowtie", "-f", "-S",
                                genome_db, pe2_path ],
                              stdout=bowtie_file,
                              stderr=stderr_file)

    elapsed = time() - start
    print 'Time elapsed for bowtie: ', elapsed

    if args.sam:
        shutil.move(bowtie_output1 ,output_path+'1.sam')
        #os.rename(bowtie_output1 ,output_path+'1.sam')
        shutil.move(bowtie_output2 ,output_path+'2.sam')
        #os.rename(bowtie_output2 ,output_path+'2.sam')
    else:
        sam_to_bam( bowtie_output1, output_path + "1.bam" )
        sam_to_bam( bowtie_output2, output_path + "2.bam" )
        # if args.sort:
        #     print 'here'
        #     # coordinate sort the file
        #     pysam.sort( bowtie_output + ".bam", output_path )
        #     pysam.index(output_path+'.bam')
        # else:
        #     os.rename(bwa_output + ".bam",output_path+'.bam')



if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Map reads with bwa." )
    parser.add_argument( 'pe1_path', type=str, help='Path to the first reads in a read pair if read pair. Just the reads if single reads' )
    parser.add_argument( 'pe2_path', type=str, nargs = '?', default='', help='Path to the second pairs. Leave unspecified if single reads.' )
    parser.add_argument( 'genome_path', type=str, help='Path to the reference genome/contigs.' )
    parser.add_argument( 'output_path', type=str, help='Output path + filename of resulting .bam file, e.g. /path/to/output/bam_file' )

    parser.add_argument( '-sort', dest='sort', action='store_true', default=False, help='Coordinate sort the reads in the bam file' )
    parser.add_argument( '-sam', dest='sam', action='store_true', default=False, help='Output a samfile (default is bam)' )
    parser.add_argument('--threads', type=str, dest='threads', default='8', required=False, help='Number of threads for bwa mem.')
    parser.add_argument('--mem', dest='mem', action="store_true", required=False,
                        help='Align with bwa mem.')
    parser.add_argument('--bowtie', dest='bowtie', action="store_true", required=False,
                        help='Align reads with bowtie.')

    args = parser.parse_args( )

    if args.mem:
        bwa_mem( args.pe1_path, args.pe2_path, args.genome_path, args.output_path, args )
    elif args.bowtie:
        bowtie( args.pe1_path, args.pe2_path, args.genome_path, args.output_path, args )

    if args.pe2_path:
        map_paired_reads( args.pe1_path, args.pe2_path, args.genome_path, args.output_path, args )
    else:
        map_single_reads(args.pe1_path, args.genome_path, args.output_path )
