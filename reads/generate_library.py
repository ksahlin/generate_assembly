import sys
import os
import random
import operator

##
# General functions

def reverse_complement(string):
    """
        Reverse complements a DNA string
        Arguments:
        string  - A DNA string

        Returns:
            A python string object that represents the
            reverse complement of the input (DNA) string.

    """
    rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def read_in_fasta_file(fasta_file):
    """
        Reads a fasta file into memory.

        Arguments:
        fasta_file - A python file object. The file should be in 
        fasta format.

        Returns:
            A python dictionary with accessions as keys and
            sequences as values.

    """  
    seqs = {}
    k = 0
    temp = []
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            temp = ''.join(temp)
            seqs[accession] = temp
            temp = []
            accession = line[1:].strip().split()[0]
        else:
            temp.append(line.strip())
    
    
    temp = ''.join(temp)
    seqs[accession] = temp
    return(seqs)


def generate_contigs(genome,min_c_len,max_c_len,min_distance,max_distance):
    position = 0
    index = 0
    while True:
        contig_len = random.randrange(min_c_len,max_c_len)
        position += random.randrange(min_distance,max_distance)
        if position + contig_len > len(genome):
            break

        rev_comp = random.randrange(0,2)
        if rev_comp:
            yield '>c{0},pos:{1}-{2},rc:1\n{3}\n'.format(index,
                position,position+contig_len, reverse_complement( genome[position:position+contig_len]))
        else:
            yield '>c{0},pos:{1}-{2},rc:0\n{3}\n'.format(index,
                position,position+contig_len, genome[position:position+contig_len])
        index += 1
        position += contig_len

def generate_contigs_two_sizes(genome,small_size,large_size,min_distance,max_distance,distr_weight_large_ctgs):
    position = 0
    index = 0
    while True:
        r = random.uniform(0,1)
        if r < distr_weight_large_ctgs:
            contig_len = large_size
        else:
            contig_len = small_size

        position += random.randrange(min_distance,max_distance)
        if position + contig_len > len(genome):
            break

        rev_comp = random.randrange(0,2)
        if rev_comp:
            yield '>c{0},pos:{1}-{2},rc:1\n{3}\n'.format(index,
                position,position+contig_len, reverse_complement( genome[position:position+contig_len]))
        else:
            yield '>c{0},pos:{1}-{2},rc:0\n{3}\n'.format(index,
                position,position+contig_len, genome[position:position+contig_len])
        index += 1
        position += contig_len



def insertion():
    return(''.join( random.choice('AGCT') for i in range( random.randint(1,10) ) ))

def deletion():
    return( abs(int(random.gauss(4,2))))

def mutation():
    return(random.choice('AGCT'))
        


# class PairedEndRead(object):
#     """docstring for PairedEndRead"""
#     def __init__(self):
#         super(PairedEndRead, self).__init__()

#     def generate(self,reference_accession, reference_sequence, read_index,mean,sigma,read_length):
#         self.fragment_length = int(random.gauss(mean,sigma))
#         if self.fragment_length >= len(reference_sequence): 
#             raise Exception("To short reference sequence length for \
#                 simulated read. \nRead fragment: {0}\nTranscript \
#                 length:{1}".format(self.fragment_length,len(reference_sequence)))
#         self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
#         self.read1 = reference_sequence[self.start_pos : self.start_pos + read_length]
#         self.read2 = reverse_complement(reference_sequence[self.start_pos + self.fragment_length - read_length : self.start_pos+self.fragment_length])
#         self.reference_accession = reference_accession
#         self.read_index = read_index
#         self.read_length = read_length
#     def fastq_format(self):
#         r1= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
#         \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
#             self.reference_accession,self.read1,'J'*self.read_length)
#         r2= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
#         \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
#             self.reference_accession,self.read2,'J'*self.read_length)        
#         yield r1
#         yield r2

#     def fasta_format(self):
#         r1= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
#             self.reference_accession,self.read1)
#         r2= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
#             self.reference_accession,self.read2)        
#         yield r1
#         yield r2

class MatePairRead(object):
    """docstring for MatePairRead"""
    def __init__(self):
        super(MatePairRead, self).__init__()

    def generate(self,reference_accession, reference_sequence, read_index,mean,sigma,read_length):
        self.fragment_length = int(random.gauss(mean,sigma))
        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reverse_complement(reference_sequence[self.start_pos : self.start_pos + read_length])
        self.read2 = reference_sequence[self.start_pos + self.fragment_length - read_length : self.start_pos+self.fragment_length]
        self.reference_accession = reference_accession
        self.read_index = read_index
        self.read_length = read_length
    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2




class Transcriptome(object):
    """docstring for Transcriptome"""
    def __init__(self):
        super(Transcriptome, self).__init__()

    def generate_splice_variants(self, nr_transcripts,genome,
        max_intron_size,min_intron_size,max_exon_size, min_exon_size):
        self.transcripts = []
        for i in range(nr_transcripts):
            self.transcripts.append(Transcript(genome,max_intron_size,min_intron_size,max_exon_size, min_exon_size))

        return self.transcripts   

    def fasta_format(self):
        for transcript in self.transcripts:
            yield str(transcript)


class Transcript(object):
    """
    Class for an Transcript object
    
    Attributes:

    genome              - A python string object that is the sequence of
                        the genome.
    max_intron_size     - Maximum size of intron
    min_intron_size     - Minimun size of intron
    max_exon_size       - Maximum size of exon
    min_exon_size       - Minimun size of exon    
    gene_distance       - Distance between genes

    """
    def __init__(self, genome_strand,max_intron_size,min_intron_size,max_exon_size, min_exon_size):
        super(Transcript, self).__init__()
        self.genome_strand = genome_strand
        self.max_intron_size = max_intron_size
        self.min_intron_size = min_intron_size
        self.max_exon_size = max_exon_size
        self.min_exon_size = min_exon_size

        self.get_sequence()
        


    def get_sequence(self):
        """
        Generates an Transcript from a genome
        Returns:
        An exon
        """

        nr_exons = random.randrange(1,5)
        self.intron_length = [random.randrange(self.min_intron_size,self.max_intron_size) for i in range(nr_exons)]
        self.exon_lengths = [random.randrange(self.min_exon_size,self.max_exon_size) for i in range(nr_exons)]
        self.reverse_complement = random.randrange(2)

        if sum(self.intron_length) + sum(self.exon_lengths) >= len(self.genome_strand):
            self.get_sequence()

        if self.reverse_complement:
            self.start_position = random.randrange(sum(self.intron_length) + sum(self.exon_lengths),len(self.genome_strand))
        else:
            self.start_position = random.randrange(0,len(self.genome_strand)-sum(self.intron_length)-sum(self.exon_lengths))


        position = self.start_position
        self.positions = []
        self.sequence = ''
        for e_len,i_len in zip(self.exon_lengths,self.intron_length):
            if self.reverse_complement:
                self.sequence += reverse_complement( self.genome_strand[position-e_len:position])
                self.positions.append((position-e_len,position))
                position -= e_len + i_len

            else:

                self.sequence += self.genome_strand[position:position+e_len]
                self.positions.append((position,position+e_len))
                position += e_len + i_len

        self.accession = '>spliced_variant{0},rc={1}'.format(self.positions,self.reverse_complement)


    def __str__(self):
        return '{0}\n{1}\n'.format(self.accession, self.sequence)
        




################################

# import os
# import subprocess
# import tempfile
# import fnmatch
# import shutil
import random

from numpy.random import lognormal
# from pkg_resources import resource_filename


##
# Simulates reads using a lognormal distribution.
#
class LogNSimulator(object):
    ##
    # Constructor.
    #
    def __init__(self, **kwargs):
        super(LogNSimulator, self).__init__()
        for key, value in kwargs.iteritems():
            print(key, value)
            setattr(self, key, value)       

    def simulate(self, genome_file, output_prefix, second_genome=None):
        if not self.read_length == 100:
            raise ValueError("Read length must be 100 for metasim.")
        if self.mean < 10 and  0 < self.std <  2: # Well defined distribution
            pass
        else: # user probably forgot that mu and std needs to be specified in log base
            raise ValueError("mu and std needs to be specified in log base (usually mu < 10 and 0 < sigma < 2).\
                    You specified mean: {0}, sigma: {1}".format(self.mean, self.std))

        all_sequences_dict = read_in_fasta_file(genome_file)
        if len(all_sequences_dict) != 1:
            print all_sequences_dict
            print len(all_sequences_dict)
            raise NotImplementedError("lognsim is only implemented for one genome sequence in the reference genome file.")

        genome_sequence = all_sequences_dict.values()[0]
        genome_length = len(genome_sequence)
        dna_library = DNAseq(self.read_length, self.coverage, mean=self.mean, stddev=self.std, distribution='lognormal')
        dna_library.simulate_pe_reads(genome_sequence)

        reads1 = open(output_prefix + "_pe1.fa",'w')
        reads2 = open(output_prefix + "_pe2.fa",'w')
        i=0
        for read in dna_library.fasta_format():
            if i%2==0:
                reads1.write(read)
            else:
                reads2.write(read)
            i+=1

        # output_dir = tempfile.mkdtemp( )

        # subprocess.call( [ "MetaSim", "cmd",
        #                    "-r", str( num_reads ), 
        #                    "-m",
        #                    "-g", self.error_model,
        #                    "-2", self.error_model,
        #                    "--empirical-pe-probability", "100",
        #                    "--clones-mean", str( self.mean ),
        #                    "--clones-param2", str( self.std ),
        #                    "-d", output_dir,
        #                    genome_path ], stdout = open( "/dev/null", "w" ) )


class PairedEndRead(object):
    """docstring for PairedEndRead"""
    def __init__(self,distribution = 'normal',mean=None,sigma=None,read_length= None, min_size=None, max_size=None):
        super(PairedEndRead, self).__init__()
        self.distribution = distribution
        self.mean = mean
        self.sigma = sigma
        self.read_length = read_length
        self.min_size = min_size
        self.max_size = max_size

    def generate(self, reference_sequence, read_index):
        if self.distribution == 'normal':
            self.fragment_length = int(random.gauss(self.mean,self.sigma))
        elif self.distribution == 'uniform':
            self.fragment_length = int(random.uniform(self.min_size,self.max_size))
        elif self.distribution == 'lognormal':
            self.fragment_length = int(lognormal(self.mean, self.sigma)) # one sample at a time to conform with the implementation...

        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reference_sequence[self.start_pos : self.start_pos + self.read_length]
        self.read2 = reverse_complement(reference_sequence[self.start_pos + self.fragment_length - self.read_length : self.start_pos+self.fragment_length])
        self.reference_accession = "reference_genome"
        self.read_index = read_index

    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2

def calculate_num_reads(coverage, read_length, genome_length):
    """
    Calculate the number of reads required to get a certain coverage.
    
    Arguments:
        coverage (int): Desired mean genome coverage.
        read_length (int): Length of each read in a pair.
        genome_length (int): Estimated genome length.
    
    Returns:
        Number of reads need to get the given coverage.
    """
    return int((genome_length * coverage) / (read_length * 2.0))

class DNAseq(object):
    """docstring for DNAseq"""
    def __init__(self,read_length, coverage, mean=None,stddev=None, min_size=None, max_size = None, distribution='normal'):
        super(DNAseq, self).__init__()
        self.distribution = distribution
        
        if self.distribution == 'normal' or self.distribution == "lognormal":
            self.mean = mean
            self.stddev = stddev
        elif self.distribution == 'uniform':
            self.min_size = min_size
            self.max_size = max_size

        self.read_length = read_length
        self.coverage = coverage

    def simulate_pe_reads(self, genome_sequence):
        """
        Arguments:
        """
        genome_length = len(genome_sequence)
        number_of_reads = calculate_num_reads(self.coverage, self.read_length, genome_length)  # Specifies the number of simulated read pairs (related to insertion size length of genome and coverage
    
        self.reads = []
        
        
        for i in range(number_of_reads):
            if self.distribution == 'normal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'lognormal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'uniform':
                read_pair = PairedEndRead(distribution=self.distribution, min_size=self.min_size,max_size=self.max_size,read_length=self.read_length)
     
            read_pair.generate(genome_sequence, i)
            self.reads.append(read_pair)


    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate





# class DNAseq(object):
#     """docstring for DNAseq"""
#     def __init__(self,lib_read_length,coverage, lib_mean,lib_std_dev,contamination_rate = 0,contamine_mean= None,contamine_std_dev= None):
#         super(DNAseq, self).__init__()

#         self.lib_mean = lib_mean
#         self.lib_std_dev = lib_std_dev
#         self.lib_read_length = lib_read_length
#         self.contamination_rate = contamination_rate
#         self.contamine_mean = contamine_mean
#         self.contamine_std_dev = contamine_std_dev
#         self.coverage = coverage
    
#     def simulate_pe_reads(self,genome):
#         """
#         Arguments:

#         """
#         genome_length = len(genome.sequence)
#         number_of_reads=(genome_length*self.coverage)/(2*self.lib_read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
    
#         self.reads = []
        
        
#         for i in range(number_of_reads):
#             read_pair = PairedEndRead()
#             read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i, self.lib_mean,self.lib_std_dev,self.lib_read_length)        
#             self.reads.append(read_pair)

#     def simulate_mp_reads(self,genome):
#         """
#         Arguments:

#         """
#         genome_length = len(genome.sequence)
#         number_of_reads=(genome_length*self.coverage)/(2*self.lib_read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
    
#         self.reads = []
        
        
#         for i in range(number_of_reads):
#             r = random.uniform(0,1)
#             if r < self.contamination_rate:
#                 read_pair = PairedEndRead()
#                 read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i, self.contamine_mean,self.contamine_std_dev,self.lib_read_length)

#             else:
#                 read_pair = MatePairRead()
#                 read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i, self.lib_mean,self.lib_std_dev,self.lib_read_length)
            
#             self.reads.append(read_pair)


#     def fastq_format(self):
#         for pe_read in self.reads:
#             for mate in pe_read.fastq_format():
#                 yield mate

#     def fasta_format(self):
#         for pe_read in self.reads:
#             for mate in pe_read.fasta_format():
#                 yield mate



class RNAseq(object):
    """docstring for RNAseq"""
    def __init__(self,lib_mean,lib_std_dev,lib_read_length):
        super(RNAseq, self).__init__()

        self.lib_mean = lib_mean
        self.lib_std_dev = lib_std_dev
        self.lib_read_length = lib_read_length

    def simulate_pe_reads(self,transcriptome,coverage):
        """
        Arguments:

        transcriptome   - A Transcriptome object
        """
        transcriptome_length = reduce(lambda x,y: x+y, [len(transcript.sequence) for transcript in transcriptome.transcripts]) # total transcriptome length
        number_of_reads=(transcriptome_length*coverage)/(2*self.lib_read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
        
        self.reads = []
        
        
        for i in range(number_of_reads):
            transcript_index = int(random.randrange(len(transcriptome.transcripts)))
            transcript = transcriptome.transcripts[transcript_index]
            pe = PairedEndRead()
            pe.generate(transcript.accession[1:].replace(" ",""), transcript.sequence, i, self.lib_mean,self.lib_std_dev,self.lib_read_length)
            self.reads.append(pe)

    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate




if __name__=='__main__':
    pass