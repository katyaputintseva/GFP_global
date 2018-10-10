import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import collections

def load_file(f):
    
    '''Load file using readlines and return the resulting list.'''
    
    with open(f, 'r') as f:
        samples = f.readlines()
    return samples

def parse_fastq(samples):
    
    '''Parse a list returned by the load_file function into a sequence dict and
    a quality dict, where keys are the corresponding barcodes.'''
    
    samples_dict = {}
    quality_dict ={}
    for i in range(0,len(samples),4):
        name = samples[i][1:].rstrip('\n')
        sq = samples[i+1].rstrip('\n')
        quality = samples[i+3].rstrip('\n')
        samples_dict[name] = sq
        quality_dict[name] = quality
    return samples_dict, quality_dict

def revcomp(sq):

    """This function returns reverse compliment of the input sequence."""

    revcomp_sq = ''
    compdict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for nt in sq[::-1]:
        revcomp_sq += compdict[nt]

    return revcomp_sq

def flatten(l):

    """This function flattens a list of lists."""

    return [item for sublist in l for item in sublist]


class Mutant:

    '''This class is for mutants extracted from fastq files.'''

    def __init__(self, sq, quality, bc):

        '''Initializes a mutant with a nt sequence,
        a quality string, and a mutant name.'''

        self.mutations = []
        self.nt = sq
        self.quality = quality
        try:
            self.bc, self.bc_coverage = bc.split(':')
        except:
            self.bc = bc
        self.length = len(self.nt)

    def slicer(self, pattern, before=True):

        """This function slices the input string (my_string) either up until the
        regexp found (before=True), or after (in the latter case, the regexp is
        also cut away from my_str). """

        index = re.search(pattern, self.nt).start()

        if index != -1:

            if before:
                self.nt = self.nt[index:]
                self.quality = self.quality[index:]
            else:
                self.nt = self.nt[:index]
                self.quality = self.quality[:index]

    def find_pattern_and_trim(self, pattern, before):

        """This function searches for a pattern (both forward and reverse complement)
        in sq and cuts it up until this pattern. It cuts the quality sequence accordingly,
        so that only the quality scores that correspond to the obtained trimmed sequence
        are stored. It also performs reverse complement of sequences, which are not
        oriented in the right way.

        In case start pattern was not found, both sequence and quality of the mutant
        become empty strings."""

        if re.search(pattern, self.nt):
            self.slicer(pattern, before)
            
        elif re.search(revcomp(pattern), self.nt):
            self.quality = self.quality[::-1]
            self.nt = revcomp(self.nt)
            self.slicer(pattern)
            
        else:
            self.quality = ''
            self.nt = ''

        self.length = len(self.nt)
        
    def find_barcode_and_trim(self, pattern, before):

        """This function extracts barcodes from sequences, searching for the input pattern and
        adjusting for the indentation of the library design. It crops the target sequence and
        the corresponding quality string up untilthe start of the pattern. It also splits the
        quality str into the barcode and the target sequence parts.

        In case start pattern was not found, both sequence and quality of the mutant
        become empty strings."""

        if re.search(pattern, self.nt):
            index = re.search(pattern, self.nt).start()
            self.slicer(pattern, before)

        else:
            self.quality = ''
            self.nt = ''

        self.length = len(self.nt)

    def extract_mutations(self, ref_seq):

        """Given a non-mutated sequence, this function converts a mutant sequence into
        a list of mutated positions in the following format:

        A10G
        where G = original_nt, 10 = position of mutation (numberation starts from 0),
        G = new_nt"""

        for i, nt in enumerate(ref_seq):
            if nt != self.nt[i]:
                self.mutations.append(nt + str(i) + self.nt[i])

        self.n_mutations = len(self.mutations)
        self.mutations_pos = [int(x[1:-1]) for x in self.mutations]
        self.mutation_type = [x[0] + x[-1] for x in self.mutations]
        self.mutations_quality = [self.quality[int(x[1:-1])] for x in self.mutations]
        
class Library:

    """This class is a collection of Mutants, extracted from a fastq file."""

    def __init__(self, samples_dict, quality_dict):

        """Initializes based on a sequences and quality dictionaries.
        Extracts sequences as instances of the Mutant class."""

        self.sequences = []

        for bc in samples_dict:
            self.sequences.append(Mutant(sq=samples_dict[bc],
                                       quality=quality_dict[bc],
                                       bc=bc))

        print('This library size is %d mutants' % len(self.sequences))

    def trim_start(self, pattern):

        """This function applies find_pattern_and_trim method of the Mutant class
        to all mutants of the library, where start pattern is used."""

        for sq in self.sequences:
            sq.find_pattern_and_trim(pattern, True)
            
    def trim_end(self, pattern):

        """This function applies find_pattern_and_trim method of the Mutant class
        to all mutants of the library, where end pattern is used."""

        for sq in self.sequences:
            sq.find_barcode_and_trim(pattern, False)

    def clean_library(self):

        """This function deletes all the empty mutants from the library."""

        print('Before cleaning we had %d mutants.' % len(self.sequences))

        self.sequences = [sq for sq in self.sequences if sq.length != 0]

        print('Now we have %d mutants left.' % len(self.sequences))
        
    def filter_by_length(self, length):

        """This function filters mutants by the specified length.
        ACHTUNG: filtering is performed inplace."""

        print('\n----\nBefore filtering we had %d mutants.' % len(self.sequences))

        self.sequences = [x for x in self.sequences if x.length == length]

        print('Now we have %d mutants left.' % len(self.sequences))

    def extract_mutations(self, ref_seq):

        """This function applies extract_mutations method of the Mutant class
        to all mutants of the library."""

        for sq in self.sequences:
            sq.extract_mutations(ref_seq)

        print('This library contains %d unique mutations' % len(set(flatten([sq.mutations for sq in self.sequences]))))
        
    def drop_mega_mutants(self, n_mutations):
        
        """Drop mutants with more than n_mutations mutations."""
        
        print('\n----\nBefore filtering we had %d mutants.' % len(self.sequences))

        self.sequences = [x for x in self.sequences if x.n_mutations <= n_mutations]

        print('Now we have %d mutants left.' % len(self.sequences))
        

def import_library(f, start_pattern, end_pattern, sq_length):
    
    '''Import a sequence library from a file.
    Trim the start pattern, trim the end pattern.
    Filter sequences by length. Exctract all the barcodes 
    in the library into a list.'''
    
    print("\nLoading the file")
    lines = load_file(f)
    samples_dict, quality_dict = parse_fastq(lines)
    
    print("\nParsing the library")
    library = Library(samples_dict, quality_dict)
    
    print("\nTrimming sequences from the beginning")
    library.trim_start(start_pattern)
    
    print("\nTrimming sequences from the end")
    library.trim_end(end_pattern)
    
    print("\nFiltering sequences by length")
    library.filter_by_length(sq_length)
        
    return library
       
        