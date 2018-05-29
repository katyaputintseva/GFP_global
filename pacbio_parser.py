import re


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


class Read:

    '''This class is for reads extracted from fastq files.'''

    def __init__(self, sq, quality, name):

        '''Initializes a read with a nt sequence,
        a quality string, and a read name.'''

        self.mutations = []
        self.nt = sq
        self.quality = quality
        self.name = name
        self.length = len(self.nt)

    def slicer(self, pattern, before=True):

        """This function slices the input string (my_string) either up until the
        regexp found (before=True), or after (in the latter case, the regexp is
        also cut away from my_str). The function also returns index of regexp
        found -- in order to be able to use it with other strings, f.e. quality."""

        index = re.search(pattern, self.nt).start()

        if index != -1:

            if before:
                self.nt = self.nt[index:]
                self.quality = self.quality[index:]
            else:
                self.nt = self.nt[:index]
                self.quality = self.quality[:index]

    def find_start_and_trim(self, pattern):

        """This function searches for pattern (both forward and reverse complement)
        in sq and cuts it up until this pattern. It cuts the quality sequence accordingly,
        so that only the quality scores that correspond to the obtained trimmed sequence
        are stored. It also performs reverse complement of sequences, which are not
        oriented in the right way.

        In case start pattern was not found, both sequence and quality of the read
        become empty strings."""

        if re.search(pattern, self.nt):
            self.slicer(pattern)

        elif re.search(revcomp(pattern), self.nt):
            self.quality = self.quality[::-1]
            self.nt = revcomp(self.nt)
            self.slicer(pattern)

        else:
            self.quality = ''
            self.nt = ''

        self.length = len(self.nt)

    def find_barcode(self, pattern, indentation, bc_length):

        """This function extracts barcodes from sequences, searching for the input pattern and
        adjusting for the indentation of the library design. It crops the target sequence and
        the corresponding quality string up untilthe start of the pattern. It also splits the
        quality str into the barcode and the target sequence parts.

        The output is: bc sequence, bc quality, target sq, target sq quality.

        In case start pattern was not found, both sequence and quality of the read
        become empty strings."""

        if re.search(pattern, self.nt):
            self.bc = re.search(pattern, self.nt).group()[indentation:]
            index = re.search(pattern, self.nt).start()
            self.bc_quality = self.quality[index + indentation:index + bc_length + indentation]
            self.slicer(pattern, before=False)

        else:
            self.quality = ''
            self.nt = ''

        self.length = len(self.nt)

    def extract_mutations(self, ref_seq):

        """Given a non-mutated sequence, this function converts a read sequence into
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

    """This class is a collection of Reads, extracted from a fastq file."""

    def __init__(self, file_path):

        """Initializes based on a fastq file. Extracts reads as instances of the
        Read class."""

        with open(file_path, 'r') as f:
            lines = f.readlines()

        self.sequences = []

        for i in range(0, len(lines), 4):
            self.sequences.append(Read(sq=lines[i + 1].rstrip('\n'),
                                       quality=lines[i + 3].rstrip('\n'),
                                       name=lines[i].rstrip('\n')))

        print('This library size is %d consensus reads' % len(self.sequences))

    def filter_by_length(self, length=733):

        """This function filters reads by the specified length.
        ACHTUNG: filtering is performed inplace."""

        print('Before filtering we had %d reads.' % len(self.sequences))

        self.sequences = [x for x in self.sequences if x.length == length]

        print('Now we have %d reads left.' % len(self.sequences))

    def trim_start(self, pattern):

        """This function applies find_start_and_trim method of the Read class
        to all reads of the library."""

        for sq in self.sequences:
            sq.find_start_and_trim(pattern)

    def get_bc_and_trim(self, pattern, indentation, bc_length):

        """This function applies find_barcode method of the Read class
        to all reads of the library."""

        for sq in self.sequences:
            sq.find_barcode(pattern, indentation, bc_length)

    def clean_library(self):

        """This function deletes all the empty reads from the library."""

        print('Before cleaning we had %d reads.' % len(self.sequences))

        for i, sq in enumerate(self.sequences):
            if sq.length == 0:
                self.sequences.pop(i)

        print('Now we have %d reads left.' % len(self.sequences))

    def extract_mutations(self, ref_seq):

        """This function applies extract_mutations method of the Read class
        to all reads of the library."""

        for sq in self.sequences:
            sq.extract_mutations(ref_seq)
