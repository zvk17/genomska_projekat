from RandomBuffer import RandomBuffer
import re

class Configuration:
    MIN_NUCLEOTIDE_QUALITY = 0.0
    MAX_NUCLEOTIDE_QUALITY = 40.0
    RANDOM_SEED = 42
    READ_FILTER_PATTERN = re.compile("^[ATCGU]+$")

    singleNucleotideVariant = {
        'A': ['T', 'C', 'G', ],
        'T': ['A', 'C', 'G', ],
        'C': ['T', 'A', 'G', ],
        'G': ['T', 'C', 'A', ],

        'U': ['A', 'C', 'G', ],  #  rna
    }
    insertionChoice = ['A', 'C', 'G', 'T']

    BUFFER_SIZE = 1000 * 100

    def __init__(self):
        # sort choices
        for key in self.singleNucleotideVariant.keys():
            self.singleNucleotideVariant[key].sort()

