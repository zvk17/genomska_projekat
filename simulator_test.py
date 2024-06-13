import os
import unittest
from tempfile import NamedTemporaryFile

from SequenceRead import SequenceRead

from simulator import (
    createCigarStringFromCigarArray,
    getSamFileChromosomeName,
    createAlignedSegmentFromRead, writeReadsToFqFile
)


class CigarGenerationCase(unittest.TestCase):
    def test_cigar_generation_1(self):
        input = ['M','M','M']
        output = "3M"
        self.assertEqual(createCigarStringFromCigarArray(input), output)  

    def test_cigar_generation_2(self):
        input = ['M', 'I', 'I', 'M']
        output = "1M2I1M"
        self.assertEqual(createCigarStringFromCigarArray(input), output)  

    def test_cigar_generation_3(self):
        input = ['M', 'D', 'I', 'M']
        output = "1M1D1I1M"
        self.assertEqual(createCigarStringFromCigarArray(input), output)  

class ChromosomeNameGenerationTestCase(unittest.TestCase):
    def test_chromosome_name_generation_1(self):
        input = ">NC_039898.1 Coffea arabica cultivar Caturra red isolate CCC135-36 chromosome 1c, Cara_1.0, whole genome shotgun sequence"
        output = "NC_039898.1"
        self.assertEqual(getSamFileChromosomeName(input), output)
    def test_chromosome_name_generation_2(self):
        input = "NC_039898.1 Coffea arabica cultivar Caturra red isolate CCC135-36 chromosome 1c, Cara_1.0, whole genome shotgun sequence"
        output = "chr_unknown"
        self.assertEqual(getSamFileChromosomeName(input), output)
    def test_chromosome_name_generation_3(self):
        input = ">NC_030384.2 Daucus carota subsp. sativus chromosome 4, DH1 v3.0, whole genome shotgun sequence"
        output = "NC_030384.2"
        self.assertEqual(getSamFileChromosomeName(input), output)




class AlignedSegmentCreationTestCase(unittest.TestCase):
    def testAlignedSegmentCreation(self):
        sequenceId = "my_id"
        cigarstring = "3M"
        read2 = True
        nucleotides = "AAA"
        input = SequenceRead(
            id = sequenceId,
            cigar=cigarstring,
            quality=[20.0,20.0,20.0],
            is_second_read=read2,
            snvCount=0,
            nucleotides=nucleotides,
            position=100
        )
        result = createAlignedSegmentFromRead(input)
        self.assertEqual(result.query_name, input.idName())
        self.assertEqual(result.cigarstring, cigarstring)
        self.assertEqual(result.is_read2, read2)
        self.assertEqual(result.query_sequence, nucleotides)

class WritingFilesTestCase(unittest.TestCase):
    FQ_LINES_NUMBER_PER_READ = 4

    def testWriting(self):
        sequenceId = "my_id"
        cigarstring = "3M"
        read2 = True
        nucleotides = "AAA"

        sequence_read = SequenceRead(
            id=sequenceId,
            cigar=cigarstring,
            quality=[20.0, 20.0, 20.0],
            is_second_read=read2,
            snvCount=0,
            nucleotides=nucleotides,
            position=100
        )


        readList = [sequence_read, sequence_read]

        i = 0
        try:
            fqFile = NamedTemporaryFile(mode="w", delete=False)
            fqFile.close()
            writeReadsToFqFile(fqFile.name, readList)
            with open(fqFile.name, "r") as readFqFile:
                i = 0
                for line in readFqFile.readlines():
                    i = i + 1
        finally:
            os.remove(fqFile.name)
        self.assertEqual(i, self.FQ_LINES_NUMBER_PER_READ * len(readList))




if __name__ == '__main__':
    unittest.main()
