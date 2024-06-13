import math
import numpy


class SequenceRead:
    def __init__(self, id, position, is_second_read, nucleotides, quality, cigar, snvCount):
        self.id = id
        self.position = position
        self.is_reverse = is_second_read
        self.nucleotides = nucleotides
        self.quality = [int(q + 0.5) for q in quality]
        self.cigar = cigar
        self.snvCount = snvCount
    def getNucleotidesString(self):
        return "".join(self.nucleotides)
    def qualityToString(self) -> str:
        qualityString = ""
        for qval in self.quality:
            asciiChar = chr(ord('!') + qval)
            qualityString = qualityString + asciiChar
        return  qualityString
    def idName(self):
        return "ID:" + str(self.id)

    def writeCigarToFile(self, file):
        file.write("@ID:")
        file.write(str(self.id))
        file.write(" ")
        file.write(self.cigar)
        file.write(" ")
        file.write(str(self.snvCount))
        file.write("\n")

    def writeToFile(self, file):
        file.write("@ID:")
        file.write(str(self.id))
        file.write("\n")
        file.write("".join(self.nucleotides))
        file.write("\n")
        file.write("+\n")
        file.write(self.qualityToString())
        file.write("\n")
