import numpy
import pysam
import scipy
from timeit import default_timer as timer
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
from pysam import AlignmentFile, AlignedRead
import SequenceRead
import RandomBuffer
import re
from ArgParser import getArgParser
from configuration import Configuration




insertSize = 1000
readLength = 300





def createSNVChoiceBufferDict(lConfig: Configuration, lRandomNumberGenerator) -> dict:
    lChoiceBufferDict: dict = {}
    for key in lConfig.singleNucleotideVariant.keys():
        lChoiceBufferDict[key] = RandomBuffer.RandomChoiceBuffer(
            randomNumberGenerator= lRandomNumberGenerator,
            bufferSize= lConfig.BUFFER_SIZE,
            choiceArray=lConfig.singleNucleotideVariant[key]
        )
    return lChoiceBufferDict

def readFasta(lFilePath: str) -> tuple[str,str]:
    lFasta: str = ""
    lineCount: int = 0
    headerLine: str = ""
    with open(lFilePath, "r") as lFastaFile:
        #line = file.readline()
        for line in lFastaFile.readlines():
            lineCount = lineCount + 1
            if (line.startswith(">")):
                if lineCount == 1:
                    headerLine = line
                continue
            else:
                lFasta += line.strip()

    return lFasta, headerLine



def truncatedNormalDistributionGenerator(
        lConfig: Configuration,
        lAverageNucleotideQuality: float,
        lStdev: float
):
    lATransformed = (lConfig.MIN_NUCLEOTIDE_QUALITY - lAverageNucleotideQuality) / lStdev
    lBTransformed = (lConfig.MAX_NUCLEOTIDE_QUALITY - lAverageNucleotideQuality) / lStdev

    lRv = scipy.stats.truncnorm(
        loc = lAverageNucleotideQuality,
        a = lATransformed,
        b = lBTransformed,
        scale = lStdev
    )
    return lRv



def createCigarStringFromCigarArray(lCigarArray: str | list) -> str:
    lCigarCondensedArray = ""

    lPreviousCount = 0
    lPreviousType = ''

    for lCurrentType in lCigarArray:

        if lCurrentType != lPreviousType:
            if lPreviousCount != 0:
                lCigarCondensedArray += str(lPreviousCount) + lPreviousType
            lPreviousType = lCurrentType
            lPreviousCount = 1
        else:
            lPreviousCount = lPreviousCount + 1

    if lPreviousCount > 0:
        lCigarCondensedArray += str(lPreviousCount) + lPreviousType
    return lCigarCondensedArray


def generateMutationErrors(nucleotidesList:str|list, buffers: dict, lArgs):
    nucleotideLength =  len(nucleotidesList)

    singleVariantMutations = buffers["randomBuffer"].take(size= nucleotideLength) 
    insertionMutation = buffers["randomBuffer"].take(size= nucleotideLength)      
    deletionMutation = buffers["randomBuffer"].take(size= nucleotideLength)       

    hasSingleVariantMutation = singleVariantMutations <= lArgs.snv_error_rate
    hasInsertionMutation = insertionMutation <= lArgs.insertion_error_rate
    hasDeletionMutation = deletionMutation <= lArgs.deletion_error_rate

    nucleotidesWithMutation = []

    cigarArray = []


    for i in range(nucleotideLength):

        if hasSingleVariantMutation[i]:
            nucleotide = nucleotidesList[i]
            mutation = buffers["choiceBufferDict"][nucleotide].take(1)[0]

            nucleotidesWithMutation.append(mutation)
            cigarArray.append('M')


        elif hasInsertionMutation[i]:
            newNucleotide = buffers["insertionChoiceBuffer"].take(1)[0] 

            nucleotidesWithMutation.append(newNucleotide)
            nucleotidesWithMutation.append(nucleotidesList[i])

            cigarArray.append('I')
            cigarArray.append('M')

        elif hasDeletionMutation[i]:
           cigarArray.append('D')
           continue
        else:
            cigarArray.append('M')
            nucleotidesWithMutation.append(nucleotidesList[i])

    cigarString = createCigarStringFromCigarArray(lCigarArray=(cigarArray))



    return nucleotidesWithMutation, hasSingleVariantMutation, hasInsertionMutation, hasDeletionMutation, cigarString



def generatePairedReads(
        nucleotideRead: str,
        readLength: int, startPos: int,
        buffers:dict,
        lArgs
):

    assert(len(nucleotideRead) > 2 * readLength)


    frontNucleotides = nucleotideRead[0:(readLength)]
    backNucleotides = (nucleotideRead[-readLength:])#[::-1]

    vf = generateMutationErrors(frontNucleotides, buffers, lArgs)
    vb = generateMutationErrors(backNucleotides, buffers, lArgs)
    finalFrontNucleotides = vf[0]
    finalBackNucleotides = vb[0]

    mutationCountFront = vf[1]
    mutationCountBack = vb[1]


    cigarFront = vf[4]
    cigarBack = vb[4]


    frontQuality = buffers["nucleotideQualityBuffer"].take(size=len(finalFrontNucleotides))#rv.rvs
    backQuality = buffers["nucleotideQualityBuffer"].take(size=len(finalBackNucleotides))#rv.rvs

    frontRead = SequenceRead.SequenceRead(
        id = startPos,
        position = startPos,
        is_second_read= False,
        nucleotides = finalFrontNucleotides,
        quality = frontQuality,
        cigar = cigarFront,
        snvCount = mutationCountFront.sum()
    )
    backRead = SequenceRead.SequenceRead(
        id=startPos,
        position=insertSize + startPos - readLength,
        is_second_read=True,
        nucleotides=finalBackNucleotides,
        quality=backQuality,
        cigar=cigarBack,
        snvCount=mutationCountBack.sum()
    )
    return frontRead, backRead


def debugGenome(
        lNucleotideGenome: str|list,
        lReadStartPositions: numpy.array,
        lNumberOfReads: int
):
    print("Genome size: " + str(len(lNucleotideGenome)))
    print("Min max randint " +
          str(lReadStartPositions.min()) +
          " " +
          str(lReadStartPositions.max()))
    print("Number of reads " + str(lNumberOfReads))





def filterBadReads(
        lReadStartPositions: list,
        lNucleotideGenome: list | str,
        lInsertSize: int,
        lConfig: Configuration
):
    print("Filtering bad reads")
    lCleanReadStartPositions = []

    for i in range(len(lReadStartPositions)):
        lStartPos = lReadStartPositions[i]
        nucleotides = lNucleotideGenome[lStartPos:(lStartPos + lInsertSize)]
        if (bool(lConfig.READ_FILTER_PATTERN.match(nucleotides))):
            lCleanReadStartPositions.append(lStartPos)
    return lCleanReadStartPositions




def generateReads(
        lCleanReadStartPositions: list[int],
        lNucleotideGenome: str | list,
        buffers: dict,
        lArgs
):
    lFrontList = []
    lBackList = []


    for lStartPos in lCleanReadStartPositions:
        lNucleotides = lNucleotideGenome[lStartPos:(lStartPos + insertSize)]

        read, reverseRead = generatePairedReads(
            lNucleotides,
            readLength,
            lStartPos,
            buffers,
            lArgs= lArgs
        )

        lFrontList.append(read)
        lBackList.append(reverseRead)


    return lFrontList, lBackList


def measureReadGenerationSpeed(lTimes: list):
    times = numpy.array(lTimes)

    print("duration total ", times.sum())
    print("duration  mean ", times.mean())
    print("duration  max ", times.max())
    print("duration  min ", times.min())

    with open("./output_files/times.txt", "w") as measureTimesFile:
        for t in times:
            measureTimesFile.write(str(t))
            measureTimesFile.write("\n")



def writeReadsToFqFile(lFilePath: str, lReadList: list[SequenceRead]):
    with open(lFilePath, "w") as lWrFile:
        for r in lReadList:
            r.writeToFile(lWrFile)
        lWrFile.flush()
        print("<{}> file closed".format(lFilePath))



def writeDebugCigar(lFrontList):
    with open("./output_files/cigar.txt", "w") as lWriteFile:
        for r in lFrontList:
            r.writeCigarToFile(lWriteFile)
        lWriteFile.flush()
        print("cigar.txt file closed")


def getSamFileChromosomeName(lHeaderLine: str) -> str:
    lSequenceName = "chr_unknown"

    if (len(lHeaderLine) > 0 and lHeaderLine.startswith(">") and len(lHeaderLine.split("-")) > 0):
        lSequenceName = lHeaderLine.split(" ")[0].strip()[1:]

    return lSequenceName

def getSamFileHeader(lNucleotideGenomeSize: int, lSequenceName: str):
    lHeader = {
        'HD': {'VN': '1.0'},
        'SQ': [{
            'LN': lNucleotideGenomeSize,
            'SN': lSequenceName
        }]
    }
    return lHeader

def createAlignedSegmentFromRead(lRead: SequenceRead)-> pysam.AlignedSegment:
    lNewSegment = pysam.AlignedSegment()


    lNewSegment.query_name = lRead.idName()
    lNewSegment.query_sequence = lRead.getNucleotidesString()
    lNewSegment.flag = 0

    lNewSegment.reference_start = lRead.position
    lNewSegment.is_paired = True
    lNewSegment.is_read1 = not lRead.is_reverse
    lNewSegment.is_read2 = lRead.is_reverse
    lNewSegment.is_mapped = True
    lNewSegment.cigarstring = lRead.cigar


    lNewSegment.query_qualities = lRead.quality
    return lNewSegment


def writeSamFile(
        lFilePath: str,
        lFrontReadList: list[SequenceRead],
        lBackReadList: list[SequenceRead],
        lNucleotideGenome: list|str,
        lHeaderLine: str,
):
    lNucleotideGenomeSize = len(lNucleotideGenome)
    lChromosomeName = getSamFileChromosomeName(lHeaderLine)

    with (
        AlignmentFile(
            lFilePath,
            "w",
            header= getSamFileHeader(lNucleotideGenomeSize, lChromosomeName)
    ) as outputSamFile):

        lTotalList = []
        assert(len(lFrontReadList) == len(lBackReadList))

        for i in range(len(lFrontReadList)):
            lTotalList.append(lFrontReadList[i])
            lTotalList.append(lBackReadList[i])

        for read in lTotalList:
            lNewSegment = createAlignedSegmentFromRead(read)
            outputSamFile.write(lNewSegment)
    print("<{}> file closed".format(lFilePath))

def generateRandomBuffers(
        lRandomNumberGenerator,
        lTruncatedNormalRandomGenerator,
        lConfig:Configuration
):
    randomBuffer = RandomBuffer.RandomBuffer(
        bufferSize=lConfig.BUFFER_SIZE * 10,
        randomNumberGenerator=lRandomNumberGenerator
    )
    insertionChoiceBuffer = RandomBuffer.RandomChoiceBuffer(
        randomNumberGenerator=lRandomNumberGenerator,
        bufferSize=lConfig.BUFFER_SIZE,
        choiceArray=lConfig.insertionChoice
    )
    nucleotideQualityBuffer = RandomBuffer.RandomNucleotideQualityBuffer(
        randomNumberGenerator=lTruncatedNormalRandomGenerator,
        bufferSize= lConfig.BUFFER_SIZE * 10,
    )
    choiceBufferDict = createSNVChoiceBufferDict(
        lConfig=lConfig,
        lRandomNumberGenerator=lRandomNumberGenerator
    )
    buffers = {
        "insertionChoiceBuffer": insertionChoiceBuffer,
        "randomBuffer": randomBuffer,
        "nucleotideQualityBuffer": nucleotideQualityBuffer,
        "choiceBufferDict": choiceBufferDict
    }
    return buffers

def main():

    argParser = getArgParser()
    args = argParser.parse_args()

    config = Configuration()
    randomNumberGenerator = numpy.random.default_rng(seed=config.RANDOM_SEED)

    nucleotideGenome, headerLine = readFasta(args.input_file)
    genomeSize = len(nucleotideGenome)

    numberOfReads = int(int((genomeSize / insertSize + 0.5) * args.coverage))



    readStartPositions = randomNumberGenerator.integers(
        low=0,
        high=int(genomeSize - insertSize),
        size=int(numberOfReads)
    )

    debugGenome(
        lNucleotideGenome=nucleotideGenome,
        lReadStartPositions=readStartPositions,
        lNumberOfReads=numberOfReads)

    truncatedNormalRandomGenerator = truncatedNormalDistributionGenerator(
        lConfig=config,
        lAverageNucleotideQuality=args.nucleotide_quality,
        lStdev=args.nucleotide_stddev
    )


    buffers = generateRandomBuffers(
        lTruncatedNormalRandomGenerator=truncatedNormalRandomGenerator,
        lRandomNumberGenerator=randomNumberGenerator,
        lConfig = config,
    )



    cleanReadStartPositions = filterBadReads(
        lReadStartPositions=readStartPositions,
        lNucleotideGenome=nucleotideGenome,
        lInsertSize=insertSize,
        lConfig=config
    )


    print("Generating reads")

    frontList, backList = generateReads(
        cleanReadStartPositions,
        nucleotideGenome,
        buffers = buffers,
        lArgs = args
    )

    print("Writing files")

    writeReadsToFqFile(args.output_fq_file1, frontList)
    writeReadsToFqFile(args.output_fq_file2, backList)

    writeSamFile(
        lFilePath=args.output_sam_file,
        lNucleotideGenome=nucleotideGenome,
        lHeaderLine=headerLine,
        lFrontReadList=frontList,
        lBackReadList=backList,
    )


if __name__ == '__main__':
    main()
