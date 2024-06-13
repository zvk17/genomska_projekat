import asyncio
import os
import subprocess
from pathlib import Path, PurePath
from tempfile import NamedTemporaryFile
from sklearn import metrics
import numpy
import seaborn
from matplotlib import pyplot as plt
from pysam import AlignmentFile, AlignedRead, AlignedSegment

simulatorFile = "./output_files/alignment.sam"
testFile = "./output_files/out_bwa.sam"

BWA_PATH = ""
BOWTIE_PATH = ""
BOWTIE_BUILD_PATH = ""
SIMULATOR_PATH = "./simulator.py"

input_genomes = [
    "./input_files/chr1c.fna",
    "./input_files/chrX.fna",
    "./input_files/daucus_carota_chr4.fna"
]




def readSimulatorSamFile(lSimulatorSam: AlignmentFile):
    frontReadDict = {}
    backReadDict = {}
    for read in lSimulatorSam.fetch():

        if (not read.is_read2):
            frontReadDict[read.qname] = read
        else:
            backReadDict[read.qname] = read
    return frontReadDict, backReadDict

def getReadWithSameName(
        lFrontReadDict: dict,
        lBackReadDict:dict,
        lRead: AlignedSegment
)->AlignedSegment:
    if lRead.is_read1 and (lFrontReadDict.get(lRead.qname) is not None):
       return lFrontReadDict[lRead.qname]

    if lRead.is_read2 and (lBackReadDict.get(lRead.qname) is not None):
        return lBackReadDict[lRead.qname]
    return None

def compareReads(lRead1: AlignedSegment, lRead2: AlignedRead,lEpsilon:float) -> bool:
    if abs(lRead1.pos - lRead2.pos) < lEpsilon:
        return True
    return False

def bwaMemQualityToProbability(lQuality: int)-> float:
    return pow(10, -lQuality * 1.0 / 10)

def getPredictionsDeprecated(
        lFrontReadDict: dict,
        lBackReadDict: dict,
        lTestSamFile: AlignmentFile,
        lEpsilon: int,
):
    lUnmappedReads = []
    lY = []
    lPred = []
    lQuals = []

    for lRead in lTestSamFile.fetch():
        if lRead.is_unmapped:
            lUnmappedReads.append(lRead)
            continue

        lSimulatorRead = getReadWithSameName(
            lFrontReadDict=lFrontReadDict,
            lBackReadDict=lBackReadDict,
            lRead=lRead
        )
        if lSimulatorRead is None:
            continue
        lResult: bool = compareReads(lRead1=lSimulatorRead, lRead2=lRead, lEpsilon=lEpsilon)
        lY.append(lResult)
        lPred.append(lRead.mapping_quality / 60.0)

        lQuals.append(lRead.mapping_quality)

    return lY, lPred, lQuals, lUnmappedReads


def printConfusionMatrixFromPredictionsDepracted(lY, lPred, lOutputFilePath:str):
    totalTrue = sum(lY)
    totalFalse = len(lY) - totalTrue

    lConfusionMatrix = metrics.confusion_matrix(lY, lPred)
    print(lConfusionMatrix)
    plt.clf()
    metrics.ConfusionMatrixDisplay(lConfusionMatrix).plot()
    plt.savefig(lOutputFilePath)

def plotAucRoc(lY: list, lPred:list, lOutputFilePath: str):
    lFpr, lTpr, _ = metrics.roc_curve(lY, lPred)
    lAuc = metrics.roc_auc_score(lY, lPred)
    plt.plot(lFpr, lTpr, label="data 1, auc=" + str(lAuc))
    plt.plot(lFpr, lTpr, label="data 1, auc=" + str(lAuc))
    plt.legend(loc=4)
    plt.savefig(lOutputFilePath)

def calcConfusionMatrix(
        frontReadDict: dict,
        backReadDict: dict,
        testSamFile: AlignmentFile,
        lEpsilon: int | float
):
    lTruePositive = 0
    lFalsePositive = 0
    lTrueNegative = 0
    lFalseNegative = 0


    for lRead in testSamFile.fetch():

        if lRead.is_unmapped:
            continue

        lSimulatorRead = getReadWithSameName(
            lFrontReadDict=frontReadDict,
            lBackReadDict=backReadDict,
            lRead=lRead
        )
        if lSimulatorRead is None:
            lFalsePositive = lFalsePositive + 1
            continue

        if compareReads(lSimulatorRead, lRead, lEpsilon):
            lTruePositive = lTruePositive + 1
        else:
            lFalsePositive = lFalsePositive + 1

    return [
        lTruePositive, lFalsePositive, lTrueNegative, lFalseNegative
    ]
def compareSamFilesDeprecated(
    lSimulatorSamFilePath: str,
    lAlignedSamFile: str
):
    lSimulatorSam = AlignmentFile(lSimulatorSamFilePath, "r", check_sq=False, ignore_truncation=True)
    frontReadDict, backReadDict = readSimulatorSamFile(lSimulatorSam)
    lAlignedSam = AlignmentFile(lAlignedSamFile, "r", check_sq=False, ignore_truncation=True)
    lEpsilon = 50
    lConfusionMatrix = (
        getPredictionsDeprecated(frontReadDict, backReadDict, lAlignedSam, lEpsilon))

    return lConfusionMatrix

def compareSamFiles(
    lSimulatorSamFilePath: str,
    lAlignedSamFile: str
):
    lSimulatorSam = AlignmentFile(lSimulatorSamFilePath, "r", check_sq=False, ignore_truncation=True)
    frontReadDict, backReadDict = readSimulatorSamFile(lSimulatorSam)
    lAlignedSam = AlignmentFile(lAlignedSamFile, "r", check_sq=False, ignore_truncation=True)
    lEpsilon = 50
    lConfusionMatrix = (
        calcConfusionMatrix(frontReadDict, backReadDict, lAlignedSam, lEpsilon))

    return lConfusionMatrix
def plotConfusionMatrix(
        truePositive,
        falsePositive,
        trueNegative,
        falseNegative,
        lOutputFilePath: str,
):
    arr = numpy.array([[truePositive, falsePositive], [trueNegative, falseNegative]])
    plt.clf()
    seaborn.heatmap(
        data=arr,
        annot=True,
        xticklabels=["true", "false"],
        yticklabels=["positive", "negative"],
        cmap = "viridis",#crest
        linewidth=2.5,
        fmt=',d'
    )
    plt.savefig(lOutputFilePath)


def printResults(lConfusionMatrix, lImageOutputFilePath: str):
    [truePositive, falsePositive, trueNegative, falseNegative] = lConfusionMatrix
    print("truePositive = " + str(truePositive))
    print("falsePositive = " + str(falsePositive))
    print("trueNegative = " + str(trueNegative))
    print("falseNegative = " + str(falseNegative))
    print("accuracy = " + str((truePositive + trueNegative) * 1.0 / (truePositive +trueNegative + falseNegative+ falsePositive)))
    print("precision = " + str((truePositive) * 1.0 / (truePositive + falsePositive)))
    print("recall = " + str((truePositive) * 1.0 / (truePositive + falseNegative)))

    plotConfusionMatrix(
        truePositive,
        falsePositive,
        trueNegative,
        falseNegative,
        lImageOutputFilePath
    )




def getSimulatorCommand(
        lSimulatorPath: str,
        lInputFile: str,
        lOutputSamFile:str,
        lFrontFqFile: str,
        lBackFqFile: str,
        lNucleotideQuality: float,
        lCoverage: float,
        lPythonPath: str = "python",
        lSnvErrorRate: float = 0.0,
        lInsertionErrorRate: float = 0.0,
        lDeletionErrorRate: float = 0.0,
        lNucleotideStddev: float = 1.0,
):
    if lCoverage < 0.0:
        raise Exception("lCoverage must be positive")
    if lSnvErrorRate < 0.0 or lInsertionErrorRate < 0.0 or lDeletionErrorRate < 0.0:
        raise Exception("error rate must be positive")

    return [
        lPythonPath,
        lSimulatorPath,
        "--input-file", lInputFile,
        "--output-sam-file", lOutputSamFile,
        "--output-fq-file1", lFrontFqFile,
        "--output-fq-file2", lBackFqFile,
        "--nucleotide-quality", str(lNucleotideQuality),
        "--nucleotide-stddev", str(lNucleotideStddev),
        "--coverage", str(lCoverage),
        "--snv-error-rate", str(lSnvErrorRate),
        "--insertion-error-rate", str(lInsertionErrorRate),
        "--deletion-error-rate", str(lDeletionErrorRate),

    ]

def executeSimulator(
        lInputFile: str,
        lNucleotideQuality: float,
        lCoverage: float,
        lSnvErrorRate: float = 0.0,
        lInsertionErrorRate: float = 0.0,
        lDeletionErrorRate: float = 0.0,
        lNucleotideStddev: float = 1.0
):
    lFq1File = NamedTemporaryFile(mode="w", delete=False)
    lFq2File = NamedTemporaryFile(mode="w", delete=False)
    lOutputSamFile = NamedTemporaryFile(mode="w", delete=False)
    lCommand = getSimulatorCommand(
        lSimulatorPath=SIMULATOR_PATH,
        lInputFile=lInputFile,
        lFrontFqFile=lFq1File.name,
        lBackFqFile=lFq2File.name,
        lOutputSamFile=lOutputSamFile.name,
        lNucleotideQuality=lNucleotideQuality,
        lCoverage=lCoverage,
        lSnvErrorRate=lSnvErrorRate,
        lInsertionErrorRate=lInsertionErrorRate,
        lDeletionErrorRate=lDeletionErrorRate,
        lNucleotideStddev=lNucleotideStddev

    )
    print(lCommand)
    print("Simulator started for file: <{}>".format(lInputFile))
    lResult = subprocess.run(lCommand, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("Simulator finished for file: <{}>".format(lInputFile))
    return lResult,lFq1File,lFq2File,lOutputSamFile

def doesBwaIndexExist(lInputFilePath: str)-> bool:
    lPath = Path(lInputFilePath + ".bwt")
    return lPath.exists() and lPath.is_file()


def doesBowtieIndexExist(lInputFilePath: str)-> bool:
    lPath = Path(lInputFilePath).absolute()
    lIndexFilePath = lPath.parent.joinpath(lPath.stem + "_bowtie_index.1.bt2")
    return lIndexFilePath.exists() and lIndexFilePath.is_file()
def getBwaIndexCommand(lBwaPath: str,lInputFile: str,):
    return [lBwaPath, "index", lInputFile]

def getBowtieIndexCommand(lBowtieBuildPath: str, lInputFile: str):
  
    lOutputFileIndexPath = getBowtieIndexFromFile(lInputFile)
    return [
        lBowtieBuildPath,
        lInputFile,
        str(lOutputFileIndexPath)
    ]

def executeBwaIndex(lInputFile: str) -> subprocess.CompletedProcess:
    print("BWA indexing {} started".format(lInputFile))
    lCommand = getBwaIndexCommand(lBwaPath=BWA_PATH, lInputFile=lInputFile)
    lResult = subprocess.run(lCommand, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if lResult.returncode == 0:
        print("BWA indexing {} finished".format(lInputFile))
    return lResult
def executeBowtieIndex(lInputFile: str) -> subprocess.CompletedProcess:
    print("Bowtie indexing {} started".format(lInputFile))
    lCommand = getBowtieIndexCommand(lBowtieBuildPath=BOWTIE_BUILD_PATH, lInputFile=lInputFile)
    print(lCommand)
    lResult = subprocess.run(lCommand, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if lResult.returncode == 0:
        print("BOWTIE indexing {} finished".format(lInputFile))
    return lResult
def getBwaMemCommand(
        lBwaPath: str,
        lInputFile: str,
        lOutputFile:str,
        lFrontFqFile:str,
        lBackFqFile:str
) -> list:
    return [
        lBwaPath,
        "mem",
        lInputFile,
        "-o",
        lOutputFile,
        lFrontFqFile,
        lBackFqFile
    ]
def getBowtieCommand(
        lBowtiePath: str,
        lInputIndexFileRootName: str,
        lOutputFile:str,
        lFrontFqFile:str,
        lBackFqFile:str
) -> list:

    return [
        lBowtiePath,
        "-x",
        lInputIndexFileRootName,
        "-1",lFrontFqFile,
        "-2",lBackFqFile,
        "-S",lOutputFile,


    ]
def getBowtieIndexFromFile(lFilePath: str):
    lPath = Path(lFilePath).absolute()
    return str(lPath.parent.joinpath(lPath.stem + "_bowtie_index"))
def executeBowtie(
        lFrontFqFile: str,
        lBackFqFile: str,
        lInputFile: str
):
    lOutputSamFile = NamedTemporaryFile(mode="w", delete=False)
    lCommand = getBowtieCommand(
        lBowtiePath=BOWTIE_PATH,
        lInputIndexFileRootName=getBowtieIndexFromFile(lInputFile),
        lOutputFile=lOutputSamFile.name,
        lFrontFqFile=lFrontFqFile,
        lBackFqFile=lBackFqFile,
    )
    print("Bowtie  {} started".format(lInputFile))
    lResult = subprocess.run(lCommand, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("Bowtie finised: ")
    print("file size ", lOutputSamFile.name)
    lOutputSamFile.close()


    return lResult, lOutputSamFile

def executeBwaMem(
    lFrontFqFile: str,
    lBackFqFile: str,
    lInputFile: str
):
    lOutputSamFile = NamedTemporaryFile(mode="w", delete=False)
    lCommand = getBwaMemCommand(
        lBwaPath=BWA_PATH,
        lInputFile=lInputFile,
        lOutputFile=lOutputSamFile.name,
        lFrontFqFile=lFrontFqFile,
        lBackFqFile=lBackFqFile,
    )
    lResult = subprocess.run(lCommand, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("bwa mem result: ", lResult)
    print("file size ", lOutputSamFile.name)

    lOutputSamFile.close()


    return lResult, lOutputSamFile



def runTest(lNucleotidePath: str, lDeleteFilesOnExit: bool = False):
    if not doesBwaIndexExist(lInputFilePath=lNucleotidePath):
        lResult = executeBwaIndex(lInputFile=lNucleotidePath)
        if lResult.returncode != 0:
            raise Exception("BWA index subprocess error")
    else:
        print("Skipping building BWA-MEM index")
    if not doesBowtieIndexExist(lInputFilePath=lNucleotidePath):
        lResult = executeBowtieIndex(lInputFile=lNucleotidePath)
        if lResult.returncode != 0:
            raise Exception("BOWTIE index subprocess error")
    else:
        print("Skipping building Bowtie index")

    lTestErrorRates = [0.00, 0.025,  0.050,  0.100, 0.150, 0.200]
    lResults = []
    for lErrorRate in lTestErrorRates:

        lSimulatorResult = executeSimulator(
            lInputFile=lNucleotidePath,
            lCoverage=0.1,
            lNucleotideQuality=20.0,
            lNucleotideStddev=20.0,
            lSnvErrorRate=lErrorRate,
        )
        if lSimulatorResult[0].returncode != 0:
            raise "Simulator crashed"


        lBwaMemResult = executeBwaMem(
            lFrontFqFile=lSimulatorResult[1].name,
            lBackFqFile=lSimulatorResult[2].name,
            lInputFile=lNucleotidePath
        )
        if lBwaMemResult[0].returncode != 0:
            raise Exception("BWA MEM crashed")

        lConfusionMatrix = compareSamFiles(
            lSimulatorSamFilePath=lSimulatorResult[3].name,
            lAlignedSamFile=lBwaMemResult[1].name,
        )
        printResults(
            lConfusionMatrix=lConfusionMatrix,
            lImageOutputFilePath="./output_files/conf_mat_bwamem_{}.png".format(int(lErrorRate * 1000))
        )
        lBowtieResult = executeBowtie(
            lFrontFqFile=lSimulatorResult[1].name,
            lBackFqFile=lSimulatorResult[2].name,
            lInputFile=lNucleotidePath
        )
        if lBowtieResult[0].returncode != 0:
            raise Exception("Bowtie crashed")

        lConfusionMatrix = compareSamFiles(
            lSimulatorSamFilePath=lSimulatorResult[3].name,
            lAlignedSamFile=lBowtieResult[1].name,
        )
        printResults(
            lConfusionMatrix=lConfusionMatrix,
            lImageOutputFilePath="./output_files/conf_mat_bowtie_{}.png".format(int(lErrorRate * 1000))
        )


def main():

    runTest(
        lNucleotidePath=input_genomes[1],
    )
    return






if __name__ == '__main__':
    main()


