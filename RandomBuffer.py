import numpy
class RandomBuffer:
    def __init__(self, randomNumberGenerator, bufferSize):
        self.randomNumberGenerator = randomNumberGenerator
        self.buffer = self.randomNumberGenerator.uniform(size = bufferSize)
        self.bufferSize = bufferSize
        self.position = 0
    def take(self, size):
        assert(size < self.bufferSize)
        if (self.bufferSize - self.position < size):
            self.buffer = self.randomNumberGenerator.uniform(size=self.bufferSize)
            self.position = 0

        part = self.buffer[self.position:self.position + size]
        self.position = self.position + size
        return part

class RandomChoiceBuffer:
    def __init__(self, randomNumberGenerator, bufferSize, choiceArray):
        self.randomNumberGenerator = randomNumberGenerator

        self.choiceArray = choiceArray
        self.buffer = self.randomNumberGenerator.choice(
            choiceArray,
            size=bufferSize
        )
        self.bufferSize = bufferSize
        self.position = 0
    def take(self, size):
        assert(size < self.bufferSize)
        if (self.bufferSize - self.position < size):
            self.buffer = self.randomNumberGenerator.choice(
                self.choiceArray,
                size=self.bufferSize
            )
            self.position = 0

        part = self.buffer[self.position:self.position + size]
        self.position = self.position + size
        return part


class RandomNucleotideQualityBuffer:
    def __init__(self, randomNumberGenerator, bufferSize):
        self.randomNumberGenerator = randomNumberGenerator


        self.buffer = self.randomNumberGenerator.rvs(
            size=bufferSize
        )
        self.bufferSize = bufferSize
        self.position = 0
    def take(self, size):
        assert(size < self.bufferSize)
        if (self.bufferSize - self.position < size):
            self.buffer = self.randomNumberGenerator.rvs(
                size = self.bufferSize
            )
            self.position = 0

        part = self.buffer[self.position:self.position + size]
        self.position = self.position + size
        return part
