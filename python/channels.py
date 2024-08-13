import sys

class Channels:
    errorChannel = sys.stderr
    stateChannel = None
    outputChannels = [sys.stdout]

    def writeStateMsg(self, msg):
        if self.stateChannel != None:
            print(msg, file= self.stateChannel)
            self.stateChannel.flush()

    def writeErrorMsg(self, msg):
        if self.errorChannel != None:
            print(msg, file= sys.stderr)
            print(msg, file= self.errorChannel)
            self.errorChannel.flush