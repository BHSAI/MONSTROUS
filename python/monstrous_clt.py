import argparse
import csv
from channels import Channels
from output import *
from mainModels import runMainModels
from preprocessing import processCompounds
from datetime import datetime

from compound import *

class InputBundle:
    names = None      
    smiles = None     
    outputPath = None 
    fileFormat = 'csv'
    errorPath = None  
    statePath = None 
    batch_size = None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", nargs=1, help="File location of input csv")
    parser.add_argument("-o", "--outputLocation", nargs=1, help="Output file path")
    parser.add_argument("-f", "--format", nargs=1, help="File format for the output. either csv, json or database. default: csv")
    # Batch size removed for public command line tool
    # parser.add_argument("-b", "--batchSize", nargs=1, help="The amount of compounds that go through Main Predictions each time in order to avoid any memory issue. By default, it is equal to the number of smiles submitted.")
    parser.add_argument("-e", "--error", nargs=1, help="The file path for error messages")
    parser.add_argument("-s", "--state", nargs=1, help="The file path for outputting state info.")
    args = parser.parse_args()

    inputs = InputBundle()
    if args.inputFile != None and args.inputFile[0] != None:
        inputs.names, inputs.smiles = parseInputCsv(args.inputFile[0])
    else:
        parser.error("Must include input file (-i 'input/file.csv')")

    if args.outputLocation != None and len(args.outputLocation) > 0:
        inputs.outputPath = args.outputLocation[0]
    else:
        inputs.outputPath = "./" #current work directory
    if args.error != None and len(args.error) > 0:
        inputs.errorPath = args.error[0]
    if args.state != None and len(args.state) > 0:
        inputs.statePath = args.state[0]

    args.batchSize = None # Added for public command line tool
    if args.batchSize != None and len(args.batchSize) > 0:
        if args.batchSize[0].isdigit():
            if int(args.batchSize[0]) <= 0:
                parser.error("'batchSize' must be > 0")
            inputs.batch_size = int(args.batchSize[0])
        else:
            parser.error("'batchSize' must be an integer")
    if args.format != None and len(args.format) > 0:
        inputs.fileFormat = args.format[0]
        if inputs.fileFormat not in ['csv', 'json','database']:
            parser.error(f"Invalid file format \"{inputs.fileFormat}\". Expecting csv, json, or database")

    return inputs

def parseInputCsv(csv_path):
    print(f'csv_path: {csv_path}')
    with open(csv_path, 'r') as file:
        csvreader = csv.reader(file)
        names = []
        smiles = []
        for row in csvreader:
            if len(row) == 1:
                smiles.append(row[0])
            else:
                names.append(row[0])
                smiles.append(row[1])
        if len(names) > 0:
            names.pop(0)
        if len(smiles) > 0:
            smiles.pop(0)
        return names, smiles

def main(inputs):
    channels = Channels()
    if inputs.outputPath != None:
        if inputs.fileFormat == "csv":
            channels.outputChannels[0] = open(inputs.outputPath + f"/Monstrous_predictions_inhibitor_with_applicability_domain_{datetime.now().strftime('%Y-%m-%d_%H%M')}.csv", 'w')
            channels.outputChannels.append(open(inputs.outputPath + f"/Monstrous_predictions_substrate_with_applicability_domain_{datetime.now().strftime('%Y-%m-%d_%H%M')}.csv", 'w'))
            channels.outputChannels.append(open(inputs.outputPath + f"/Monstrous_predictions_inhibitor_{datetime.now().strftime('%Y-%m-%d_%H%M')}.csv", 'w'))
            channels.outputChannels.append(open(inputs.outputPath + f"/Monstrous_predictions_substrate_{datetime.now().strftime('%Y-%m-%d_%H%M')}.csv", 'w'))
        elif inputs.fileFormat == "json":
            channels.outputChannels[0] = open(inputs.outputPath + f"/Monstrous_predictions_{datetime.now().strftime('%Y-%m-%d_%H%M')}.json", 'w')
        elif inputs.fileFormat == "database":
            channels.outputChannels[0] = open(inputs.outputPath + f"/Monstrous_predictions_{datetime.now().strftime('%Y-%m-%d_%H%M')}.csv", 'w')
    if inputs.errorPath != None:
        channels.errorChannel = open(inputs.errorPath, 'w')
    if inputs.statePath != None:
        channels.stateChannel = open(inputs.statePath, 'w')

    compounds = getCompounds(inputs.names, inputs.smiles)

    if inputs.batch_size == None:
        inputs.batch_size = len(inputs.smiles)

    processCompounds(compounds, channels)

    runMainModels(compounds, inputs.batch_size, channels)

    write_output(inputs.fileFormat, compounds, channels.outputChannels)

    channels.writeStateMsg("Execution finished")
    
    if inputs.errorPath != None:
        channels.errorChannel.close()
    if inputs.statePath != None:
        channels.stateChannel.close()
    if inputs.outputPath != None:
        for channel in channels.outputChannels:
            channel.close()


   

if __name__ == "__main__":
    inputs = parse_args()
    main(inputs)