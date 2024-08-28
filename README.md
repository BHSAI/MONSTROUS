# MONSTROUS

MOlecular traNSporT inhibitoR and substrate predictOr Utility Server (MONSTROUS) is a computational transporter profiler that predicts the potential of a chemical to interact with transporters recommended for testing in drug development by regulatory agencies. Currently, these transporters are considered to be a major player in determining the safety and efficacy of drugs. MONSTROUS utilizes either graph convolutional neural networks or similarity-based cheminformatics approaches to screen query chemicals against 12 transporters widely expressed in various tissues, including liver, brain, and kidney, and makes predictions as to their potential to be inhibitors as well as substrates.

##### Supporting information for paper:
Title: 
Authors: 
Journal: 
### Intro
This repository contains the data and models used to make MONSTROUS's predictions and has sorted this data into 4 sections:
- A GCNN folder containing the data for our GCNN transporters. This data includes csv files containing lists of compounds for each GCNN transporter as well as their GCNN models.
- A Similarity Approach folder containing the data for our non-GCNN transporters. This data includes csv files containing lists of compounds for each transporter.
- A python folder containing the code that will run the MONSTROUS command line tool.
- An examples folder containing an example input as well as output files for each output format the tool supports.

### GCNN
The GCNN folder contains two subfolders: compounds and models. In the compounds folder are CSV files for each transport protein, containing a list of reference compounds that are known inhibitors or substrates for the given transporter. These compounds are used in generating the applicability domain for that transporter. In the models folder we hold the models for each GCNN transporter. These models are used to generate the values for GCNN transporters

### Similarity Approach
The similarity approach folder contains CSV files for each transport protein, containing a list of reference compounds that are known inhibitors or substrates for the given transporter. These compounds are used in generating the applicability domain for that transporter, as well as generating the values for similarity approach transporters

## MONSTROUS Command Line Tool

### Setting up the anaconda environment
To run the python script, you will first need Anaconda installed. From an Anaconda prompt, set up a new environment using the following commands:

`conda create -n monstrous python=3.9`

`conda activate monstrous`

Next, navigate to this repository's folder and enter the following command to install MONSTROUS's dependencies:

`pip install -r requirements.txt`

### Running the MONSTROUS command line tool

Once everything is installed, you can then run the script by running `python python/monstrous_clt.py` followed by any of the following tags (and must include the `-i` , input file and tag):
- `-h` or `--help`: Shows a help message explaining these tags.
- `-i [INPUT]` or `--input [INPUT]`: The file location of a .CSV file whose first column is 'Name' and whose second is 'SMILES' and contains the list of SMILES to be submitted.
- `-o [OUTPUT]` or `--output [OUTPUT]`: The output file path
- `-f [FORMAT]` or `--format [FORMAT]`: File format for the output. The options are 'json', 'csv', or 'database'. The JSON format produces a JSON file, the CSV format produces 4 CSV files (for inhibitor and substrate with and without the applicability domain applied), and the database format produces a csv file with 1 result per row. 
- `-e [ERROR]` or `--error [ERROR]`: The file path of a .txt file for outputting error messages.
- `-s [STATE]` or `--state [STATE]`: The file path of a .txt file for outputting state info.

Here are some example prompts:

- `python python/monstrous_clt.py -i examples/example_input.csv -o examples -f json`
	- [output](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_json_ouput.json) 
- `python python/monstrous_clt.py -i examples/example_input.csv -o examples -f database`
	- [output](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_database_output.csv)
- `python python/monstrous_clt.py -i examples/example_input.csv -o examples -f csv`
	- [inhibitor output](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_csv_output_inhibitor.csv)
	- [inhibitor output with applicability domain](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_csv_output_inhibitor_with_applicability_domain.csv)
	- [substrate output](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_csv_output_substrate.csv)
	- [substrate output with applicability domain](https://github.com/BHSAI/MONSTROUS/blob/main/examples/example_csv_output_substrate_with_applicability_domain.csv)

These prompts take the input file from the `examples` folder and output back to the example folder.
