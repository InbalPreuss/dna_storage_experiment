# DNA Storage experiment
In this code we have two different options we can choose to run:
1. Analyze the sequencing data after the sequencing of the experiment. The sequencing was done with illumina sequencing
2. Produce the molecules with the right properties for the experiment

## Step 1
Activating a virtual environment (venv)
A virtual environment is a separate Python environment, where you can install packages without affecting the global Python installation on your system.

To activate a venv on Windows:

```
cd myproject
path\to\python -m venv venv
venv\Scripts\activate
```

To activate a venv on Unix or macOS:
```
cd myproject
source venv/bin/activate
```

## Step 2
Installing requirements from requirements.txt
To install the required packages for your project, you can use the following command

```
pip install -r requirements.txt
```
This command reads the requirements.txt file and installs the listed packages and their dependencies in the active virtual environment.


## Step 3

### Run 1. Analyze the sequencing data after the sequencing of the experiment

#### Files to have inorder to run the full main.py pipline:

Make sure you have the following files:  
Sequence design file:
* sequence_design_file.dna  

Experiment results files:
* DS_S1_L001_R1_001.fastq
* DS_S1_L001_R2_001.fastq
1. Used Pear, pair-end read merger to merge R1 and R2. https://cme.h-its.org/exelixis/web/software/pear/
2. Put the output file after using Pear in this directory "analyze_sequencing_data/data/" and under this name "output_prefix.assembled.fastq"
3. Put the sequence_design_file.dna in this directory "analyze_sequencing_data/data/"
3. Make sure the config files have the correct configuration:
   4. config.py
      ```commandline
      'payload_pos': [60, 100, 140, 180],
              'amount_of_bc': 167,
              'design_len': 220,
              'payload_len': 20,
              'universal_len': 20,
              'barcode_len': 20,
              'subset_size': 5,
              'amount_of_payloads': 16,
      ```
           
   5. Make sure design.py, payload_design.py, barcodes_design.py contains the correct sequences and their positions in the sequence

In the output folder you can find all the files with the analysis output results.
The main result would be in the compare_design_to_experiment_results.csv which provides the analysis if the design matches the experiment sequencing output.
#### Run
```
cd analyze_sequencing_data
python main.py
```
Output folder will appear with all the necessary analysis:
* compare_design_to_experiment_results.csv - contains the analysis results, how many bc we recoverd 

### Run 2. Produce the molecules with the right properties for the experiment
This code produce the sequences for the combinatorial algo for a gibson assembly or for sticky ends experiment.

Make sure the config.py has the correct properties:
```commandline
config_1 = {
    'number_of_cycles': 1,
    'data_in_kb': 1,
    'bits_per_synthesis_cycle': 12,
    'subset_size': 5,
    'num_of_building_blocks': 16,
    'sequence_len': 20,
    'num_hamming_distance': 10,
    'forbidden_x_letters_in_a_row': 3,
    'each_letter_has_to_be_x_percent_occurrences_from': 0.18,
    'each_letter_has_to_be_x_percent_occurrences_to': 0.32,
    'stick_ends_overhang': 8
}
```
Run:
```commandline
cd produce_seqs_for_experiment
python main.py
```
In the sequences folder you can find the sequences that fit your criteria.
