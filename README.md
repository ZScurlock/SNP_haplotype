# SNP_haplotype
A python script for outputting haplotype networks in NEXUS file format for use in Popart haplotyping network software.
Inputs a multi-aligned fasta file for SNP and haplotype creation between samples. Within this file, samples may belong to specific outbreaks noted in the outbreak metadata text file. Certain outbreaks for analysis can be specified to analyse only specific combinations of outbreaks.

Inputs: Aligned Fasta file, outbreak metadata text file (see template_outbreak_meta.txt)

Outputs: SNP.csv, pseudohaplotype.txt, haplotype_network.png, haplotype.nex, outbreak.csv, outbreak_pseudohaplotype.txt, outbreak_network.png, outbreak.nex

### Dependencies:
Python 3.6.9,
Pandas 0.22.0,
Numpy 1.19.4,
Networkx 2.5,
matplotlib 2.1.1,
pairsnp 0.07,


## Usage
```python3
pipeline.py input_path output_path metadata_path outbreak_only outbreaks n_threshold remove_n
```
## Arguments:
Example_format

Argument = Description

(type)

Example

-----
input_path = The path to the fasta file

(str)

'/zscurlock/Documents/data.fasta'

------

output_path = The path to the output directory

(str)

'/zscurlock/Documents/results'

------

metadata_path = The path to the outbreak metadata txt file in the format below (str)

(str)

'/zscurlock/Documents/outbreak_metadata.txt

------

outbreak_only = A T/F statement choosing whether to obtain results for just the samples in the outbreak_metadata.txt.
		True: only outbreak samples
		False: All samples in data.fasta and outbreak samples
		
(bool)

True

------

outbreaks = Which outbreaks should be run e.g
	    -a    All outbreaks
	    None  No outbreaks
	    1,7	  Outbreaks 1 and 7
	    
(str)

-a / None / 1,7,8,19

------

n_threshold = The number of N's permitted in a sample sequence before exclusion. Samples with more than the n_threshold will be excluded from the run

(int)

5000

------

remove_n = A T/F statement choosing to exclude/include SNP positions where there is an N in any sample.

(bool)

True


