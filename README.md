# SNP_haplotype
A python script for outputting haplotype networks in NEXUS file format for use in Popart.

Inputs: Aligned Fasta file, outbreak metadata text file

Outputs: SNP csv, pseudohaplotype.txt, rough haplotype network.png, haplotype.NEX

## Usage
```python3
pipeline.py input_path output_path metadata_path outbreak_only outbreaks n_threshold remove_n
```
## Arguments:
Description

(type)

Example

-----
input_path = The path to the fasta file

(str)

'/mnt/bio-tarako-home/zgls201/Documents/data.fasta'

output_path = The path to the output directory

(str)

'/mnt/bio-tarako-home/zgls201/Documents/results'

metadata_path = The path to the outbreak metadata txt file in the format below (str)

(str)

'/mnt/bio-tarako-home/zgls201/Documents/outbreak_metadata.txt

outbreak_only = A T/F statement choosing whether to obtain results for just the samples in the outbreak_metadata.txt.
		True: only outbreak samples
		False: All samples in data.fasta and outbreak samples
(bool)

True

outbreaks = Which outbreaks should be run e.g
	    -a    All outbreaks
	    None  No outbreaks
	    1,7	  Outbreaks 1 and 7
	    
(str)

-a / None / 1,7,8,19

n_threshold = The number of N's permitted in a sample sequence before exclusion. Samples with more than the N_threshold will be excluded from the run

(int)

5000

remove_N = A T/F statement choosing to exclude/include SNP positions where there is an N in any sample.

(bool)

True


