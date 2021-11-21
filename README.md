# FastQC-analogue
Implementation of an analogue of the FastQC program <br>

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. You can find more information here: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ <br>

The goal is to implement the programm which is analoge of the FastQC <br>

My team name is Lonely Snail. The only participant is Kravchuk Ekaterina (telegram: @kate_the_snail). <br>

![snail](./pictures/snail.jpg "snail")


The script was tested on Ubuntu 21.04, Python 3.9.7

## The instructions for the programm

The file requirements.txt contains a list of libraries, necessary to launch the file snail_fastgc.py. <br>
snail_fastqc works only with .fastqc files.

The programm does not support colorspace data. Please, convert it to base calls before use. <br>
<br>
1) Download snail_fastqc.py and requirements.txt <br>

2) If you don't have the library for creating virtual environments, install it 
```bash
python -m pip install virtualenv
```

3) Create the virtual environment 
```bash
python -m virtualenv any_name
```

4) Activate the virtual environment
```bash
source any_name_that_you_like/bin/activate
```

5) Install packages listed in the file requirements.txt 
```bash
python -m pip install -r requirements.txt
```

6) Finally launch snail_fastqc.py

```bash
python snail_fastqc.py --i (--input) input_file.fastq -o (--outdir) dir_to_save_output/
```
### Here is the help message for the programm:
![help](./pictures/help.jpg "help")

## The output

The programm outputs the following files:
1) input_file_name_Basic_Statistics.tsv
2) input_file_name_per_seq_qscore.png
3) input_file_name_seq_len_distribution.png
4) input_file_name_duplicate_seq.png
5) input_file_name_Overrepresented_sequences.tsv

### Basic Statistics
The data contains in input_file_name_Basic_Statistics.tsv. <br>
In the file you can find information about the input file name, file type, encoding, sequence length ang gc content.

### Per Sequence Quality Scores
input_file_name_per_seq_qscore.png shows per sequence Q-score destribution <br>
The status is 'Green (entirely normal)' if the per sequence quality score is satisfactory, 'Orange (slightly abnormal)' if the most frequently observed mean quality is below 27 - this equates
to a 0.2% error rate, 'Red (very unusual)' if the most frequently observed mean quality is below 20 - this equates to a 1% error rate.

### Sequence Length Distribution
input_file_name_seq_len_distribution.png shows sequence length destribution <br>
The status is 'Green (entirely normal)' if all the reads have the same length, 'Orange (slightly abnormal)' if all sequences are not the same length, 'Red (very unusual)' any of the sequences have zero length.

### Duplicate Sequences
input_file_name_duplicate_seq.png contains sequence duplication level distribution <br>
The status is 'Green (entirely normal)' if the feature is satisfactory, 'Orange (slightly abnormal)' if non-unique sequences make up more than 20% of the total, 'Red (very unusual)' if non-unique sequences make up more than 50% of the total.

### Overrepresented Sequences
This module lists all of the sequence which make up more than 0.1% of the total. To conserve memory only sequences which appear in the first 20000 sequences are tracked to the end of the file.It is therefore possible that a sequence which is overrepresented but doesn't appear at the start of the file for some reason could be missed by this module.<br>
The status is 'Green (entirely normal)' if there is no overrepresented sequences, 'Orange (slightly abnormal)' if any sequence is found to represent more than 0.1% of the total, 'Red (very unusual)' if any sequence is found to represent more than 1% of the total.
