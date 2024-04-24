# Background Frequency
Randomly selecting 1,000, 10,000, 100,000 and 1,000,000 windows of 33 base pairs from the human genome reference sequence (GRCh37/ hg19) and calculated the background frequency 

## Requirements
1. Install [Samtools](http://www.htslib.org/download/), version 1.20 or up
2. Download [hg19_ref.fa.gz](link_will_be_provided)
3. Python3
4. Perl

```
usage: background_frequency.py [-h] [-w WINDOW] [-n NUMBER] [-s STRAND]

Background frequency calculation

options:
  -h, --help            show this help message and exit
  -w WINDOW, --window WINDOW
                        33 bases or 31 bases
  -n NUMBER, --number NUMBER
                        number of sequences
  -s STRAND, --strand STRAND
                        0 for postive strand, 1 for negative strand
```

## Example Command
```
python background_frequency.py -w 33 -n 100000 -s 0
```

## Output File: (Frequency_output.csv)

