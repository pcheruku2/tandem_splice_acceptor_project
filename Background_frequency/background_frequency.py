import os,sys
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Background frequency calculation')
parser.add_argument("-w","--window", type=str, help="33 bases or 31 bases")
parser.add_argument("-n","--number", type=str, help="number of sequences")
parser.add_argument("-s","--strand", type=int, help="0 for postive strand, 1 for negative strand")
args = parser.parse_args()

window = args.window
number = args.number
strand = args.strand
if strand not in [0,1]:
    print('Please select 0 or 1 for strand. 0 for postive strand, 1 for negative strand')
    sys.exit()

data = pd.DataFrame(index=['AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT','AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT','CAG','CCG','CGG','CTG','GAG','GCG','GGG','GTG','TAG','TCG','TGG','TTG'])

for i in range(10):
    cmd = "perl random_bed.pl confirmed_chr_lengths.txt "+window+" "+number
    os.system(cmd)
    cmd = "samtools faidx hg19_ref.fa.gz -r "+"random_w"+window+"_n"+number+".sam.txt -o random.fa"
    os.system(cmd)
    cmd = "perl filter.pl random.fa > random_filtered.fa"
    os.system(cmd)
    if strand ==0:
        cmd = "perl trimer_count.pl uniq_trimers_AG_gain_loss_neutral.txt random_filtered.fa random"
        os.system(cmd)
    else:
        cmd = "perl reverse_complement.pl random_filtered.fa > random_filtered_revComp.fa"
        os.system(cmd)
        cmd = "perl trimer_count.pl uniq_trimers_AG_gain_loss_neutral.txt random_filtered_revComp.fa random"
        os.system(cmd)  
    data_tmp = pd.read_csv("counts_random.out",sep="\t",index_col=0,names=["Count"+str(i)])
    if strand == 0:
        num_lines = sum(1 for _ in open('random_filtered.fa'))
        total_trimers = (int(window)-2)*(num_lines/2) 
    else:
        num_lines = sum(1 for _ in open('random_filtered_revComp.fa'))
        total_trimers = (int(window)-2)* (num_lines/2)       
    data_tmp["Freq"+str(i)] = data_tmp["Count"+str(i)] / total_trimers
    data = data.merge(data_tmp,how="left",left_index=True, right_index=True)
data.to_csv("Frequency_output.csv")