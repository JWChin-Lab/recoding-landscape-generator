#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 14:28:32 2023
Script based on the previous landscape script
#@TS 2017
#@AC 9/4/2018
#@JF 29/8/2018
#@DC updated in July 2020
#@YG updated in Jan 2023
#LF: Run alignment

@author: ygu
"""
#==============================================================================
# New analysis by codon rather than nucleotide. Align to wt-template containing 
# potential synthetic inserts. Parse through each individual good reads.
# This script might be available for a recoded genome, so be careful when using
# it 
#==============================================================================
import pysam  
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
import glob
import time
import itertools
from operator import itemgetter
import matplotlib.pyplot as plt

start_time = time.time()
#==============================================================================
# Input
#==============================================================================

# =============================================================================
# input_bam_file = sys.argv[1]
# ref_fasta = sys.argv[2] 
# inputgenbankfile = sys.argv[3]  # only for extracting recoded position and fasta length
# =============================================================================
#Prepare for folders, file and all the others
#os.chdir("/Users/ygu/Desktop/code_test")
current_directory = os.getcwd()
current_directory, folders, files = next(os.walk(current_directory))

Input_GenBank_File = glob.glob('*.gb')[0]
switch = 1 #set to 1 if the genbank file is recoded and 0 if wt.


Filename = os.path.splitext(Input_GenBank_File)[0]
#Open the genbank file and generate a fasta output

with open(Input_GenBank_File, "r") as input_handle:
    with open(Filename+".fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")        
        SeqIO.write(sequences, output_handle, "fasta")
        
fasta_length = len(next(SeqIO.parse(Input_GenBank_File, 'genbank')).seq)

###Input_Bam_File     = '/Users/wesr/Desktop/Landscapes/Ala_schemes/vA5/A5_6/vA5_20kb_REXER_MDS42_only.fa_sorted.bam'
###ref_fasta          = '/Users/wesr/Desktop/Landscapes/Ala_schemes/vA5/A5_6/vA5_20kb_REXER_MDS42_only.fa' 

# This is for Ala scheme 7 serine and alanine recoding schemes
recoding_scheme = [{'target': 'TCT', 'recode': 'AGC'},
                   {'target': 'TCC', 'recode': 'AGC'},
                   {'target': 'TCG', 'recode': 'AGC'},
                   {'target': 'TCA', 'recode': 'AGT'},
                   {'target': 'GCG', 'recode': 'GCT'},
                   {'target': 'GCA', 'recode': 'GCT'},
                   {'target': 'TAG', 'recode': 'TAA'}]
# NGS filter
Min_Read_Length = 10 # filter out read with read length lower than indicated
MAPQ = 30 # filter out reads with mapping quality less than indicated. 0 if problems with synthetic insertion mismatching, higher to remove poor reads.

#Smoothing and artifact removal
Name_SI = 's.i.' #root name of the synthetic insertions in the genbank file, e.g. for s.i.1, s.i.2 root is s.i.
Artifact_removal_count = 2 #how many codons before the SI to not plot - adjust till artifacts vanish
Min_cov_threshold = 10 # if coverage below this do not calculate frequency. hmmm
#==============================================================================
#This part of the code is modified so that only recoding events are captured
def get_recoded_codon(genbank):
    recoded_codon_pos_list=[]
    genbank.features=sorted(genbank.features,key=lambda sf: sf.location.start.position)
    #Record only codons with recoded events
    recoded_codons=[feature for feature in genbank.features if feature.type=="misc_feature" and "to" in str(feature.qualifiers) and len(feature.location)==3]
    recoded_codons=sorted(recoded_codons,key=lambda sf: sf.location.start.position)
    #position_counter=0
    for feature in genbank.features:
        #record cds
        if feature.type=="CDS":
            current_cds_start=feature.location.start.position
            current_cds_end=feature.location.end.position
            if feature.strand==1:
                for recoded_codon in recoded_codons:
                    if recoded_codon.location.end.position>current_cds_end:
                        break
                    else:
                        if recoded_codon.location.start.position>=current_cds_start:
                            codon_pos = {'codon': recoded_codon.qualifiers["label"][0].split()[2],
                                          'start': recoded_codon.location.start.position,
                                          'end': recoded_codon.location.end.position,
                                          'position': sorted(list(range(recoded_codon.location.start.position,recoded_codon.location.end.position))), # set(range(start, end))
                                          'strand': 'f',
                                          'recode':recoded_codon.qualifiers["label"][0].split()[0]}
                            recoded_codon_pos_list.append(codon_pos)
            elif feature.strand==-1:
                for recoded_codon in recoded_codons:
                    if recoded_codon.location.end.position>current_cds_end:
                        break
                    else:
                        if recoded_codon.location.start.position>=current_cds_start:
                            codon_pos = {'codon': recoded_codon.qualifiers["label"][0].split()[2],
                                          'start': recoded_codon.location.end.position,
                                          'end': recoded_codon.location.start.position,
                                          'position': sorted(list(range(recoded_codon.location.start.position,recoded_codon.location.end.position))), # set(range(start, end))
                                          'strand': 'r',
                                          'recode':recoded_codon.qualifiers["label"][0].split()[0]}
                            recoded_codon_pos_list.append(codon_pos)
    #return sorted(list(itertools.chain.from_iterable(recoded_codon_pos_list)), key = itemgetter('start'))
    return(recoded_codon_pos_list)
 
    
def read_bam_by_row(CodonPosition, TargetCodon, RecodedCodon, FastaFilename, InputBamFile, MinReadLength, mapq, strand, swap):
    ''' 
    return output that is a dictionary containing:
                       ({'Position'       : position investigated,
                        'ID list'         : ID of all corespponding ngs read in 'Seq list',
                        'Seq list'        : all the sequence containing the targeted codon,
                        'freq'            : frequency of read at this particular position,
                        'Depth'           : ngs read depth at this position,
                        'NoneRecodedCodon': NoneRecodedCodon})
    
    == Arg ==
    TargetList         : a list of targeted positions (remember to -1 for Python 0-indexing)
        Eg. targetlist = [29,30,31]
    NoneRecodedCodon   : string of the wt codon
    RecodedCodon       : string of the recoded codon
    MinReadLength      : read less than this number will be discarded from further analysis
    MAPQ               : MAPQ quality below this value will be discarded from further analysis

    '''
    #TargetCodonseqrc = Seq(TargetCodon).reverse_complement()
    #RecodedCodonseqrc = Seq(RecodedCodon).reverse_complement()
    codon_record = []
    with pysam.AlignmentFile(InputBamFile,"rb") as bamfile: 
        #The CodonPosition is fetched from the recoded codon list file generated from get_codon function
        for read in bamfile.fetch(FastaFilename, CodonPosition[0], CodonPosition[-1]):  # take only the seq that contain targeted region, make sure targeted region is arranged by order
            try:  # file 1 sequence number 9 raised ValueError
                sequence = read.query_alignment_sequence # string the seqeunce of the read
                position = read.get_reference_positions() # list this is the positions of the read
                quality = read.mapping_quality
                length = read.reference_length #This in in principle should match the sequence length and the position length, assuming each read is equal
            except ValueError:
                pass   
                     
            if len(position) == len(sequence) and length > MinReadLength: # no indel & correct position and min length = 200
                if quality > mapq and len(position) > 0: 
                    try:     
                        target_seq = [sequence[position.index(i)].upper() for i in CodonPosition] # exact index for targeted position for each read 
                        seq = Seq(''.join(target_seq)) # Eg ['T','C','G'] are joined as 'TCG' 
                        #This part is unncessary now, because the recoded 
                        if strand =='f':
                            if seq in [TargetCodon,RecodedCodon]: # check that the seq is either nrc or rc, else is ngs errors and discarded
                                codon_record.append(str(seq))
                        else:  #reverse complements on the reverse strand
                            if seq.reverse_complement() in [TargetCodon,RecodedCodon]: # check that the seq rc is either the codon or target codon
                                codon_record.append(str(seq.reverse_complement()))
                                #print("A")
                    except ValueError:  # ValueError raised when there is a deletion in the sequence and targetlist position is not present
                        pass
            

        if swap < 1: #This is for WT check     
            try: 
                TotalRecoded = codon_record.count(RecodedCodon)
                freq = round((float(TotalRecoded)/len(codon_record)),2)
                ambig = round((float(TotalRecoded)/len(codon_record)),2)#Not sure what this is for
            except ZeroDivisionError:
                freq = 0
                ambig = 0
            output = {'Position'        : CodonPosition[0],
                      'Frequency'       : freq,
                      'Ambiguity'       : ambig,
                      'Depth'           : len(codon_record),
                      'Target Codon'    : TargetCodon,
                      'Recoded Codon'   : RecodedCodon,
                      'Entries'         : codon_record,
                      'Strand'          : strand,}
        else:      
            try: 
                TotalRecoded = codon_record.count(TargetCodon)
                freq = round((float(TotalRecoded)/len(codon_record)),2)
                ambig = round((float(TotalRecoded)/len(codon_record)),2)
            except ZeroDivisionError:
                freq = 0
                ambig = 0
            output = {'Position'        : CodonPosition[0],
                      'Frequency'       : freq,
                      'Ambiguity'       : ambig,
                      'Depth'           : len(codon_record),
                      'Target Codon'    : RecodedCodon,
                      'Recoded Codon'   : TargetCodon,
                      'Entries'         : codon_record,
                      'Strand'          : strand,}
    return output

#==============================================================================
# Calculate codon frequency, main function part 
#==============================================================================

mds = SeqIO.read(Input_GenBank_File, "genbank")

recoded_codon_list = get_recoded_codon(mds)

summary_data_frame = {}
cumulative_frequency = []
total_clones = 0


for subfolder in folders:
    InputBamFile = glob.glob(subfolder + '/*.bam')
    if InputBamFile:
        codon_frequency = [read_bam_by_row(CodonPosition = i['position'], 
                                           TargetCodon = i['codon'], 
                                           RecodedCodon = i['recode'], 
                                           FastaFilename = Filename,
                                           InputBamFile = InputBamFile[0], 
                                           MinReadLength = Min_Read_Length, 
                                           mapq = MAPQ,
                                           strand = i['strand'], 
                                           swap = switch) for i in recoded_codon_list]
            
          
        df = pd.DataFrame.from_dict(codon_frequency)
        if 'Position' in df: # and 'Depth' in df and 'Frequency' in df:
            position  = list(df['Position']+1) #adjusts for python indexing to corresponds to the genbank file
            depth     = list(df['Depth'])
            frequency = list(df['Frequency'])
            ambiguity = list(df['Ambiguity'])
        
        #smooth frequencies for low coverage  
            for i in range(len(position)):
                if i > 1:
                    if depth[i] < Min_cov_threshold:
                        frequency[i] = frequency[i-1]  
                        ambiguity[i] = 0.5 
            
        # export data to a csv file
            df.to_csv(subfolder + '/' + Filename + '_' + subfolder + '.csv')
        
            # =============================================================================
            # # backup for manual entry
            # df = pd.read_csv(output_csv_name)
            # =============================================================================
            
            #==============================================================================
            # Plot
            #==============================================================================
            
            # get the length of fasta to set x-axis limit
            ###for record in SeqIO.parse(Input_GenBank_File, 'genbank'):
            ###    print('record: ', record)
            ###    fasta_length = len(record.seq)
            
            # Arrange Data to start and end from 0
            ###position.insert(0,0) # start of ORF = 685 (684 in python)
            ###position.append(fasta_length) # end = 18485 (18484 in python)
            ###depth.insert(0,0)
            ###depth.append(0)
            ###frequency.insert(0,1)
            ###frequency.append(1)
            ###ambiguity.insert(0,1)
            ###ambiguity.append(1)

            if 'Position' not in summary_data_frame:
                summary_data_frame['Position'] = position
            summary_data_frame['Frequency ' + subfolder] = frequency

            if not cumulative_frequency:
                cumulative_frequency = frequency.copy()
            else:
                for i in range(len(cumulative_frequency)):
                    cumulative_frequency[i] += frequency[i]

            total_clones += 1
            
            #==============================================================================
            # Plot line chart 
            #==============================================================================
                        
            plt.rcParams["figure.figsize"] = [4,1]
            plt.rcParams["figure.dpi"] = 300
            
            fig,(ax1,ax2) = plt.subplots(2,sharex = True)
            
            ax1.plot(position, depth, linestyle = "-", color = "blue", linewidth = 0.5, markerfacecolor="none")
            ax1.axis([0, fasta_length, 0, max(depth)])    # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene
            
            ax2.plot(position, frequency,"s", linestyle = "-", color = "red", linewidth = 0.5, markerfacecolor="none", markeredgecolor = "none")  # Draw line connecting dots
            #ax2.plot(position, frequency,"s", ms = 0.5, markerfacecolor="red", markeredgecolor = "none")  # Draw red dots then black dots 
            ax2.axis([0, fasta_length, -0.1, 1.1])  # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene 
            
            # Show and save image as
            fig.savefig(subfolder + '/' + Filename + '_' + subfolder + '.svg', dpi=300)
        else:
            print('Could not process data')
            
        # Print processing time
#Filter the false positive
replace_value=0.0
freq_threshold=1/total_clones
summary_data_frame['Frequency Average'] = [replace_value if i/total_clones <freq_threshold else i/total_clones for i in cumulative_frequency]
df_summary = pd.DataFrame.from_dict(summary_data_frame)
df_summary.to_csv(Filename + '_summary.csv')

target_codons_positions = summary_data_frame['Position']

plt.rcParams["figure.figsize"] = [4,1]
plt.rcParams["figure.dpi"] = 300

fig,(ax1,ax2) = plt.subplots(2,sharex = True)

for item in target_codons_positions[1: -1]:
    ax1.axvline(x = (item+1), c = 'red', linewidth = 0.15)

ax1.axis([0, fasta_length, -0.1, 1.1])    # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene
ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)

ax2.plot(target_codons_positions, summary_data_frame['Frequency Average'], "s", linestyle = "-", color = "red", linewidth = 0.5, markerfacecolor="none", markeredgecolor = "none")  # Draw line connecting dots
ax2.axis([1,fasta_length,-0.1,1.1])  # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene 

# Show and save image as
fig.savefig(Filename + '_summary.svg', dpi=300)


end_time = time.time() - start_time
print('Total processing time = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end_time))))


