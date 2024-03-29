#!/usr/bin/env python

#@TS 2017
#@AC 9/4/2018
#@JF 29/8/2018
#@DC updated in July 2020
#LF: Run alignment

#==============================================================================
# New analysis by codon rather than nucleotide. Align to wt-template containing 
# potential synthetic inserts. Parse through each individual good reads.
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

current_directory = os.getcwd()
current_directory, folders, files = next(os.walk(current_directory))

Input_GenBank_File = glob.glob('*.gb')[0]
switch = 1 #set to 1 if the genbank file is recoded and 0 if wt.

Filename = os.path.splitext(Input_GenBank_File)[0]

with open(Input_GenBank_File, "r") as input_handle:
    with open(Filename+".fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")        
        SeqIO.write(sequences, output_handle, "fasta")

fasta_length = len(next(SeqIO.parse(Input_GenBank_File, 'genbank')).seq)

###Input_Bam_File     = '/Users/wesr/Desktop/Landscapes/Ala_schemes/vA5/A5_6/vA5_20kb_REXER_MDS42_only.fa_sorted.bam'
###ref_fasta          = '/Users/wesr/Desktop/Landscapes/Ala_schemes/vA5/A5_6/vA5_20kb_REXER_MDS42_only.fa' 

# This is for Ala scheme 8
recoding_scheme = [{'target': 'TCG', 'recode': 'AGC'},
                   {'target': 'TCA', 'recode': 'AGT'},
                   {'target': 'TAG', 'recode': 'TAA'}]

# NGS filter
Min_Read_Length = 10 # filter out read with read length lower than indicated
MAPQ = 30 # filter out reads with mapping quality less than indicated. 0 if problems with synthetic insertion mismatching, higher to remove poor reads.

#Smoothing and artifact removal
Name_SI = 's.i.' #root name of the synthetic insertions in the genbank file, e.g. for s.i.1, s.i.2 root is s.i.
Artifact_removal_count = 2 #how many codons before the SI to not plot - adjust till artifacts vanish
Min_cov_threshold = 10 # if coverage below this do not calculate frequency, instead use last frequency

# =============================================================================
# Define output and functions
# =============================================================================
###Fasta_Filename = os.path.splitext(os.path.basename(ref_fasta))[0] # isolate filename           
###output_image_filename = ref_fasta + '.svg'
###output_csv_name = ref_fasta + '.csv'

# get a list of codon properties: 0 = codon, 1 = codon start, 2 = codon end, 3 = codon position set; 4 = strand ('f' or 'r'), 5= new genbank annotation
def get_codon(genbank, target_codon, recode_scheme, swap):  
    codon_pos_list = []
    if type(target_codon) != list:
        target_codon = [target_codon] 
    for feature in genbank.features:
        if feature.type == "CDS":
            start=feature.location.start.position
            end=feature.location.end.position
            sense=feature.strand   
            if sense == 1:                                                  
                orf = genbank.seq[start:end].upper()   
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]     
                if swap < 1:
                    codon_pos = [{'codon': str(x),
                                  'start': (i*3)+start,
                                  'end': (i*3)+start+3,
                                  'position': sorted(list(range((i*3)+start, (i*3)+start+3))), # set(range(start, end))
                                  'strand': 'f',
                                  'recode':[recode_scheme[a].get('recode') for a in range(len(recode_scheme)) if str(x) == recode_scheme[a].get('target')][0]} for i, x in enumerate(codon) if x in target_codon]     
                else:
                    codon_pos = [{'codon': str(x),
                                  'start': (i*3)+start,
                                  'end': (i*3)+start+3,
                                  'position': sorted(list(range((i*3)+start, (i*3)+start+3))), # set(range(start, end))
                                  'strand': 'f',
                                  'recode':[recode_scheme[a].get('target') for a in range(len(recode_scheme)) if str(x) == recode_scheme[a].get('recode')][0]} for i, x in enumerate(codon) if x in target_codon]     
            else:
                orf = genbank.seq[start:end].reverse_complement().upper()              
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]      
                if swap < 1:                    
                    codon_pos = [{'codon': str(x),
                                  'start': ((len(codon)-i)*3)+start,
                                  'end': ((len(codon)-i)*3)+start-3,
                                  'position': sorted(list(range(((len(codon)-i)*3)+start-3, ((len(codon)-i)*3)+start))), 
                                  'strand': 'r',
                                  'recode':[recode_scheme[a].get('recode') for a in range(len(recode_scheme)) if str(x) == recode_scheme[a].get('target')][0]} for i, x in enumerate(codon) if x in target_codon] 
                else:                    
                    codon_pos = [{'codon': str(x),
                                  'start': ((len(codon)-i)*3)+start,
                                  'end': ((len(codon)-i)*3)+start-3,
                                  'position': sorted(list(range(((len(codon)-i)*3)+start-3, ((len(codon)-i)*3)+start))), 
                                  'strand': 'r',
                                  'recode':[recode_scheme[a].get('target') for a in range(len(recode_scheme)) if str(x) == recode_scheme[a].get('recode')][0]} for i, x in enumerate(codon) if x in target_codon] 
            codon_pos_list.append(codon_pos)    
        elif feature.type == "misc_feature":
            if 'label' in feature.qualifiers:
                if (feature.qualifiers['label'][0]).find(Name_SI) > -1:
                    last_element = codon_pos_list[-1]
                    count = 1
                    while (count < Artifact_removal_count) and (len(last_element) > 1):
                        count = count + 1
                        last_element = last_element[:-1]
                    if Artifact_removal_count > 0:
                        if len(last_element) > 1:
                            last_element = last_element[:-1]
                            codon_pos_list[-1] = last_element
                        else:
                            codon_pos_list = codon_pos_list[:-1]
    return sorted(list(itertools.chain.from_iterable(codon_pos_list)), key = itemgetter('start'))

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
    TargetCodonseqrc = Seq(TargetCodon).reverse_complement()
    RecodedCodonseqrc = Seq(RecodedCodon).reverse_complement()
    codon_record = []
    with pysam.AlignmentFile(InputBamFile,"rb") as bamfile: 
        
        for read in bamfile.fetch(FastaFilename, CodonPosition[0], CodonPosition[-1]):  # take only the seq that contain targeted region, make sure targeted region is arranged by order
            try:  # file 1 sequence number 9 raised ValueError
                sequence = read.query_alignment_sequence # string
                position = read.get_reference_positions() # list
                quality = read.mapping_quality
                length = read.reference_length
            except ValueError:
                pass   
                     
            if len(position) == len(sequence) and length > MinReadLength: # no indel & correct position and min length = 200
                if quality > mapq and len(position) > 0: 
                    try:     
                        target_seq = [sequence[position.index(i)].upper() for i in CodonPosition] # exact index for targeted position for each read 
                        seq = ''.join(target_seq) # Eg ['T','C','G'] are joined as 'TCG' 
                        if strand is 'f':
                            if seq in [TargetCodon,RecodedCodon]: # check that the seq is either nrc or rc, else is ngs errors and discarded
                                codon_record.append(seq)
                        else:  #reverse complements on the reverse strand
                            if seq in [TargetCodonseqrc,RecodedCodonseqrc]: # check that the seq is either nrc or rc, else is ngs errors and discarded
                                codon_record.append(seq)   
                                #print("A")
                    except ValueError:  # ValueError raised when there is a deletion in the sequence and targetlist position is not present
                        pass
            

        if swap < 1:      
            try: 
                TotalRecoded = codon_record.count(RecodedCodon) + codon_record.count(RecodedCodonseqrc)
                freq = round((float(TotalRecoded)/len(codon_record)),2)
                ambig = round((float(TotalRecoded)/len(codon_record)),2)
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
                TotalRecoded = codon_record.count(TargetCodon) + codon_record.count(TargetCodonseqrc)
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
# Calculate codon frequency
#==============================================================================

mds = SeqIO.read(Input_GenBank_File, "genbank")

if switch < 1:      
    target_codon = [i['target'] for i in recoding_scheme]   
else:
    target_codon = [i['recode'] for i in recoding_scheme]

codon_list = get_codon(mds, target_codon, recoding_scheme, switch)

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
                                           swap = switch) for i in codon_list]
            
          
        df = pd.DataFrame.from_dict(codon_frequency)
        if 'Position' in df: # and 'Depth' in df and 'Frequency' in df:
            position  = [0] + list(df['Position'])  + [fasta_length] #adjusts for python 0 indexing and appends total legth
            depth     = [0] + list(df['Depth'])     + [0]
            frequency = [1] + list(df['Frequency']) + [1]
            ambiguity = [1] + list(df['Ambiguity']) + [1]
        
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

summary_data_frame['Frequency Average'] = [i/total_clones for i in cumulative_frequency]
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