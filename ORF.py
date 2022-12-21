import csv
import os
from unittest import result
import pandas as pd
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import re
from textwrap import wrap


replacements = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y','TAA':'','TAG':'',
                'TGT':'C','TGC':'C','TGA':'','TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
                }
replacer = replacements.get  # For faster gets.

# only accept A,T,C,G in sequence
def FoolProof(seq):
    anti_atcg = re.compile('[^ATCG]')
    error = set(anti_atcg.findall(seq))
    return error

# python commond setting
def ArgSet():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help='Path of a input fasta')
    parser.add_argument("output", help='Enter the name of output csv')
    parser.add_argument('min_protein_len', default=20, nargs='?', const=1, type=int ,help='Minmum length(a.a.) of ORF (default = 20)')
    parser.add_argument('seq_repeat', default=4, nargs='?', const=1, type=int ,help='How many times sequence repeat to find ORF (default = 4 times)')

    args = parser.parse_args()
    return args

# ORF analyze
def main(min_orf_length,seq_repeat):
    orf_data = []
    # get iso and seq in fasta file 
    for iso in SeqIO.parse(input_file, "fasta"):
        iso_name = iso.id
        result = str(iso.seq)

        if len(FoolProof(result)) != 0:
            print('Sequence only accept A,T,C,G')
            print('"{}" include "{}" are not available.'.format(iso_name,'","'.join(FoolProof(result))))
            continue
        else:
            pass

        oorf = []
        # seq 頭尾相接四次用於尋找ORF
        result_4times = result*seq_repeat
        # seq 尾加上seq 開頭兩nucl.，尋找ATG
        circle_seq = result + result[0] + result[1]
        ATG_pos = [m.start() for m in re.finditer('ATG', circle_seq)]

        maybe_ORF_pos = []
        bsj = len(result)
        for atg in ATG_pos:
            taa = 0
            tag = 0
            tga = 0
            codon_type = wrap(result_4times[atg:], 3)
            if len(codon_type[len(codon_type)-1]) != 3:
                codon_type = codon_type[:-1]

            stop_co = []
            try:
                TAA_pos = codon_type.index('TAA')
                stop_co.append(TAA_pos)
                taa = 1
            except:
                #如果未找到TAA TAG TGA則定義結尾在codon_type最後
                stop_co.append(len(codon_type)-1)
            
            try:
                TAG_pos = codon_type.index('TAG')
                stop_co.append(TAG_pos)
                tag = 1
            except:
                stop_co.append(len(codon_type)-1)

            try:
                TGA_pos = codon_type.index('TGA')
                stop_co.append(TGA_pos)
                tga = 1
            except:
                stop_co.append(len(codon_type)-1)
            #print(atg,stop_co)

            #取ORF位置
            if (taa + tag + tga != 0):
                orf_type = 'ORF'
            else:
                orf_type = 'cORF'
            ORF_start = atg
            ORF_end = (atg + min(stop_co)*3+3)-1
            bsj = len(result)
            #判斷長度大於輸入的min_orf_length a.a.(cORF必定成立，所以不須另外判斷)
            if min(stop_co) >= min_orf_length:
                if seq_repeat == 1:
                    maybe_ORF_pos.append([orf_type,ORF_start,ORF_end])
                    continue
                elif seq_repeat > 1:
                    pass
                #如果整條序列橫跨bsj
                if ((ORF_start<bsj) and (ORF_end>=bsj)):
                    maybe_ORF_pos.append([orf_type,ORF_start,ORF_end])
            
        #去redundant_相同ORF end代表相同reading frame且為同一條seq
        if len(maybe_ORF_pos) >= 2:
            orf_df = pd.DataFrame(maybe_ORF_pos,columns=['type','start','end'])
            orf_df['end_on_circle'] = orf_df.end.apply(lambda x: x % bsj)
            orf_df['ORF_len'] = orf_df['end'] - orf_df['start'] + 1        
            orf_df = orf_df.sort_values(by='ORF_len')
            orf_df = orf_df.drop_duplicates(subset='end_on_circle',keep='last')
            orf_df = orf_df.drop(['end_on_circle'],axis=1)
            orf_df = orf_df.drop(['ORF_len'],axis=1)
            orf_df = orf_df.sort_values(by='start')
            ORF_pos = orf_df.values.tolist()
            for iorf in ORF_pos:
                oorf.append([iorf[0],iorf[1],iorf[2]])
            
        elif len(maybe_ORF_pos) == 1:
            ORF_pos = maybe_ORF_pos[0]
            oorf.append([ORF_pos[0],ORF_pos[1],ORF_pos[2]])
            
        else:
            oorf = []
        # orf_data record isoform name, orf position of isoform, isoform sequence
        orf_data.append([iso_name,oorf,result])
    return orf_data

# ORF nucl. translate to peptide
def translate(seq,orf_type):
    codon = wrap(seq, 3)
    peptide = ''.join([replacer(n, n) for n in codon])
    pep_len = len(peptide)
    if orf_type == 'ORF':
        stop_codon = seq[-3]+seq[-2]+seq[-1]
    elif orf_type == 'cORF':
        peptide = peptide + '*'
        stop_codon = '*'
    return peptide, pep_len ,stop_codon

# determine reading frame of orf on isoform sequence
def readframe(orf_start):
    frame_test = orf_start % 3
    if frame_test == 0:
        frame = 'frame1'
    elif frame_test == 1:
        frame = 'frame2'
    elif frame_test == 2:
        frame = 'frame3'
    return frame
    
# using isoform name, orf position of isoform, isoform sequence to ORF dataframe
def ORF_result(orf_temp):
    result_data = []
    for io in orf_temp:
        iso_id = io[0]
        orf_pos = io[1]
        iso_seq = io[2]
        iso_len = len(iso_seq)
        iso_4seq = iso_seq*4
        if orf_pos == 0:
            continue
        else:
            for orfp in orf_pos:
                orf_count = len(orf_pos)
                orf_type = orfp[0]
                orf_start = orfp[1]+1
                orf_end = orfp[2]+1
                orf_len = orfp[2] - orfp[1] + 1
                reading_frame = readframe(orfp[1])
                peptide, peptide_len ,stop_codon = translate(iso_4seq[orfp[1]:orfp[2]+1],orfp[0])
                result_data.append([iso_id,iso_len,orf_count,orf_type,orf_start,orf_end,orf_len,peptide_len,reading_frame,peptide,stop_codon])
    
    result_df = pd.DataFrame(result_data,columns=['Isoform_ID','len','ORF_count','ORF_type','ORF_start','ORF_end','ORF_len(nt)','peptide_len(a.a.)','reading_frame','peptide','stop_codon'])
    return result_df

if __name__ == '__main__':
    args = ArgSet()
    input_file = args.input
    output_name = args.output
    min_orf_length = args.min_protein_len
    seq_repeat = args.seq_repeat

    orf_temp = main(min_orf_length,seq_repeat)
    result_df = ORF_result(orf_temp)

    print(result_df)
    output_dir = './ORF_output/'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    result_df.to_csv('./ORF_output/{}.csv'.format(output_name),index_label='index')

