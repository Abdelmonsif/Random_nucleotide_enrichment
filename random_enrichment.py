#!/usr/bin/env python3
import numpy as np
import os,sys,argparse
import pandas as pd
import time
import glob
import matplotlib.pyplot as plt
from functools import partial
import random as rand

def Nt_real(SamFile, features, reall):
    df100 = pd.read_csv(SamFile, delim_whitespace=True, names=['Qname', 'flag', 'transcript_id', 'start', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'lol', 'a7a'], usecols=['Qname', 'transcript_id', 'start', 'SEQ'])
    df100['length'] = df100['SEQ'].str.len()
    df100['end'] = df100['start'] + df100['length']
    df2 = pd.read_csv(args.features, delim_whitespace=True)
    df10 = df100.merge(df2, on=['transcript_id'])
    df10['transcript_length'] = df10['utr5_length'] + df10['cds_length'] + df10['utr3_length']
    print(len(df10))
    df10['G'] = df10['SEQ'].str.count('G')
    df10['C'] = df10['SEQ'].str.count('C')
    df10['A'] = df10['SEQ'].str.count('A')
    df10['T'] = df10['SEQ'].str.count('T')

    G = df10['G'].sum()
    C = df10['C'].sum()
    A = df10['A'].sum()
    T = df10['T'].sum()
    total = G + C + A + T
    real = pd.DataFrame(columns=['G', 'C', 'A', 'T'])
    real = real.append({'G': G, 'C': C, 'A': A,'T': T}, ignore_index=True)
    real.to_csv(args.reall, index=False, sep='\t')

def Nt_random(SamFile, features, transcripts, randomization, RandomFolder):
    df100 = pd.read_csv(SamFile, delim_whitespace=True, names=['Qname', 'flag', 'transcript_id', 'start', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'lol', 'a7a'], usecols=['Qname', 'transcript_id', 'start', 'SEQ'])
    df100['length'] = df100['SEQ'].str.len()
    df100['end'] = df100['start'] + df100['length']
    df2 = pd.read_csv(features, delim_whitespace=True)
    df10 = df100.merge(df2, on=['transcript_id'])
    df10['transcript_length'] = df10['utr5_length'] + df10['cds_length'] + df10['utr3_length']
    seq =[]
    random_tr = []
    with open(transcripts) as f:
        for line in f:
            if line.startswith('>'):
                hobba = next(f)
                seq.append(hobba)
                zobry = line.split(' ')[0].strip('>')
                random_tr.append(zobry)
    all_tr = pd.DataFrame(seq, columns=['transcript_seq'])
    all_tr['transcript_id'] = random_tr
    df_final = df10.merge(all_tr, on=['transcript_id'])
    df_final['Mid'] = (df_final['start'] + df_final['end'])/2 
    random_times = pd.DataFrame(columns=['G', 'C', 'A', 'T'])
    for i in range(int(randomization)):
        random_start = []
        random_seq = []
        for index, row in df_final.iterrows():
            a7a = []
            z = row['utr5_length'] + row['cds_length']
            middle = int(row['Mid'])
            x = row['length']
            ur5 = int(row['utr5_length'])
            cd = int(row['cds_length'] + ur5)
            ur3 = int(row['transcript_length'] - x)
            if middle <= row['utr5_length']:
                if x <= row['utr5_length']:   
                    row['randomm'] = rand.randint(0, ur5)
                    a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                    random_seq.append(a7a)
                elif x > row['utr5_length']: 
                    row['randomm'] = rand.randint(0, cd)
                    a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                    random_seq.append(a7a)
            elif middle > row['utr5_length']:
                if middle <= int(z):
                    if x <= row['cds_length']:
                        row['randomm'] = rand.randint(int(row['utr5_length']), cd)
                        a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                        random_seq.append(a7a)
                    elif x > row['cds_length']:
                        row['randomm'] = rand.randint(0, cd)
                        a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                        random_seq.append(a7a)
                elif middle > z:
                    if x <= row['utr3_length']:
                        row['randomm'] = rand.randint(int(z), ur3)
                        a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                        random_seq.append(a7a)
                    elif x > row['utr3_length']:
                        row['randomm'] = rand.randint(int(row['utr5_length']), ur3)
                        a7a = row['transcript_seq'][row['randomm']:(row['randomm'] + row['length'])]
                        random_seq.append(a7a)
        random_df = pd.DataFrame(random_seq, columns=['random_seq'])
        random_df = random_df.apply(lambda x: x.astype(str).str.upper())
        random_df['G'] = random_df['random_seq'].str.count('G')
        random_df['C'] = random_df['random_seq'].str.count('C')
        random_df['A'] = random_df['random_seq'].str.count('A')
        random_df['T'] = random_df['random_seq'].str.count('T')
        random_G = random_df['G'].sum()
        random_C = random_df['C'].sum()
        random_A = random_df['A'].sum()
        random_T = random_df['T'].sum()
        random_total = random_G + random_C + random_A + random_T
        random_times = pd.DataFrame(columns=['G', 'C', 'A', 'T'])
        random_times = random_times.append({'G': random_G, 'C': random_C, 'A': random_A,'T': random_T}, ignore_index=True)
        random_times.to_csv(RandomFolder + '/' + str(i), index=False, sep='\t')
def plotting(reall, RandomFolder, title, sora):
    real = pd.read_csv(reall, delim_whitespace=True)
    dfs=[]
    all_files = glob.glob(RandomFolder + "*")
    for filename in all_files:
        df = pd.read_csv(filename, delim_whitespace=True)
        dfs.append(df)
        frame = pd.concat(dfs, axis=0, ignore_index=True)
    random = pd.DataFrame(columns=['G', 'C', 'A', 'T'])
    random_G = frame['G'].mean()
    random_C = frame['C'].mean()
    random_A = frame['A'].mean()
    random_T = frame['T'].mean()
    random = random.append({'G': random_G, 'C': random_C, 'A': random_A,'T': random_T}, ignore_index=True)
    G_more = 0
    G_less = 0
    C_more = 0
    C_less = 0
    A_more = 0
    A_less = 0
    T_more = 0
    T_less = 0
    for index, row in frame.iterrows():
        x = row['G'] - real['G']
        if x[0] >= 0:
            G_more += 1
        else:
            G_less += 1
        y = row['C'] - real['C']
        if y[0] >= 0:
            C_more += 1
        else:
            C_less += 1
        z = row['A'] - real['A']
        if z[0] >= 0:
            A_more += 1
        else:
            A_less += 1
        n = row['T'] - real['T']
        if n[0] >= 0:
            T_more += 1
        else:
            T_less += 1 
    significance = pd.DataFrame(columns=['A_more', 'A_less', 'C_more', 'C_less', 'T_more', 'T_less', 'G_more', 'G_less'])
    significance = significance.append({'A_more': A_more/len(frame), 'A_less': A_less/len(frame), 'C_more': C_more/len(frame),'C_less': C_less/len(frame), 'T_more': T_more/len(frame), 'T_less': T_less/len(frame), 'G_more': G_more/len(frame), 'G_less': G_less/len(frame)}, ignore_index=True)
    if significance['A_more'][0] <= 0.001 or significance['A_more'][0] >=0.999:
        a7a_A = '***'
    elif significance['A_more'][0] <= 0.01 or significance['A_more'][0] >=0.99:
        a7a_A = '**'
    elif significance['A_more'][0] <= 0.05 or significance['A_more'][0] >= 0.95:
        a7a_A = '*'   
    elif significance['A_more'][0] >= 0.05 or significance['A_more'][0] <= 0.95:
        a7a_A = 'ns'

    if significance['C_more'][0] <= 0.001 or significance['C_more'][0] >=0.999:
        a7a_C = '***'
    elif significance['C_more'][0] <= 0.01 or significance['C_more'][0] >=0.99:
        a7a_C = '**'
    elif significance['C_more'][0] <= 0.05 or significance['C_more'][0] >= 0.95:
        a7a_C = '*'   
    elif significance['C_more'][0] >= 0.05 or significance['C_more'][0] <= 0.95:
        a7a_C = 'ns'

    if significance['T_more'][0] <= 0.001 or significance['T_more'][0] >=0.999:
        a7a_T = '***'
    elif significance['T_more'][0] <= 0.01 or significance['T_more'][0] >=0.99:
        a7a_T = '**'
    elif significance['T_more'][0] <= 0.05 or significance['T_more'][0] >= 0.95:
        a7a_T = '*'   
    elif significance['T_more'][0] >= 0.05 or significance['T_more'][0] <= 0.95:
        a7a_T = 'ns'

    if significance['G_more'][0] <= 0.001 or significance['G_more'][0] >=0.999:
        a7a_G = '***'
    elif significance['G_more'][0] <= 0.01 or significance['G_more'][0] >=0.99:
        a7a_G = '**'
    elif significance['G_more'][0] <= 0.05 or significance['G_more'][0] >= 0.95:
        a7a_G = '*'   
    elif significance['G_more'][0] >= 0.05 or significance['G_more'][0] <= 0.95:
        a7a_G = 'ns'
    ratio = real/random
    ratio1 = ratio[['A', 'C', 'T', 'G']]
    lol = ratio1.T
    ax = lol.plot(kind='bar', rot=360, figsize=(12,8), color='grey')
    plt.axhline(y=1, color='r', linestyle='-.')
    ax.legend([])
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    plt.xlabel('', fontsize=16)
    plt.text(-0.05, 1.04, a7a_A, fontsize=18)
    plt.text(0.95, 1.03, a7a_C, fontsize=18)
    plt.text(1.95, 0.92, a7a_T, fontsize=18)
    plt.text(2.95, 1.04, a7a_G, fontsize=18)
    plt.ylabel('Reads enrichment', fontsize=22)
    ax.set_title(title + "\n",fontsize=24)

    fig = ax.get_figure()
    fig.savefig(sora, dpi=100)        
def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-SamFile', required=True)
    parser.add_argument('-features', required=True)
    parser.add_argument('-reall', required=True)
    parser.add_argument('-RandomFolder', required=True)
    parser.add_argument('-transcripts', required=True)
    parser.add_argument('-title', required=True)
    parser.add_argument('-sora', required=True)
    parser.add_argument('-randomization', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    metnaka = Nt_real(args.SamFile, args.features, args.reall)
    sharmoota = Nt_random(args.SamFile, args.features, args.transcripts, args.randomization, args.RandomFolder)
    wes5a = plotting(args.reall, args.RandomFolder, args.title, args.sora)
    end = time.time()
    print ('time elapsed:' + str(end - start))
