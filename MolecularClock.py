#!/usr/bin/env python
import pandas as pd
from Bio import SeqIO
import numpy as np
from itertools import permutations
def FromFASTAtoDataFrame(input_file):
    with open(input_file) as fasta_file: 
        identifiers = []
        lengths = []
        seq =[]
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  
            identifiers.append(seq_record.id)
            lengths.append(len(seq_record.seq))
            seq.append(seq_record.seq)
    s1 = pd.Series(identifiers, name='ID')
    s2 = pd.Series(lengths, name='length')
    s3 = pd.Series(seq, name='sequence')
    df = pd.DataFrame(dict(ID=s1, length=s2, sequence=s3))
    return df
def HammingDistance(seq1,seq2):
    if len(seq1)==len(seq2):
        dist=sum(nucleotide1 != nucleotide2 for nucleotide1, nucleotide2 in zip(seq1, seq2)) #wikipedia
        return dist

def LevenshteinDistance(seq1, seq2):
    matriz=np.empty((len(seq1)+1, len(seq2)+1))
    for i in range(1, len(seq1)+1):
        matriz[i, 0]=i
    for j in range(1, len(seq2)+1):
        matriz[0, j]=j
    for j in range(1, len(seq2)+1):
        for i in range(1, len(seq1)+1):
            if seq1[i-1]==seq2[j-1]:
                cost=0
            else:
                cost=1
            matriz[i, j]= min(matriz[i-1, j]+1,matriz[i, j-1]+1, matriz[i-1, j-1]+cost)
    #print(matriz)
    return int(matriz[len(seq1), len(seq2)])
def matrizDistancias(df):
    matriz=np.empty((len(df), len(df)))
    for index, row in df.iterrows():
        seq1=row.sequence
        for i in range(len(df)):
                seq2=df.loc[i]['sequence']
                if len(seq1)==len(seq2):
                    dist=HammingDistance(seq1, seq2)
                else:
                    dist=LevenshteinDistance(seq1, seq2)
                matriz[index, i]=dist
    return matriz
if __name__ == "__main__":
    df=FromFASTAtoDataFrame("cytb_Birds.fasta")
    print(df.head())
    matriz=matrizDistancias(df.head())
    dfm = pd.DataFrame(matriz, columns=df.head().ID, index=df.head().ID)
    print(dfm)

