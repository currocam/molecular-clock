#!/usr/bin/env python
import pandas as pd
from Bio import SeqIO
import numpy as np


def FromFASTAtoDataFrame(input_file):
    """Procesar archivo .fasta con las secuencias de interés

    Parameters
    ----------
    input_file : str
        Dirección del archivo .fasta

    Returns
    -------
    df : DataFrame
        DataFrame con la información de las secuencias a estudiar.

    """
    with open(input_file) as fasta_file:
        identifiers = []
        lengths = []
        seq = []
        nombre = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            identifiers.append(seq_record.id)
            lengths.append(len(seq_record.seq))
            seq.append(seq_record.seq)
            # Nombre especie son 2 y 3 palabra de la descripción.
            nombre.append(' '.join(seq_record.description.split(' ')[1:3]))
    s1 = pd.Series(identifiers, name='ID')
    s2 = pd.Series(lengths, name='length')
    s3 = pd.Series(seq, name='sequence')
    s4 = pd.Series(nombre, name='nombre')
    df = pd.DataFrame(dict(ID=s1, length=s2, sequence=s3, nombre=s4))
    return df


def HammingDistance(seq1, seq2):
    """Calcula distancia Hamming entre dos secuencias, es decir, el mínimo nº de sustituciones para que coincidan ambas secuencias.

    Parameters
    ----------
    seq1 : str
        Secuencia 1 a estudiar.
    seq2 : str
        Secuencia 2 a estudiar.

    Returns
    -------
    dist : int
        Distancia de Hamming.

    """
    if len(seq1) == len(seq2):
        dist = sum(nucleotide1 != nucleotide2 for nucleotide1,
                   nucleotide2 in zip(seq1, seq2))  # wikipedia
        return dist


def LevenshteinDistance(seq1, seq2):
    """Calcula distancia Levenshtein entre dos secuencias, es decir, el mínimo nº de sustituciones, adiciones o deleciones para que ambas secuencias sean iguales.

    Parameters
    ----------
    seq1 : str
        Secuencia 1 a estudiar.
    seq2 : str
        Secuencia 2 a estudiar.

    Returns
    -------
    dist : int
        Distancia de Levenshtein.

    """
    matriz = np.empty((len(seq1) + 1, len(seq2) + 1))
    for i in range(1, len(seq1) + 1):
        matriz[i, 0] = i
    for j in range(1, len(seq2) + 1):
        matriz[0, j] = j
    for j in range(1, len(seq2) + 1):
        for i in range(1, len(seq1) + 1):
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            else:
                cost = 1
            matriz[i, j] = min(
                matriz[i - 1, j] + 1, matriz[i, j - 1] + 1, matriz[i - 1, j - 1] + cost)
    return int(matriz[len(seq1), len(seq2)])


def matrizDistancias(df):
    """Construye matriz con las distancias entre cada una se las secuencias. .

    Parameters
    ----------
    df : DataFrame
        DataFrame con la información de las secuencias a estudiar.
    Returns
    -------
    type
        Description of returned object.

    """
    matriz = np.empty((len(df), len(df))).astype(int)
    n = len(df)
    for i in range(n):
        seq1 = df.loc[i]['sequence']
        for j in range(n):
            seq2 = df.loc[j]['sequence']
            if matriz[j, i] != 0:
                matriz[i, j] = matriz[j, i]
            else:
                if len(seq1) == len(seq2):
                    dist = HammingDistance(seq1, seq2)
                else:
                    dist = LevenshteinDistance(seq1, seq2)
                matriz[i, j] = int(dist)
    print(type(matriz))
    return matriz


def matrizDistancias(df):
    """Construye matriz con las distancias entre cada una se las secuencias. .

    Parameters
    ----------
    df : DataFrame
        DataFrame con la información de las secuencias a estudiar.
    Returns
    -------
   matriz : numpy.ndarray
       Matriz con las distancias entre las distintas secuencias

    """
    matriz = np.zeros((len(df), len(df))).astype(int)
    n = len(df)
    for i in range(n):  # Para evitar calcular distancias innecesarias
        seq1 = df.loc[i]['sequence']
        for j in range(n):
            seq2 = df.loc[j]['sequence']
            if i == j:
                pass
            elif matriz[j, i] != 0:
                matriz[i, j] = matriz[j, i]
            else:
                if len(seq1) == len(seq2):
                    dist = HammingDistance(seq1, seq2)
                else:
                    dist = LevenshteinDistance(seq1, seq2)
                matriz[i, j] = int(dist)
    return matriz


if __name__ == "__main__":
    # cargamos secuencias del archivo.fasta
    df = MolecularClock.FromFASTAtoDataFrame("cytb_Birds.fasta")
    print("Los genes a estudiar son")
    print(df.head())

    # Generamos matriz distancias entre las distintas secuencias
    matriz = MolecularClock.matrizDistancias(df)
    print("la matriz con las cambios (sustituciones, inserciones o deleciones son: )")
    dfm = pd.DataFrame(matriz, columns=df.nombre, index=df.nombre)
    print(dfm)
