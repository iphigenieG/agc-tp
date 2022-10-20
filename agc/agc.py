#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
#import statistics
import textwrap
#from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Iphigenie Gonnet"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Iphigenie Gonnet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Iphigenie Gonnet"
__email__ = "iphigenie.g@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = f"{path} is a directory"
        else:
            msg = f"{path} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     f"{sys.argv[0]} -h")
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int,
                        default = 400,
                        help="Minimum sequence length for dereplication"
                        "(default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int,
                        default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int,
                        default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int,
    default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """read_fasta
    Builds sequence generator from compressed fasta file

    Parameters
    ----------
    amplicon_file : str
        path to the .fasta.gz file containing sequences
    minseqlen : int
        minimum sequence length

    Yields
    ------
    str
        sequences from the fasta file longer than minseqlen
    """
    seq = ""
    with gzip.open(amplicon_file, "rt") as file:
        for line in file:
            if line.startswith('>'):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq+=line.strip()
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """dereplication_fulllength
    Builds generator of unique seqs in reverse order of occurrence from an
    amplicon file (compressed .fasta.gz file)


    Parameters
    ----------
    amplicon_file : str
        path to the .fasta.gz file containing sequences
    minseqlen : int
        minimum sequence length
    mincount : int
        minimum number of occurrences
    Yields
    ------
    list
        two item lists with first item the sequence and second its occurrence
    """
    dic = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        dic[seq] = dic.get(seq,0)+1
    sorted_seqs = sorted(dic, key = dic.get, reverse=True)
    for seq in sorted_seqs :
        if dic[seq]>= mincount :
            yield [seq,dic[seq]]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    seq1 = alignment_list[0]
    seq2 = alignment_list[1]
    nb_id = 0
    for i, nucleotide in enumerate(seq1):
        if nucleotide == seq2[i]:
            nb_id += 1
    return nb_id/len(seq1)*100

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount,
                                chunk_size, kmer_size):
    """abundance_greedy_clustering
    performs greedy clustering algorithm on an amplicon file, returns OTU list

    Parameters
    ----------
    amplicon_file : str
        path to the .fasta.gz file containing sequences
    minseqlen : int
        minimum sequence length
    mincount : int
        minimum number of occurrences
    chunk_size, kmer_size : unused

    Returns
    -------
    list
        list of OTUs with their number of occurrences
    """
    gen = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    list_OTU = []
    list_OTU.append(next(gen))
    for seq in gen :
        flag = True
        for otu in list_OTU :
            alignment_list = nw.global_align(seq[0], otu[0],
                                            gap_open=-1, gap_extend=-1,
                                            matrix=os.path.abspath(
                                            os.path.join(
                                            os.path.dirname(__file__),
                                            "MATCH")))
            if get_identity(list(alignment_list)) >= 97:
                otu[1] += seq[1]
                flag = False
                break
        if flag :
            list_OTU.append(seq)
    return list_OTU

def write_OTU(OTU_list, output_file):
    """write_OTU writes fasta file from OTU list

    Parameters
    ----------
    OTU_list : list
        list of OTUs whith each OTU a 2 item list with first item the seq and 
        second item the number of occurences of the OTU
    output_file : str
        path to the file to write to
    """
    with open(output_file,"w") as file :
        for i,otu in enumerate(OTU_list):
            sequence = otu[0]
            file.write(f">OTU_{i+1} occurrence:{otu[1]}\n")
            file.write(textwrap.fill(sequence, width=80))
            file.write("\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU_list = abundance_greedy_clustering(args.amplicon_file,args.minseqlen,
                                            args.mincount,args.chunk_size,
                                            args.kmer_size)
    write_OTU(OTU_list,args.output_file)

#==============================================================
# Chimera removal section
#==============================================================

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):
    """Split sequences in a least 4 chunks
    """
    pass

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


if __name__ == '__main__':
    main()
