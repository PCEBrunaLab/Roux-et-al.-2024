#! /usr/bin/env Python3
'''
Process barcoding experiments by looping over groups of FASTQ files

Each set of concordant lines in I1, I2, R1 and R2 represent a single molecule from a single cell.

The Cellecta barcode format is such that there are 14nt and 4nt spacer, 'TGGT', and a 30nt barcode. We assume 
that there will be some mutations in any of these 3 components, so using the whitelist will help enormously here.

File contents:
I1 - Sample dual index i5 - 10nts
I2 - Sample dual index i7 - 10 nts
R1 - 10X Cell barcode - 16nts + UMI - 12nts + any variable length TruSeq Read1
R1 - Cellecta barcode - Xnts + any variable length TruSeq Read2
'''

import argparse
import re
import os
import sys
import itertools
import logging
import numpy as np
import pandas as pd
import itertools
#Cloned from Github repository https://github.com/MikeDMorgan/Genomics
sys.path.append("/homes/mdmorgan/src/Genomics")
import Genomics.Classes.Fastq as Fastq
import Genomics.Classes.Kmer as Kmer
import multiprocessing
import functools
import math


def revComp(sequence):
    '''
    Return the reverse complement of an input sequence
    '''

    nt_dict = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}

    return("".join([nt_dict[NT] for NT in sequence])[::-1])


def computeEditDistance(seq1, seq2):
    '''
    Count the number of like-for-like differences between 2 strings of the same length
    '''

    if len(seq1) != len(seq2):
        raise TypeError("String lengths do not match")
    else:            
        return(sum([0 if seq1[ix] == seq2[ix] else 1 for ix in range(len(seq1))]))


def matchKmers(seq1, known_seqs, k=14, epsilon=0.1):
    '''
    Match a set of kmers to a list of known sequences,
    return the Kmers that have an edit distance/k < epsilon
    '''

    seq_kmers = Kmer.findKmers(seq1, k)
    best_edits = []
    best_kmer = []
    for kmer in seq_kmers:
        for seq2 in known_seqs:
            # window over the sequence to find matches
            known_kmers = np.array(list(Kmer.findKmers(seq2, k)), dtype=str)
            known_edits = np.empty(known_kmers.shape)

            for i in range(known_kmers.shape[0]):
                seq2_kmer = known_kmers[i]
                known_edits[i] = computeEditDistance(kmer, seq2_kmer)

            sort_edits = np.argsort(known_edits)
            sorted_kmer = known_kmers[sort_edits]

            if any(known_edits[sort_edits]/float(k) < epsilon):
                best_kmer.append(sorted_kmer[0])
                best_edits.append(sort_edits[0])

    return(np.array(best_edits), np.array(best_kmer))


def bestKmer(seq1, whitelist, krange, epsilon):
    '''
    Given a kmer length, find the nearest matching kmer
    '''

    longest_kmer = 0
    best_edit = None
    best_kmer = None
    
    for k in krange:
        edits, kmer = matchKmers(seq1, whitelist, k=k, epsilon=epsilon)
        
        if edits.shape[0] > 0:
            longest_kmer = k
            best_edit = edits.min()
            if kmer.shape[0] == 1:
                best_kmer = kmer[0]
            elif kmer.shape[0] > 1 and (k == 14 | k == 30):
                logging.warning("Multiple best kmers: {}".format(len(kmer)))

    return(best_kmer, best_edit)


def parseFileSet(file_dict):
    '''
    Parse a set of 4 FASTQ files: I1, I2, R1 and R2 
    '''

    out_dict = {}
    for fkey in file_dict.keys():        
        i1_file = [i1 for i1 in file_dict[fkey] if re.search("I1", i1)][0]
        i1_records = Fastq.Fastq(i1_file)
        
        i2_file = [i2 for i2 in file_dict[fkey] if re.search("I2", i2)][0]
        i2_records = Fastq.Fastq(i2_file)

        r1_file = [r1 for r1 in file_dict[fkey] if re.search("R1", r1)][0]
        r1_records = Fastq.Fastq(r1_file)

        r2_file = [r2 for r2 in file_dict[fkey] if re.search("R2", r2)][0]
        r2_records = Fastq.Fastq(r2_file)

        out_dict[fkey] = [i1_records, i2_records, r1_records, r2_records]

    return(out_dict)


## file handles can't be pickled, so need to extract the tuple of sequences first
## can this be parallelised??
def tupleSequences(fastq_set):
    '''
    Extract the tuple of sequences from I1, I2, R1, R2
    return a list of tuples
    '''

    fcount = 0 # count the number of matched lines so far
    is_next = True # need this to check if EOF is reached

    if len(fastq_set) != 4:
        raise IOError("Missing FASTQ file records - only {} present".format(len(fastq_set)))

    f1_iter, f2_iter, f3_iter, f4_iter = enumerate(fastq_set)
    f1_name = f1_iter[1].file_handle.name
    f2_name = f2_iter[1].file_handle.name
    f3_name = f3_iter[1].file_handle.name
    f4_name = f4_iter[1].file_handle.name

    while is_next:
        fcount += 1
        try:
            f1_rec = f1_iter[1].next()
            f2_rec = f2_iter[1].next()
            f3_rec = f3_iter[1].next()
            f4_rec = f4_iter[1].next()
        except StopIteration:
            is_next = False

        i1_seq = f1_rec.sequence
        i2_seq = f2_rec.sequence
        if fcount % 100000 == 0:
            logging.info("Scanned {} sequences".format(fcount))

        # make sure the reads have been processed properly
        if f1_rec.length != f2_rec.length:
            is_valid = False
        else:
            is_valid = True
            yield((f1_rec, f2_rec, f3_rec, f4_rec))


def extractFastqSeq(fastq_set, wl_14nt=None, wl_30nt=None):
    '''
    Given an input tuple of Fastq.Fastq objects - loop over all 4 files 
    simultaneously to extract the matched information across them.
    '''

    f1_rec, f2_rec, f3_rec, f4_rec = enumerate(fastq_set)
    
    # process the UMI and cell barcode sequences
    cb = f3_rec[1].sequence[:16]
    umi = f3_rec[1].sequence[16:]
        
    # process the celecta barcode sequences
    # find the closest matching 14nt sequence
    best_kmer, best_edit = bestKmer(f4_rec[1].sequence, whitelist=wl_14nt, krange=range(14, 15), epsilon=3/14.0)
    # try the rev comp
    is_rev = False
    if best_kmer is None:
        best_kmer, best_edit = bestKmer(revComp(f4_rec[1].sequence), whitelist=wl_14nt,
                                        krange=range(14, 15), epsilon=3/14.0)
        is_rev = True

    if best_kmer is not None:
        if is_rev:
            ksplit = re.compile(best_kmer).split(revComp(f4_rec[1].sequence))
        else:
            ksplit = re.compile(best_kmer).split(f4_rec[1].sequence)
            # the split is on the best-matching kmer
            # it does not return the actual splitting sequence
            
        if len(ksplit) == 1:
            valid_bc14 = False
            bc14 = None
            spacer = None
            bc30 = None
            rec_dict = None
        else:
            valid_bc14 = True
            bc14 = best_kmer
            spacer = ksplit[1][:4]
            bc30 = ksplit[1][4:34]
                
            # find the best 30mer now
            if len(bc30) < 30:
                rec_dict = None
            else:
                best_30mer, best_30edit = bestKmer(bc30, whitelist=wl_30nt, krange=range(30, 31), epsilon=6/30)
                is_30rev = False
                if best_30mer is None:
                    best_30mer, best_30edit = bestKmer(revComp(bc30), whitelist=wl_30nt, krange=range(30, 31), epsilon=6/30)
                    is_30rev = True

                if is_rev == is_30rev:
                    bc30 = best_30mer

                    rec_dict = {"CB":cb, "UMI":umi, "BC.14":bc14, "Valid.BC14":valid_bc14,
                                "Edit.14":best_edit, "BC.30":bc30, "Spacer":spacer}
                else:
                    rec_dict = None
    else:
        rec_dict = None

    return(rec_dict)                


def mapBarcodeProcess(seq_tuple, wlist_14nt, wlist_30nt):
    '''
    Convenience wrapper for mapping in chunks
    '''  

    ores = extractFastqSeq(seq_tuple, wl_14nt=wlist_14nt, wl_30nt=wlist_30nt)
    if ores is not None:
        return(ores)


def mapExtractSequences(fq_key, fq_dict):
    '''
    Convenience wrapper to extract sequences using map()
    '''

    return(tupleSequences(fq_dict[fq_key]))


def dropNone(element):
    '''
    Filter out elements that aren't dict and are None
    '''

    if element is None:
        return False
    else:
        return True
    

parser = argparse.ArgumentParser(description="Extract and process Cellecta barcodes from 10X Chromium libraries")
parser.add_argument("--I1", dest="I1_file", type=str,
                    help="The path to the index 1 file(s) for the input library - separated by commas for multiple files")

parser.add_argument("--I2", dest="I2_file", type=str,
                    help="The path to the index 2 file(s) for the input library - separated by commas for multiple files")

parser.add_argument("--R1", dest="R1_file", type=str,
                    help="The path to the read 1 file(s) for the input library - separated by commas for multiple files")

parser.add_argument("--R2", dest="R2_file", type=str,
                    help="The path to the read 2 file(s) for the input library - separated by commas for multiple files")

parser.add_argument("--whitelist-14nt", dest="wl_14nt", type=str,
                    help="The whitelist of known 14nt Cellecta barcodes", default=None)

parser.add_argument("--whitelist-30nt", dest="wl_30nt", type=str,
                    help="The whitelist of known 30nt Cellecta barcodes", default=None)

parser.add_argument("--output", dest="output", type=str,
                    help="The output file path to save barcode X cell barcode matrix", default=sys.stdout)

parser.add_argument("--processes", dest="proc_workers", type=int,
                    help="Number of processes/CPUs to use", default=4)

parser.add_argument("--chunk-size", dest="chunk_size", type=int,
                    help="Number of sequences to evalutate on each pool of workers at once", default=1000000)

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination", default=sys.stdout)

args = parser.parse_args()

# setup the logger   
if type(args.logfile) == str: 
    logging.basicConfig(level=logging.INFO,
                        filename=args.logfile) 
else:
    logging.basicConfig(level=logging.INFO)

logging.info("Reading in 14nt barcode whitelist: {}".format(args.wl_14nt))
with open(args.wl_14nt, "rt") as short_ifile:
    allowed_14nt = [bx.strip("\n") for bx in short_ifile.readlines()]
    
logging.info("Reading in 30nt barcode whitelist: {}".format(args.wl_30nt))
with open(args.wl_30nt, "rt") as short_30ifile:
    allowed_30nt = [bx.strip("\n") for bx in short_30ifile.readlines()]

# input are the 4 sets of FASTQ files - need to check names for matching lanes etc
# need to check they all have the same number of files

## NB: The Cellecta barcodes will be the reverse complement, i.e. BC30-ACCA-BC14
i1_files = args.I1_file.split(",")
i2_files = args.I2_file.split(",")
r1_files = args.R1_file.split(",")
r2_files = args.R2_file.split(",")

n_file = set([len(i1_files), len(i2_files), len(r1_files), len(r2_files)])
if len(n_file)  == 1:
    pass
else:
    raise IOError("Input FASTQ file numbers should match")

logging.info("Grouping FASTQ files by sequencing lane")
name_regex = re.compile("([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9]+)\.fastq\.gz")
if len(i1_files) == 1:
    name_form  = str(set(["_".join(name_regex.search(IF).group(3,4)) for IF in (i1_files[0], i2_files[0], r1_files[0], r2_files[0])]))
    file_dict = {name_form : (i1_files[0], i2_files[0], r1_files[0], r2_files[0])}
else:
    # flatten all of the files first
    file_flatten = [x for x in itertools.chain.from_iterable([i1_files, i2_files, r1_files, r2_files])]
    name_form = list(set(["_".join(name_regex.search(IF).group(3,4)) for IF in file_flatten]))
    file_dict = {}
    for nm in name_form:
        file_dict[nm] = tuple([fx for fx in file_flatten if re.search(nm, fx)])
             
logging.info("Parsing FASTQ file sets")
# now we have a dictionary for the set of 4 files for a given sequencing lane
fq_records = parseFileSet(file_dict) # dictionary containing generators of FASTQ file contents

## need to extract the reads first as file handles can't be pickled.
tup_func = functools.partial(mapExtractSequences, fq_dict=fq_records)
all_fqrecs = list(map(tup_func, fq_records.keys()))

all_rec = 0
invalid_rec = 0
big_dict = {}

logging.info("Processing sequences for cell barcodes, UMIs and Cellecta barcodes")
mapping_func = functools.partial(mapBarcodeProcess, wlist_14nt=allowed_14nt,
                                 wlist_30nt=allowed_30nt)

logging.info("Starting multiprocessing with {} CPUs".format(args.proc_workers))
work_pool = multiprocessing.Pool(args.proc_workers)

chunk_list = []
fcount = 0
chnk_n = 0
for FASTQ in all_fqrecs:
    # generator object contains all of the FASTQ records - this is what we need to chunk up!
    for fq in FASTQ:
        # stack up a list of sequence tuples here - then submit that iterator to the pool
        chunk_list.append(fq)
        chnk_n += 1

        if chnk_n == args.chunk_size:
            logging.info("Running sequence processing over {} sequences".format(chnk_n))
            chunk_res = work_pool.imap_unordered(mapping_func, chunk_list, chunksize=10000) # this is the bottleneck
            # filter out the None elements
            chunk_res = list(filter(dropNone, chunk_res))
            chunk_list = []
            chnk_n = 0

            for rq in chunk_res:
                all_rec += 1

                # check any None types
                if any([NT is None for NT in rq.values()]):
                    logging.info("None Type found in output dict - ignoring record")
                else:
                    if rq:
                        if all_rec % 1000 == 0:
                            logging.info("Found {} valid Cellecta barcodes with cell barcodes".format(all_rec))
                        bk = "_".join([rq["CB"], rq["UMI"], rq["BC.14"], rq["BC.30"]])

                        try:
                            big_dict[bk] += 1
                        except KeyError:
                            big_dict[bk] = 1
            logging.info("Found {} valid barcodes in {} sequences".format(all_rec, fcount))
            
        fcount += 1

# need to: 1) error correct cell barcodes, 2) error correct UMIs, 3) error correct barcodes, 4) count/check for
# valid cellecta barcodes
# this scale of this processing can be reduced by only considering fragments with valid barcodes
logging.info("Collating UMIs and Cellecta barcodes over {} valid fragments".format(len(big_dict)))
big_df = pd.DataFrame.from_dict(big_dict, orient='index', columns=['count'])
big_df.reset_index(inplace=True)

big_df["CB"] = big_df["index"].apply(lambda X: X.split("_")[0])
big_df["UMI"] = big_df["index"].apply(lambda X: X.split("_")[1])
big_df["BC.14"] = big_df["index"].apply(lambda X: X.split("_")[2])
big_df["BC.30"] = big_df["index"].apply(lambda X: X.split("_")[3])
big_df.drop(axis=1, labels="index", inplace=True)

# I'll leave the error correction for a different script
# always save compressed format
logging.info("Writing large data file to {}".format(args.output))
# compresslevel : 1 gives the fastest compression - worth balancing due to the large amount of data here?
# compresslevel : 6 apparently gives a good trade-off between speed and compression
big_df.to_csv(args.output, sep="\t", index=False, compression={'method':'gzip', 'compresslevel':6},
              chunksize=10000)

#raise IOError("stop here")

## The memory ends up exploding here from having to constantly grow the list
## I need a way to chunk this up better - can I run 1million reads and summarise those first, then do the rest?

# logging.info("Running sequence processing over {} pools".format(npools))
# all_res = list(itertools.chain.from_iterable([pool_list[PX].map(mapping_func, all_fqrecs[PX]) for PX in range(npools)]))

# for rq in all_res:
#     if rq["CB"] is None:
#         pass
#     else:    
#         all_rec += 1
#         if all_rec % 1000 == 0:
#             logging.info("Found {} valid Cellecta barcodes with cell barcodes".format(all_rec))
#             bk = "_".join([rq["CB"], rq["UMI"], rq["BC.14"], rq["BC.30"]])

#         try:
#             big_dict[bk] += 1
#         except KeyError:
#             big_dict[bk] = 1

# for fk in fq_records.keys():
#     logging.info("Processing data from {} FASTQ files".format(fk))
#     for rq in extractFastqSeq(fq_records[fk], allowed_14nt, allowed_30nt):
#         all_rec += 1
#         if all_rec % 1000 == 0:
#             logging.info("Found {} valid Cellecta barcodes with cell barcodes".format(all_rec))
#         bk = "_".join([rq["CB"], rq["UMI"], rq["BC.14"], rq["BC.30"]])

#         try:
#             big_dict[bk] += 1
#         except KeyError:
#             big_dict[bk] = 1





# logging.info("Processing sequences for cell barcodes, UMIs and Cellecta barcodes")
# mapping_func = functools.partial(mapBarcodeProcess, mydict=fq_records,
#                                  wlist_14nt=allowed_14nt, wlist_30nt=allowed_30nt)

# logging.info("Starting multiprocessing with {} CPUs".format(args.proc_workers))
# wpool = multiprocessing.Pool(args.proc_workers)
# fq_keys = list(fq_records.keys())
# all_res = wpool.map(mapping_func, fq_keys)

# for rq in all_res:
#     all_rec += 1
#     if all_rec % 1000 == 0:
#         logging.info("Found {} valid Cellecta barcodes with cell barcodes".format(all_rec))
#     bk = "_".join([rq["CB"], rq["UMI"], rq["BC.14"], rq["BC.30"]])

#     try:
#         big_dict[bk] += 1
#     except KeyError:
#         big_dict[bk] = 1

# for fk in fq_records.keys():
#     logging.info("Processing data from {} FASTQ files".format(fk))
#     for rq in extractFastqSeq(fq_records[fk], allowed_14nt, allowed_30nt):
#         all_rec += 1
#         if all_rec % 1000 == 0:
#             logging.info("Found {} valid Cellecta barcodes with cell barcodes".format(all_rec))
#         bk = "_".join([rq["CB"], rq["UMI"], rq["BC.14"], rq["BC.30"]])

#         try:
#             big_dict[bk] += 1
#         except KeyError:
#             big_dict[bk] = 1

# a list containing dicts of each read set with CB, UMI and barcodes - includes invalid barcodes
