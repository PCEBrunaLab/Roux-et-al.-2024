#! /usr/bin/env Python3
'''
Extract the cellranger corrected barcodes from the alignment bam files - map to the observed/sequenced 
cellranger barcodes
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
import pysam

parser = argparse.ArgumentParser(description="Extract error-corrected cell barcodes from 10X Chromium libraries")
parser.add_argument("--bam", dest="bam_files", type=str,
                    help="A list of paths to the input .bam files")

parser.add_argument("--output", dest="output", type=str,
                    help="The output file path to save the mapping between corrected and observed/sequenced cell barcodes", default=sys.stdout)

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination", default=sys.stdout)

args = parser.parse_args()

# setup the logger   
if type(args.logfile) == str: 
    logging.basicConfig(level=logging.INFO,
                        filename=args.logfile) 
else:
    logging.basicConfig(level=logging.INFO)

nfiles = len(args.bam_files.split(","))
logging.info("Found {} bam files".format(nfiles))

bam_files = [pysam.AlignmentFile(bfx, "rb") for bfx in args.bam_files.split(",")]

bc_dict = {}
ncount = 0
nfile = 0
for bam_file in bam_files:
    logging.info("Extracting barcodes from {}".format(bam_file))
    for read in bam_file.fetch():
        try:
            try:
                xbcs = bc_dict[read.get_tag("CB")].split(",")
                xbcs.append(read.get_tag("CR"))
                
                bc_dict[read.get_tag("CB")] = ",".join(list(set(xbcs))
            except KeyError:
                bc_dict[read.get_tag("CB")] = read.get_tag("CR")
                ncount +=1
        except KeyError:
            pass
    
    bam_file.close()

bc_map = pd.DataFrame.from_dict(bc_dict, orient='index')
bc_map.columns = ["Original"]

bc_map.to_csv(args.output, sep="\t", compression="gzip", index_label="Corrected")

