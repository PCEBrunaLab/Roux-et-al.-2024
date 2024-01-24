#! /usr/bin/env Python3

'''
Run over the cellbarcodes and compress them down to the unique observed barcodes
'''

import argparse
import re
import os
import sys
import logging
import pandas as pd

parser = argparse.ArgumentParser(description="Compress Cell barcode mapping to just unique observations")
parser.add_argument("--file", dest="infile", type=str,
                    help="Input gzipped file to process")

parser.add_argument("--output", dest="outfile", type=str,
                    help="Output file name")

parser.add_argument("--log", dest="logfile", type=str,
                    help="Logging file destination", default=sys.stdout)

args = parser.parse_args()
# setup the logger   
if type(args.logfile) == str: 
    logging.basicConfig(level=logging.INFO,
                        filename=args.logfile) 
else:
    logging.basicConfig(level=logging.INFO)

logging.info("Reading {}".format(args.infile))
in_df = pd.read_table(args.infile, sep="\t", compression="gzip", header=0, iterator=True, chunksize=10000)
out_list = []
xcount = 1
for x in in_df:
    # read chunk at a time
    logging.info("Reading chunk {}".format(xcount))
    x.loc[:, "Original"] = x["Original"].apply(lambda X: ",".join(set(X.split(","))))
    out_list.append(x)
    xcount += 1

out_df = pd.concat(out_list)
logging.info("Saving DF to {}".format(args.outfile))
out_df.to_csv(args.outfile, sep=",", compression="gzip", index=False)
    

