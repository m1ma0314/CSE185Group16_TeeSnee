#!/usr/bin/env python

"""
Command-line script to run t-SNE
"""

import argparse
import pandas as pd

def main():
  parser = argparse.ArgumentParser(
      prog="tsne.py",
      description="command-line script to generate tsne plots"
  )

  # Input
  parser.add_argument("filename", help="gene data file (specify if zipped or not)",type=str)
  parser.add_argument("-p","--target_perplexity",help="user specificed perplexity",type=int,metavar="PERPLEXITY",required=False)
  parser.add_argument("-z","--zipped",help="unzip file if input is zipped",action="store_true")

  args = parser.parse_args()
  if args.zipped:
      print("File is zipped. Extracting and reading as CSV:", args.filename)
      unzipped_filename = args.filename
      unzipped_filename = pd.read_csv(unzipped_filename,compression='gzip', delimiter='\t')
      tsne(unzipped_filename,args.target_perplexity)
  else:
      tsne(args.filename,args.target_perplexity)


if __name__ == "__main__":
  main()
