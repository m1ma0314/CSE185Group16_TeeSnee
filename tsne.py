import pandas as pd
import numpy as np
import math

#read in cell by gene count matrix as pandas dataframe
count_matrix_df = pd.read_csv("GSE164073_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", delimiter = "\t") 

#convert dataframe (cell by gene count matrix) to 2D numpy array 
counts = count_matrix_df.to_numpy() 

# initialize an empty matrix to store pairwise distances 
pairwise_distances = np.zeros(shape =(len(counts), len(counts)))

# calculate euclidean distance between data points 
for cell in range(len(counts)): 
    for cell2 in range(len(counts)):
        distance = 0 
        for gene in range(len(counts[0])):
            distance += (counts[cell][gene] - counts[cell2][gene])**2
        distance = math.sqrt(distance)
        pairwise_distances[cell][cell2] = distance
