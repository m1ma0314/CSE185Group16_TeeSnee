# CSE185Project

Data (count matrix) taken from: 
https://singlecell.broadinstitute.org/single_cell/study/SCP1526/functional-metabolic-and-transcriptional-maturation-of-human-pancreatic-islets-derived-from-stem-cells#study-download
* need to make a Broad Institute account (Free!) to download the data


After downloading the data: 

```
gunzip endocrine.counts.mtx.gz
```

Functions Draft: 
1. Open & organize data into dataframe:
  * Read in data as a pandas dataframe 
  * Convert data into high dimensional space 
  * Use numpy 2D array
2. Calculate pairwise distances between all data points:
  * for each vector (representing the cell
