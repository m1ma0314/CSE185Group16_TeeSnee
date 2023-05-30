# CSE185Project t-SNE implementation
This is group 16 final project for CSE185. It follows the algorithm description in  the [original t-SNE parper](https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf) to implement t-SNE tool.

# Table of contents <a name="toc"></a >
- [Sample Data Download](#data)
- [Installation Instructions](#install)
- [Basic Usage](#usage)
- [Complete Usage](#instruction)
- [Contributors](#credit)


# Sample Data Download <a name="data"></a>
[BACK TO TABLE OF CONTENTS](#toc)

Sample data (count matrix) taken from: 
[https://singlecell.broadinstitute.org/single_cell/study/SCP1526/functional-metabolic-and-transcriptional-maturation-of-human-pancreatic-islets-derived-from-stem-cells#study-download](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164073)
* need to make a Broad Institute account (Free!) to download the data

After downloading the data: 
```
gunzip endocrine.counts.mtx.gz
```

# Installation Instruction <a name="install"></a>
[BACK TO TABLE OF CONTENTS](#toc)


Installation requires the `numpy`, `pandas`, `matplotlib`, `scikit-learn` libraries to be installed. You can install these with `pip`:

```
pip install -r requirements.txt
```
Once required libraries are installed, you can install `tsne` with the following command:
```

```

Functions Draft:

Algorithm as described in original t-sne paper https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf

1. Open & organize data into dataframe:
  * Read in data as a pandas dataframe 
  * Convert data into high dimensional space 
  * Use numpy 2D array
  
2. Calculate pairwise distances between all data points:
  * for each vector (representing the cell), look at each index (representing the gene counts) and subtract to a second cell at the same index
  * calculate distance between each cell, euclidean distance = sqrt(sum(subtracted values for each index))
  * repeat for all cells 
  * store distances in 2D numpy array 
  
3. Calculate similarity matrix
  *  For each cell, subtract all other cell x cell 1D matrices from current cell’s 1D matrix
  *  Get sum of values calculated in previous step
  *  *how to calculate sigma???*
  *  sigma is determined by binary search using a perplexity value 
  *  the perplexity value is related to Shannon Entropy value
  *  normalize similarity matrix 
<img width="350" alt="Screenshot 2023-05-18 at 1 29 39 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/997a4eea-2650-4d4b-b39d-e2b5193a27a3">

4. Calculate low dimensional counterpart to high dimensional similarity matrix 
<img width="350" alt="Screenshot 2023-05-18 at 1 30 28 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/7b6c75da-2372-4cd9-b092-86b08590d5a6">


  * *what is y?* 

5. Calculate symmetrical probabilities
<img width="100" alt="Screenshot 2023-05-18 at 1 32 37 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/0883fc13-748b-46c3-aefd-fd969cc2fd74">

  * Check that sum of probabilities of i from j is greater than 1/2n
  
6. Calculate gradient function of symmetrical probabilities 
<img width="350" alt="Screenshot 2023-05-18 at 1 34 39 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/a500f4cf-aac3-466a-bda7-8150a9477d2b">

7. Sample initial solution to find gamma 
<img width="100" alt="Screenshot 2023-05-18 at 1 35 22 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/5475e934-fd2f-4303-b330-918ca6023875">

  * gamma = 
  
8. Find gamma(t) 
<img width="200" alt="Screenshot 2023-05-18 at 1 36 47 PM" src="https://github.com/m1ma0314/CSE185Project/assets/59674595/f245809c-2f02-4c2d-ba6d-d1afa489ad5d">

  * Once gamma(t) converges (no longer changes after x number of iterations) we have identified the low dimensional map 
   
9. Plot 
