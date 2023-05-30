import pandas as pd
import numpy as np
import math

#if num_iterations is not given by user, set to 1000
# add parameter to set the number of dimmensions, if not given by user, set to 2 

def calculate_tSNE(cellxgene_mat, output, num_iterations):
    '''
     Parameters
     ----------
     cellxgene_mat: 
     
     output:
     
     num_iterations: 
     
     num_dimensions: 
     
     References
     ----------
     [1] van der Maaten, L.J.P.; Hinton, G.E. Visualizing Data using t-SNE.
     Journal of Machine Learning Research 9:2579-2605, 2008.
     https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf
     
    '''
    
    #STEP 1: parse cell x gene matrix  
    if cellxgene_mat.endswith('.tsv'): #check if file is tab vs. comma separated file 
        count_matrix_df = pd.read_csv(cellxgene_mat, delimiter = "\t") 
    else: 
        count_matrix_df = pd.read_csv(cellxgene_mat) 
    
    #convert dataframe (cell by gene count matrix) to 2D numpy array 
    counts = count_matrix_df.to_numpy() 
    
    #STEP 2: calculate pairwise distances between each cell 
    pairwise_distances = calculate_pairwise_distances(counts)
    
    #STEP 3: calculate similarity matrix 
    target_perplexity = 40 #define target perplexity -- 40 is the value used in the paper section 4.2  
    similarities = find_similarities(pairwise_distances, target_perplexity)
    
    #STEP 4: calculate symmetrical conditional probabilities
    symmetrical_probabilities = symmetrical_probabilities(similarities)
    
    #STEP 5: calculate low dimensional counterpart to similarity matrix 
    rand_low_dim = np.random.normal(loc=0.0, scale=1.0, size=(len(counts), num_dimensions)) # generate random low dimensional space (2 dimensions) of Euclidean distances by drawing random samples from a normal (Gaussian) distribution
    
    #STEP 6: sample initial solution (gamma[0])
    gamma = np.zeros(num_iterations) #gamma is the low dimensional embedings 
    gamma[0] = np.zeros_like(low_dim)
    gamma[1] = low_dim
    
    learning_rate = 1000 #initially set learning rate to 1000, this will be updated after each iteration of optimization
    
    #most time-consuming part of function
    for t in range(1, num_iterations - 1): #iterate num_iterations -1 number of times for optimization
        
        if(t < 250): #set momentum rate, we used values as described in the paper section 3.4
            momentum = 0.5
        else:
            momentum = 0.8
        
        #STEP 7: calculate low dimensional counterpart to similarity matrix with new low dimensional embeddings
        low_dim_affinities = calculate_low_dimension(gamma[t])
        
        #STEP 8: calculate gradients between low-dimensional datapoints -- a function of pairwise Euclidean distances in the high-dimensional and low-dimensional space
        gradients = calculate_gradient(symmetrical_probabilities, low_dim_affinities, gamma[t])
        
        #STEP 9: calculate gamma[t]
        gamma[t+1] = gamma[t] + learning_rate * gradients + momentum*(gamma[t] - gamma[t-1]) #modified equation to substitute gamma[t-2] to gamma[t-1]     
    
    #STEP 10: perform k-means clustering to cluster data 
    max_num_clusters = 10
    clustered_data = calculate_clusters(gamma[len(gamma)-1], max_num_clusters)
    
    #STEP 11 (FINAL STEP !!!)
    plot(clustered_data, output)
    
