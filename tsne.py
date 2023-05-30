import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
from sklearn.cluster import KMeans
import warnings
import os

# Suppress all warnings
warnings.filterwarnings("ignore")

def calculate_pairwise_distances(counts):
    # initialize an empty matrix to store pairwise distances 
    pairwise_distances = np.zeros(shape =(len(counts), len(counts)))

    # calculate euclidean distance between data points 
    for cell in range(len(counts)): 
        for cell2 in range(len(counts)):
            distances = math.sqrt(np.sum((counts[cell] - counts[cell2])**2))
            pairwise_distances[cell][cell2] = distances
    return pairwise_distances

def find_similarities(distance_matrix, target_perplexity):
#sigma will produce the perplexity close to the user defined perplexity
#try a sigma each time and calculate the corresponding perplexity
#if calculated perplexity is larger; set sigma_max to the sigma, otherwise,
#set sigma_min to the sigma

  count = distance_matrix.shape[0]
  denominators = np.zeros(count)
  numerators = np.zeros(shape=(count,count))

  p_ij = np.zeros(shape=(count,count))
  sigmas = np.zeros(count)

  #array storing original sigmas values
  sigmas = np.std(distance_matrix,axis=1)

  # iterate through each value and calculate the probabilities
  iterations = 10

  for row in range(count):  
    sigma_min = 0
    sigma_max = np.inf
    sigma = sigmas[row]

    for sigma_search in range(iterations):
      for i in range(count):
        for j in range(count):
          numerators[i][j] = np.exp((-np.linalg.norm(distance_matrix[i]-distance_matrix[j])**2)/(2*sigmas[row]**2))
        denominators[i]=np.sum(numerators[i,:])
  
      # calculate the conditional probabilities
      np.fill_diagonal(p_ij, 0)
      for j in range(count):
          p_ij[i][j] = numerators[i][j]/denominators[i]


      # calculate entropy and perplexity corresponding to p_i we have
      entropy = 0
      ε = np.nextafter(0,1)
  
      for i in range(count):
        for j in range(count):
          p_new = np.maximum(p_ij[i][j],ε)
          entropy += (p_new*np.log2(p_new))
      curr_perp = np.power(2,-entropy)

      if math.isnan(curr_perp) or np.abs(curr_perp-target_perplexity) < 1e-5:
        break

      # binary search for sigma
      if curr_perp < target_perplexity:
        sigma_min = sigma
      else:
        sigma_max = sigma

      sigma = (sigma_max+sigma_min)/2
      
    sigmas[row] = sigma
  
  #after finding the real sigmas, calculate the real similarities matrix
  similarities = np.zeros(shape=(count,count))
  for i in range(count):
    for j in range(count):
      numerators[i][j] = np.exp((-np.linalg.norm(distance_matrix[i]-distance_matrix[j])**2)/(2*sigmas[i]**2))
    denominators[i]=np.sum(numerators[i,:])
  
  for i in range(count):
    for j in range(count):
      similarities[i][j] = numerators[i][j]/denominators[i]

  return similarities

def symmetrical_probabilities(similarities):
    symmetrical_prob = np.zeros_like(similarities)
    n = len(similarities)**2
    for i in range(len(similarities)):
        for j in range (len(similarities[i])):
            symmetrical_prob[i][j] = (similarities[i][j] + similarities[j][i])/2*n
    return symmetrical_prob 

def calculate_low_dimension(low_dim):
    ld_pairwise_sim = np.zeros(shape =(len(low_dim), len(low_dim)))
    denominators = np.zeros(len(low_dim))
    numerators = np.zeros(shape =(len(low_dim), len(low_dim)))
    
    for i in range(0,len(low_dim)):
        for j in range(len(low_dim)):
            if(i == j): numerators[i][j] = 0
            else: 
                numerators[i][j] = ((1+(np.linalg.norm(low_dim[i] - low_dim[j])**2))**-1)    
        denominators[i]=np.sum(numerators[i,:])
        
    ld_pairwise_sim= numerators/np.sum(denominators)
    ε = np.nextafter(0,1)
    ld_pairwise_sim = np.maximum(ld_pairwise_sim, ε)
    
    return ld_pairwise_sim 
def calculate_gradient(symmetrical_probabilities, ld_sim, low_dim):
    gradients = np.zeros(shape = (len(low_dim),len(low_dim[0])))
    sym_minus_ld_sim = np.zeros(shape = (len(low_dim),len(low_dim)))
    for i in range(len(low_dim)):
        for j in range(len(low_dim)):
            sym_minus_ld_sim[i][j] = symmetrical_probabilities[i][j] - ld_sim[i][j]
            gradients[i] += 4 * sym_minus_ld_sim[i][j] * (low_dim[i] - low_dim[j]) * (1 + np.linalg.norm(low_dim[i]-low_dim[j])**2)**-1
            
    return gradients

def calculate_clusters(low_dim_embeddings, max_num_clusters):
    inertia = []
    # code based on example https://www.geeksforgeeks.org/elbow-method-for-optimal-value-of-k-in-kmeans/#
    
    n_clusters = range(1, max_num_clusters)
    for i in range(1,max_num_clusters):
        kmeans = KMeans(i)
        kmeans.fit(low_dim_embeddings)
        inertia.append(kmeans.inertia_)
    
    #find elbow point (maximum second derivative val
    _2nd_derivatives = np.diff(np.diff(inertia))
    elbow_index = np.argmax(_2nd_derivatives) + 1
    elbow_point = n_clusters[elbow_index]

    #calculate kmeans based on final number of clusters 
    kmeans = KMeans(elbow_point)
    kmeans.fit(low_dim_embeddings)
    labels = kmeans.labels_
    return labels

def plot(low_dim_embeddings, clustered_data, output):
    #check if output is a directory or filename 
    xcoords = [vals[0] for vals in low_dim_embeddings]
    ycoords = [vals[1] for vals in low_dim_embeddings]
    
    plt.scatter(xcoords, ycoords, c= clustered_data)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('tSNE 1')
    plt.ylabel('tSNE 2')
    if os.path.isdir(output):  
       plt.savefig(output + "tsneplot.png") 
    else: plt.savefig(output) 
    
#if num_iterations is not given by user, set to 1000
# add parameter to set the number of dimmensions, if not given by user, set to 2 

def calculate_tSNE(cellxgene_mat, output, target_perplexity):
    '''
     Parameters
     ----------
     cellxgene_mat: 
     
     output:
     
     References
     ----------
     [1] van der Maaten, L.J.P.; Hinton, G.E. Visualizing Data using t-SNE.
     Journal of Machine Learning Research 9:2579-2605, 2008.
     https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf
     
    '''
    num_iterations = 100
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
    # target_perplexity = 40 #define target perplexity -- 40 is the value used in the paper section 4.2  
    similarities = find_similarities(pairwise_distances, target_perplexity)
    
    #STEP 4: calculate symmetrical conditional probabilities
    sym_probabilities = symmetrical_probabilities(similarities)
    
    #STEP 5: calculate low dimensional counterpart to similarity matrix 
    rand_low_dim = np.random.normal(loc=0.0, scale=1.0, size=(len(counts), 2)) # generate random low dimensional space (2 dimensions) of Euclidean distances by drawing random samples from a normal (Gaussian) distribution
    
    #STEP 6: sample initial solution (gamma[0])
    gamma = [0] * num_iterations #gamma is the low dimensional embedings 
    gamma[0] = np.zeros_like(rand_low_dim)
    gamma[1] = rand_low_dim
    
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
        gradients = calculate_gradient(sym_probabilities, low_dim_affinities, gamma[t])
        
        #STEP 9: calculate gamma[t]
        gamma[t+1] = gamma[t] + learning_rate * gradients + momentum*(gamma[t] - gamma[t-1]) #modified equation to substitute gamma[t-2] to gamma[t-1]     
    
    #STEP 10: perform k-means clustering to cluster data 
    max_num_clusters = 10
    clustered_data = calculate_clusters(gamma[len(gamma)-1], max_num_clusters)
    #STEP 11 (FINAL STEP !!!)
    plot(gamma[len(gamma)-1], clustered_data, output)
    
"""
Command-line script to run t-SNE
"""

def main():
  parser = argparse.ArgumentParser(
      prog="tsne.py",
      description="command-line script to generate tsne plots"
  )

  # Input
  parser.add_argument("filename", help="gene data file (specify if zipped or not)",type=str)
  parser.add_argument("-o", "--output", help="output directory",type=str)
  parser.add_argument("-p","--target_perplexity", help="user specificed perplexity",type=int,metavar="PERPLEXITY",required=False)
  parser.add_argument("-z","--zipped",help="unzip file if input is zipped",action="store_true")

  args = parser.parse_args()

  if not args.target_perplexity:
    args.target_perplexity = 40
  if args.zipped:
    print("File is zipped. Extracting and reading as CSV:", args.filename)
    unzipped_filename = args.filename
    #check if zipped file is a csv
    if unzipped_filename.endswith('.tsv'): #check if file is tab vs. comma separated file 
        unzipped_filename = pd.read_csv(unzipped_filename,compression='gzip', delimiter = "\t")
    else: 
        unzipped_filename = pd.read_csv(unzipped_filename,compression='gzip')
        calculate_tSNE(unzipped_filename,args.output, args.target_perplexity)
  else:
    calculate_tSNE(args.filename,args.output, args.target_perplexity)
      
if __name__ == "__main__":
  main()
    
