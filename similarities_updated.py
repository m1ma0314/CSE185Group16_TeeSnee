import numpy as np
import math

def find_similarities(distance_matrix, target_perplexity):
#sigma will produce the perplexity close to the user defined perplexity
#try a sigma each time and calculate the corresponding perplexity
#if calculated perplexity is larger; set sigma_max to the sigma, otherwise,
#set sigma_min to the sigma

  count = distance_matrix.shape[0]
  print(distance_matrix)
  print(distance_matrix.shape)
  denominators = np.zeros(count)
  nominators = np.zeros(shape=(count,count))
  
  p_ij = np.zeros(shape=(count,count))
  sigma = 1.0
  # iterate through each value and calculate the probabilities
  #norm = np.linalg.norm(distance_matrix)
  for i in range(count):
    for j in range(count):
      nominators[i][j] = np.exp((-distance_matrix[i][j]**2)/(2*sigma**2))
    denominators[i]=np.sum(nominators[i,:])
  
  for i in range(count):
    for j in range(count):
      p_ij[i][j] = nominators[i][j]/denominators[i]
  
  # calculate perplexity corresponding to p_i we have
  entropy = np.zeros(count)
  perp = np.zeros(count)
  for i in range(count):
    entropy[i] = -np.sum(p_ij[i,:]*np.log2(p_ij[i,:]))
    perp[i] = np.power(2,entropy[i])
  
  #calculate sigmas for each row
  sigmas = np.ones(count)
  for i in range(count):
    #initial sigma value 
    sigma2 = 1.0
    sigma_min = -np.inf
    sigma_max = np.inf

    while True:
      if math.isnan(perp[i]):
        break
    
      if perp[i]<target_perplexity:
        sigma_max = sigma2
      else:
        sigma_min = sigma2

      if np.abs(perp[i]-target_perplexity) < 1e-5:
        break

      sigma2 = np.average(sigma_max,sigma_min)
    
    sigmas[i] = sigma2
  
  #after finding the real sigmas, calculate the real similarities matrix
  similarities = np.zeros(shape=(count,count))
  for i in range(count):
    for j in range(count):
      nominators[i][j] = np.exp((-distance_matrix[i][j]**2)/(2*sigmas[i]**2))
    denominators[i]=np.sum(nominators[i,:])

  for i in range(count):
    for j in range(count):
      similarities[i][j] = nominators[i][j]/denominators[i]

  return similarities

find_similarities(pairwise_distances,10)

