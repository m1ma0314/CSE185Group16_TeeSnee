import numpy as np
import math

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
