import numpy as np

def find_sigma(diff_i, i, target_perplexity):
  result = np.inf

  norm = np.linalg.norm(diff_i,axis=1)
  std_norm = np.std(norm)

  #I am not sure what values to put so I just used the same as the paper
  for s in np.linspace(0.01*std_norm,5*std_norm,150):
    p = np.exp(-norm**2/(2*s**2))
    p[i] = 0

    epsilon = np.nextafter(0,1)
    p_new = np.maximum(p/np.sum(p),epsilon)

    # calculate Shanoon entropy
    H = np.sum(-p_new*np.log2(p_new))

    #find sigma that produces P with a fixed perplexity that is specified by the user
    if np.abs(np.log(target_perplexity)-H*np.log(2)) < np.abs(result):
      result = np.log(target_perplexity)-H*np.log(2)
      sigma = s

  return sigma


def similarities_matrix(distance_matrix, target_perplexity):
  # set to this value for ease of similarities calculation
  sigma = 1.0

  size = distance_matrix.shape[0]
  similarity = np.zeros(shape=(size,size))
  
  

  for i in range(size):
    diff = distance_matrix[i] - distance_matrix

    sigma = find_sigma(diff, i, target_perplexity)
    norm = np.linalg.norm(diff,axis=1)
    
    
    for j in range(size):
      if i==j:
        similarity[i][j]==0
      else:similarity[i][j]= np.exp(-norm[j]**2/(2*sigma**2))
    
    for j in range(size):
      if i==j:
        similarity[i][j]==0
      else:similarity[i][j] = np.exp(-norm[j]**2/(2*sigma**2))/np.sum(similarity[i,:])
      
  
  return similarity