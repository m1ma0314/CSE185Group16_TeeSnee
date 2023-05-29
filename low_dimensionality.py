import numpy as np 

low_dim = np.random.randn(len(counts), 2)

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
