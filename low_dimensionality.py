import numpy as np 

low_dim = np.random.randn(len(counts), 2)

def calculate_low_dimension(counts, low_dim):
    ld_pairwise_sim = np.zeros(shape =(len(counts), len(counts)))
    denominators = np.zeros(len(counts))
    numerators = np.zeros(shape =(len(counts), len(counts)))
    
    for i in range(len(low_dim)):
        for j in range(len(low_dim[0])):
            print(i, j)
            if(i == j): numerators[i][j] = 0
            else: 
                numerators[i][j] = (1+(np.linalg.norm(low_dim[i] - low_dim[j])**2))**-1
                
        denominators[i]=np.sum(numerators[i,:])
    print(numerators.shape)  
    
    for i in range(len(counts)):
        for j in range(len(counts)):
            ld_pairwise_sim[i][j] = numerators[i][j]/denominators[i]
            
    return ld_pairwise_sim 
