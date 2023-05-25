import numpy as np 
def calculate_low_dimension(low_dim, counts):

    denominators = np.zeros(len(counts))
    numerators = np.zeros(shape=(len(counts),len(counts)))
    
    for i in range(len(counts)):
        for j in range(len(counts)):
            if(i ==j): numerators[i][j] = 0
            else: 
                numerators[i][j] = (1+(y[i] - y[j])**2)**-1
        denominators[i]=np.sum(numerators[i,:])

    for i in range(len(counts)):
        for j in range(len(counts)):
            low_dim[i][j] = numerators[i][j]/denominators[i]
    return low_dim 
