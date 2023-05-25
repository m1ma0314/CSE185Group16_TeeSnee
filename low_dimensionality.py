low_dim = np.random.randn(len(counts), 2)

def calculate_low_dimension(low_dim):

    denominators = np.zeros(len(counts))
    numerators = np.zeros_like(counts)

    for i in range(len(counts)):
        for j in range(len(counts)):
            if(i == j): numerators[i][j] = 0
            else: 
                print((1+(low_dim[i] - low_dim[j])**2)**-1)
                numerators[i][j] = (1+(np.linalg.norm(low_dim[i] - low_dim[j])**2))**-1
                
        denominators[i]=np.sum(numerators[i,:])

    for i in range(len(counts)):
        for j in range(len(counts)):
            low_dim[i][j] = numerators[i][j]/denominators[i]
    return low_dim 
