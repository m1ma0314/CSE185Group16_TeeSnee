import numpy as np
def calculate_gradient(symmetrical_probabilities, ld_sim, low_dim):
    gradients = np.zeros(shape = (len(low_dim),len(low_dim[0])))
    sym_minus_ld_sim = np.zeros(shape = (len(low_dim),len(low_dim)))
    for i in range(len(low_dim)):
        for j in range(len(low_dim)):
            sym_minus_ld_sim[i][j] = symmetrical_probabilities[i][j] - ld_sim[i][j]
            gradients[i] += 4 * sym_minus_ld_sim[i][j] * (low_dim[i] - low_dim[j]) * (1 + np.linalg.norm(low_dim[i]-low_dim[j])**2)**-1
            
    return gradients
