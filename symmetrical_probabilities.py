import numpy as np 
def symmetrical_probabilities(similarities):
    symmetrical_prob = np.zeros_like(similarities)
    n = len(similarities)**2
    for i in range(len(similarities)):
        for j in range (len(similarities[i])):
            symmetrical_prob[i][j] = (similarities[i][j] + similarities[j][i])/2*n
    return symmetrical_prob 
