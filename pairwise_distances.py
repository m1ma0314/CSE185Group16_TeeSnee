def calculate_pairwise_distances(counts):
    # initialize an empty matrix to store pairwise distances 
    pairwise_distances = np.zeros(shape =(len(counts), len(counts)))

    # calculate euclidean distance between data points 
    for cell in range(len(counts)): 
        for cell2 in range(len(counts)):
            distances = math.sqrt(np.sum((counts[cell] - counts[cell2])**2))
            pairwise_distances[cell][cell2] = distances
    return pairwise_distances
