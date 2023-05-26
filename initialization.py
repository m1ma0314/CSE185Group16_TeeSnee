def initialization(distance_matrix, n_dimensions=2):
  return np.random.normal(loc=0, scale=1e-4, size=(len(distance_matrix),n_dimensions))
