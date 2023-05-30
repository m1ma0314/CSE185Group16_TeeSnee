import matplotlib.pyplot as plt
def plot(clustered_data, output):
    #check if output is a directory or filename 
    xcoords = [vals[0] for vals in clustered_data]
    ycoords = [vals[1] for vals in clustered_data]
    clusters = [vals[2] for vals in clustered_data]
    
    plt.scatter(xcoords, ycoords, c= clusters)
