'''
Script to implement the bionj algorithm based on raw sequence distances.
I have not tested this extensively.
'''

import numpy as np

from sticcs import sticcs

## a distance matrix method that uses numerical arrays (i.e. alleles coded as 0, 1, 2,and 3, with negative numbers indicating missing data)

def numHamming(numArrayA, numArrayB):
    dif = numArrayA - numArrayB
    return np.mean(dif != 0)

def distMatrix(haplotype_array): #assumes each row is a different haplotype
    nanMaskTotal = haplotype_array>=0
    N,ln = haplotype_array.shape
    distMat = np.zeros((N,N))
    for i in range(N - 1):
        for j in range(i + 1, N):
            nanMask = nanMaskTotal[i,:] & nanMaskTotal[j,:]
            distMat[i,j] = distMat[j,i] = numHamming(haplotype_array[i,:][nanMask], haplotype_array[j,:][nanMask])
    return distMat


def bionj(distance_matrix):
    """
    Implements the BIONJ (Balanced Iterative Optimal Neighbor-Joining) algorithm,
    but returns a Tree object instead of a Newick string.
    
    Parameters:
    - distance_matrix: A 2D numpy array representing the pairwise distances between the species.
    
    Returns:
    - A sticc Tree object
    """
    n = len(distance_matrix)
    
    #labels for the indices of the distance matrix. These are kept separate because some will be removed anbd replaced
    labels = list(range(n))
        
    # Initialize tree object
    tree = sticcs.Tree(n_leaves=n)
    
    # Calculate the initial Q matrix (based on distance matrix)
    def calculate_Q_matrix(D):
        n = D.shape[0]
        Q = np.zeros_like(D)
        for i in range(n-1):
            for j in range(i + 1, n):
                Q[i,j] = Q[j,i] = (n-2) * D[i,j] - D[i,:].sum() - D[j,:].sum()
        #make the diagonal values infinite so that they are not ever minimal
        np.fill_diagonal(Q, np.inf)
        return Q
    
    # Build tree by iteratively merging nodes
    while n > 2:
        Q = calculate_Q_matrix(distance_matrix)
        
        # Find the pair of nodes with the minimum Q value
        min_Q_index = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = min_Q_index
        
        #Create the new distance row for the merged taxon (new node)
        new_distance_row = (distance_matrix[i, :] + distance_matrix[j, :] - distance_matrix[i, j]) / 2
        
        #Remove the columns corresponding to i and j from the new row
        new_distance_row = np.delete(new_distance_row, [i, j])
        
        #Create the new distance matrix by removing rows and columns for i and j
        new_distance_matrix = np.delete(distance_matrix, [i, j], axis=0)  # Remove rows for i and j
        new_distance_matrix = np.delete(new_distance_matrix, [i, j], axis=1)  # Remove columns for i and j
        
        #Add the new row (but without the self-distance 0)
        new_distance_matrix = np.vstack((new_distance_matrix, new_distance_row))  # Add new row
        
        #Add the new column with the self-distance 0 (representing the new node's distance to itself)
        new_distance_column = np.append(new_distance_row, 0)  # Add the self-distance (0) to the new column
        new_distance_matrix = np.hstack((new_distance_matrix, new_distance_column.reshape(-1, 1)))  # Add new column
        
        # Set the new distance matrix and update the number of nodes
        distance_matrix = new_distance_matrix
        
        # Add a new internal node
        new_node_id = tree.add_node(children=[labels[i],labels[j]])
        
        #update the labels
        labels = [labels[k] for k in range(n) if k not in [i,j]]
        labels.append(new_node_id)
        
        n -= 1
    
    return tree
