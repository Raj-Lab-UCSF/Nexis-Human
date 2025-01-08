import os
import numpy as np
import scipy as sp
import scipy.io
from scipy.linalg import expm
import pandas as pandas
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from scipy.stats import pearsonr

# Function to generate heat map of NexIS output 
def heatmap(init_vec_method, seeding_indices,Y):
    if init_vec_method == 'baseline':
        plt.imshow(Y, cmap='viridis', interpolation='none', aspect='auto')

    else: 
        # Exclude binary seeding location for binary initial vector so it does not drown out the signal in other regions 
        Y_modified = np.delete(Y, [seeding_indices], axis=0)
        plt.imshow(Y_modified, cmap='viridis', interpolation='none', aspect='auto')
        
    plt.colorbar()  # Add a color bar to map colors to values
    plt.title('Nexis Heatmap of Diffusion Over Time')
    return plt


# Function to normalize by L2 norm
def normalize_by_l2_norm(matrix):
    l2_norms = np.linalg.norm(matrix, axis=1, keepdims=True)  # Calculate L2 norm for each row
    normalized_matrix = matrix / l2_norms  # Normalize each row by its L2 norm
    return normalized_matrix


# Function to calculate mean squared error
def mse_matrix(matrix1,matrix2):
    # Ensure the matrices have the same shape
    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same dimensions")
    return np.mean((matrix1 - matrix2) ** 2) 


# Cost function - once you've determined your cost function, add it here and import it to your running script

    