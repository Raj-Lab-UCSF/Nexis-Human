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

# Function to generate heat map of NEXIS output 
def heatmap(init_vec_method, Y):
    if init_vec_method == 'baseline':
        plt.imshow(Y, cmap='viridis', interpolation='none', aspect='auto')

    else: 
        # Exclude binary seeding location for binary initial vector so it does not drown out the signal in other regions (EDIT SEEDING LOCATION HERE)
        Y_modified = np.delete(Y, [14,48], axis=0) # NEED TO CHANGE for different seeding regions or different list of total regions
        plt.imshow(Y_modified, cmap='viridis', interpolation='none', aspect='auto')
        
    plt.colorbar()  # Add a color bar to map colors to values
    plt.title('Nexis Heatmap of Tau Time Series Across Regions')
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


# Cost function
def Nexis_error(params, patient_tau, stages, nexis_model):
    
    param1, param2, param3, param4, param5, param6, param7 = params 
    # param1 = alpha, param2 = beta, param3 = gamma, param4 = s, param5 = b, param6 = p, param7 = k

    # Parameters for simulate_nexis method
    parameters = [param1, param2, param3, param4, param5, param6, param7]  # [alpha, beta, gamma, s, b, p , k] 

    # Call the simulate_nexis method with the parameters
    Y = nexis_model.simulate_nexis(parameters)

    # For optimization, only take stages from Y that correspond to patient's stages 
    Y_edited = Y[:, stages]
    # Calculate R
    corr_coeff, p_value = pearsonr(patient_tau.flatten(), Y_edited.flatten())
    error = mse_matrix(patient_tau, Y_edited) + (1- corr_coeff)

    if np.isnan(error) or np.isinf(error):
            print(f"NaN or inf encountered with parameters: {params}")
            print(f"Result: {Y}")
        
    return error

   

    