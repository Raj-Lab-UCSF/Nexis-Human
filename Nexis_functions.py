import os
import numpy as np
import scipy as sp
import scipy.io
from scipy.linalg import expm
import pandas as pandas
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize

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


# Function to plot total tau over time
def total_tau_plot(data1, data2, name1, name2):

    total_tau1 = np.sum(data1, axis=0)
    total_tau2 = np.sum(data2, axis=0)

    time_points = np.linspace(0, 99, 100)

    plt.figure(figsize=(10, 6))

    # Plot total_tau_Y
    plt.plot(time_points, total_tau1, marker='o', linestyle='-', color='b', label= name1)  

    # Plot total_tau_EBM on the same graph
    plt.plot(time_points, total_tau2, marker='x', linestyle='--', color='r', label= name2)

    plt.title('Total tau over time') 
    plt.xlabel('Time')  # Label the x-axis
    plt.ylabel('Total tau')  # Label the y-axis
    plt.xticks(rotation=45)  # Optional: rotate x-axis labels for better readability
    plt.legend()  # Add a legend to distinguish the two lines

    # Show the plot
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