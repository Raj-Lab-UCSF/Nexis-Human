
import os
import numpy as np
import scipy as sp
import scipy.io
from scipy.linalg import expm
import pandas as pandas
import matplotlib.pyplot as plt
from scipy.integrate import odeint

class run_Nexis:
    def __init__(self,C_,U_,init_vec_,t_vec_,w_dir_=0,use_baseline_=0, logistic_term_=0):
        self.C = C_ # Connectivity matrix, nROI x nROI
        self.U = U_ # Matrix or vector of cell type or gene expression, nROI x nTypes
        self.init_vec = init_vec_ # Binary vector indicating seed location OR array of baseline pathology values, nROI x 1
        self.t_vec = t_vec_ # Vector of time points to output model predictions, 1 x nt
        self.w_dir = w_dir_ # Binary flag indicating whether to use directionality or not 
        self.use_baseline = use_baseline_ # Binary flag indicating whether you are using baseline or a binary seed to initialize the model
        self.logistic_term = logistic_term_ 
        
    def simulate_nexis(self, parameters):
        """
        Returns a matrix, Y, that is nROI x nt representing the modeled Nexis pathology
        given the provided parameters. alpha, beta, and gamma should be nonnegative scalars;
        s should be bounded between 0 and 1; b and p can be negative or positive
        """
        # Define parameters
        ntypes = np.size(self.U,axis=1)
        alpha = parameters[0] # global growth/aggregation rate 
        beta = parameters[1] # global diffusivity/spread rate 
        if self.use_baseline:
            gamma = 1
        else:
            gamma = parameters[2] # seed rescale value 
        if self.w_dir==0:
            s = 0.5
        else:
            s = parameters[3] # directionality (0 = anterograde, 1 = retrograde)
        b = np.transpose(parameters[4:(ntypes+4)]) # mediator-dependent spread/diffusion modifier 
        p = np.transpose(parameters[(ntypes+4):6]) # mediator-dependent growth/aggregation modifier 
        k = parameters[6] # Carrying capacity ADDED
        
        # Define starting pathology x0
        x0 = gamma * self.init_vec
        
        # Define diagonal matrix Gamma containing spread-independent mediator terms
        s_p = np.dot(self.U,p)
        Gamma = np.diag(s_p) + (alpha * np.eye(len(s_p))) 

        # Define Laplacian matrix L (which contains connectivity matrix and mediator terms)
        C_dir = (1-s) * np.transpose(self.C) + s * self.C
        coldegree = np.sum(C_dir,axis=0)
        L_raw = np.diag(coldegree) - C_dir
        s_b = np.dot(self.U,b)
        s_b = np.reshape(s_b,[len(s_b),1])
        S_b = np.tile(s_b,len(s_b)) + np.ones([len(s_b),len(s_b)])
        L = np.multiply(L_raw,np.transpose(S_b))

        # Define system dydt = Ax
        A = Gamma - (beta * L)

        # Solve 
        if self.logistic_term:
            y = self.sim_logistic(self.t_vec,x0,A,Gamma,k) 
        else:
            y = self.forward_sim(A,self.t_vec,x0)

        return y
    
    # Solve via analytic method (no logistic term)
    def forward_sim(self,A_,t_,x0_):
        y_ = np.zeros([np.shape(A_)[0],len(t_)])
        for i in list(range(len(t_))):
            ti = t_[i]
            y_[:,i] = np.dot(expm(A_*ti),np.squeeze(x0_)) 
        return y_
    
    # Solve via odeint with logistic term
    def sim_logistic(self,t_,x0_,A_,Gamma_,k_):

        # Define ODE function with a logistic term
        def ode_func(y, t, A, Gamma, k):
            dydt = np.dot(A, y) - np.dot(Gamma,np.square(y)) / k
            return dydt

        # Initial condition
        y0 = x0_

        # Solve ODE using odeint
        sol = odeint(ode_func, y0, t_, args=(A_,Gamma_,k_))

        # Transpose so that sol is an array with dim nROI x time points
        sol = sol.T

        return sol