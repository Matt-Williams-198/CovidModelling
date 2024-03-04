# -*- coding: utf-8 -*-
###############################################################################
'''                Evaluation of the SIR-DDFT model in 2D                   '''
###############################################################################


__author__ = "Michael te Vrugt, Jens Bickmann, Raphael Wittkowski"
__copyright__ = "Copyright (C) 2020 Michael te Vrugt, Jens Bickmann, Raphael Wittkowski"
__license__ = "MIT"  
__version__ = "1.0"
__email__ = "raphael.wittkowski@uni-muenster.de"  


'''
    Imports
'''

import argparse
import os
from joblib import Parallel, delayed
import multiprocessing
from functools import partial
from timeit import default_timer as timer
import time
import joblib.parallel
import numpy as np
import matplotlib.pyplot as plt

# Import the solver
import solver as sol


###############################################################################


'''
    Constants and general settings
'''

# Parse commandline args
parser = argparse.ArgumentParser(
    description="Compute the SIR-DDFT model.",  
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "--notex", action="store_true",
    help="""Disable TeX in plotting""")
parser.add_argument(
    "--font", default="serif",
    help="""Set font type for plots""")
parser.add_argument(
    "-N", default=400,
    help="""Points in space domain""")
parser.add_argument(
    "-L", default=10.0,
    help="""Length of space domain""")
parser.add_argument(
    "--plot-interval", default=0.1,
    help="""Interval in which plots are generated""")
parser.add_argument(
    "--data-interval", default=1,
    help="""Interval in which plots are generated""")
parser.add_argument(
    "--nodatasave", action="store_true",
    help="""Disable data savings of every timestamp""")
args = parser.parse_args()


# Use LaTeX fonts in the plot
plt.rc('text', usetex = not args.notex)
plt.rc('font', family = args.font)

# Space discretization
N = args.N  # Do not go below N = 100

# Space domain
L = args.L  # Length of the space domain
dx = L/N
x_plot = np.linspace(0, L, N+1)
y_plot = x_plot # Quadratic domain
x_grid,y_grid = np.meshgrid(x_plot, y_plot) # Define grid

dt = 0.001 # Time step size
plot_interval = args.plot_interval  # Interval of data extraction
plot_interval2 = args.data_interval # Interval of data extraction

data_save = not args.nodatasave
eps = 1E-6 # A small value larger than zero

'''
    The field 'I_source' is the optional source term of I 
'''
global I_source
I_source_amplitude = 0.05 # Amplitude of the source
I_source = I_source_amplitude * np.exp(-1/(L * 2*(L/1000)) * ((L/2 - x_grid)**2 + (L/2 - y_grid)**2))  # Gaussian at (L/2,L/2)

'''                               
    The field 'bc_array' is the background array. It can be used for imposing Dirichlet boundary conditions (=0) in a separate dgl function.
    bc_array[i][j] = 0 corresponds to the fact that all relevant fields vanish at point {i,j}.
    bc_array[i][j] = 1 corresponds to no influence at point {i, j}.
'''
global bc_array
bc_array = np.ones_like(I_source)

                    
###############################################################################


def initializeState():
    # Initial conditions: a Gaussian for S and I, R=0 
    S = np.exp(-1/(L * 2*(L/50)) * ((L/2.0 - x_grid)**2+ (L/2.0 - y_grid)**2))  # Gaussian at (L/2,L/2)
    S = S / (dx**2 * np.sum(S) / L**2)*0.3543165399952919 # 0.3543165399952919 \approx \sqrt{\pi}/5, the \approx (instead of =) comes from the 1D density we mimic here
    I = 0.001 * S
    S = S - I
    R = np.zeros_like(S)  # No recovered persons at t=0
    return S, I, R


'''
    Solution for the convolution terms of the model
'''

# Pre-computation of the convolution integral
sig = 100 # Sets the parameter sigma of the model, cannot be changed later

saf = int(N // L) # Defines a safe zone at the boundary
matLen = len(x_plot) + 2 * saf # Length of the complete domain we want to convolve

# Initialize and compute the matrix used for the integration
matConv = np.zeros((matLen, matLen))
for i in range(matLen):
    for j in range(matLen):
        matConv[i, j] = np.exp(-sig * ((i - j)*dx)**2 ) # Equation given by model
        if j == 0 or j == matLen-1:
            matConv[i, j] *= 0.5 # Midpoints

def ConvPhi1dFF(phi):
    '''
        Calculate the convolution in 1D of a vector phi and return the convolved vector.
    '''
    rb = phi[:saf]
    lb = phi[-saf:]
    phi3 = np.concatenate((lb, phi, rb), axis = 0) # Extend the vector phi artificially with periodic boundaries 
    # Compute the convolution
    foo = np.zeros_like(phi3)
    foo = dx * np.dot(matConv, phi3) 
    return  foo[saf:N+saf+1] # Only return the middle part (from 0 to L)

def ConvPhi2dF(p):
    '''
        Calculate the convolution in 2D of a vector phi and return the convolved vector.
    '''
    # The different patches (analogy: numeric keypad with the center domain corresponding to the 5)
    q1 = p[-saf:, -saf:]
    q2 = p[:, -saf:]
    q3 = p[:saf, -saf:]
    q4 = p[-saf:, :]
    q6 = p[:saf, :]
    q7 = p[-saf:, :saf]
    q8 = p[:, :saf]
    q9 = p[:saf:, :saf]

    # Combine the different patches, first to lines (1-3, 4-6, 7-9) and then to one big domain (1-9)
    p9_123 = np.concatenate((q1,q2,q3), axis = 0)
    p9_456 = np.concatenate((q4,p,q6), axis = 0) 
    p9_789 = np.concatenate((q7,q8,q9))
    p9 = np.concatenate((p9_123,p9_456,p9_789), axis = 1)
    
    # Compute the convolution of that domain (1-9)
    foo = np.zeros_like(p9)
    foo = dx**2 * np.transpose(np.dot(
            matConv, np.transpose(np.dot(matConv, p9))
            ))
    bar = foo[saf:N + 1 + saf, saf:N + 1 + saf] # Extract the '5' from the domain (1-9)
    return bar # Only return the middle part


###############################################################################


'''
    Differential equations
'''

# Set the PERMANENT parameters (these cannot be varied by means of parameter scans)
Csigmasd = sig # DO NOT CHANGE THIS
Csigmasi = sig # DO NOT CHANGE THIS

def dgl_SIR_DDFT2D(state, para):
    '''
        Our SIR-DDFT model
        state: is the complete state of the system (S, I, R = state)
        para:  is a vector that contains the values for the parameters of the
        model that we want to investigate
    '''
    S, I, R = state # Extract the current state of the system
    
    # Extract the parameters of the model from para
    Di, Ccsd, Ccsi, CgammaI, Ds, Dr, Cc, Cw, Cm, CgammaS, CgammaR = para
    
    # SIR-DDFT model
    foo = sol.grad2d(Ccsd * ConvPhi2dF(S + 0 * R), dx) + sol.grad2d(Ccsi * ConvPhi2dF(I), dx) # Needed twice 
    dgl_S   =   Ds * sol.laplace2dr9(S, dx) - Cc * S * I - CgammaS * sol.div2d(np.array((S, S)) * foo, dx)
    dgl_I   =   Di * sol.laplace2dr9(I, dx) + Cc * S * I - (Cw + Cm) * I  - CgammaI * sol.div2d(np.array((I, I)) * Ccsi * sol.grad2d(ConvPhi2dF( S + I + 0 * R ), dx), dx)
    dgl_R   =  Cw * I + Dr * sol.laplace2dr9(R, dx)  - CgammaR * sol.div2d(np.array((R, R)) * foo, dx)
    
    # Return a 'state' vector of the RHS
    return np.array((dgl_S, dgl_I, dgl_R))

def dgl_SIR_DDFT2D_bc_w_source(state, para):
    '''
        Our SIR-DDFT model
        state: is the complete state of the system (S, I, R = state)
        para:  is a vector that contains the values for the parameters of the
        model that we want to investigate
    '''
    S, I, R = state # Extract the current state of the system
    
    # Apply the Dirichlet boundary conditions to the next time step
    S = bc_array * S
    I = bc_array * I
    R = bc_array * R
    
    # Extract the parameters of the model from para
    Di, Ccsd, Ccsi, CgammaI, Ds, Dr, Cc, Cw, Cm, CgammaS, CgammaR = para
    
    # SIR-DDFT model
    foo = sol.grad2d(Ccsd * ConvPhi2dF(S + R), dx) + sol.grad2d(Ccsi * ConvPhi2dF(I), dx) # Needed twice 
    dgl_S   =   Ds * sol.laplace2dr9(S, dx) - Cc * S * I - CgammaS * sol.div2d(np.array((S, S)) * foo, dx)
    dgl_I   =   I_source + Di * sol.laplace2dr9(I, dx) + Cc * S * I - (Cw + Cm) * I  - CgammaI * sol.div2d(np.array((I, I)) * Ccsi * sol.grad2d(ConvPhi2dF( S + R + I), dx), dx)
    dgl_R   =   Dr * sol.laplace2dr9(R, dx) + Cw * I - CgammaR * sol.div2d(np.array((R, R)) * foo, dx)
    
    # Apply the Dirichlet boundary conditions to the next time step
    dgl_S = bc_array * dgl_S
    dgl_I = bc_array * dgl_I
    dgl_R = bc_array * dgl_R
    
    # Return a 'state' vector of the RHS
    return np.array((dgl_S, dgl_I, dgl_R))


###############################################################################


'''
    Auxiliary functions
'''

# Construct the parameter vector
def makeParaVec(Di, Ccsd, Ccsi, CgammaI):
    # Diffusion constants
    Ds = 0.01
    Dr = Ds
    
    # SIR-model constants
    Cc = 1
    Cw = 0.1
    Cm = 0.0
    
    # Interaction constants
    CgammaS = 1 # 1
    CgammaR =CgammaS
    return np.array((Di, Ccsd, Ccsi, CgammaI, Ds, Dr, Cc, Cw, Cm, CgammaS, CgammaR))

def plot_flatten_the_curve_plots(name, time_array, sum_arrayS, sum_arrayI, sum_arrayR, ges_pop = 1):
    # Plot the flatten-the-curve plots
    figEnd = plt.figure(figsize=(3.4, 2.5))
    plt.plot(time_array, (sum_arrayS)/ges_pop, label=r'$\bar{S}$')
    plt.plot(time_array, sum_arrayI/ges_pop, label=r'$\bar{I}$')
    plt.plot(time_array, sum_arrayR/ges_pop, label=r'$\bar{R}$')
    plt.plot(time_array, (1 - (sum_arrayS + sum_arrayI + sum_arrayR) / ges_pop), label=r'$\bar{D}$')

    plt.legend(loc=5, frameon=False)
    plt.margins(0, 0)
    plt.subplots_adjust(top=0.96, right=0.96, left=0.1, bottom=0.17)
    plt.xlabel(r'$t / d$')
    # Save plot with an additional name
    add_name = 'flatten_the_curve_plot'  # Arbitrary
    plt.savefig(name + '/' + add_name + '.pdf', format='pdf', dpi=600)
    #plt.show()
    plt.close()
    return

def create_dir(name):
    # Create target directory if it does not exist
    if not os.path.exists(name):
        os.mkdir(name)
    else:
        return
    return

def make_circle_bc_array(radius):
    # Makes the global bc_array a circle with radius 'radius' (1 inside, 0 outside)
    for i in range(N+1):
        for j in range(N+1):
            x = x_plot[i]
            y = y_plot[j]
            if ((L/2 - x)**2 + (L/2 - y)**2) >= radius**2:
                bc_array[i, j] = 0
                        
def scan_input(input_val):
    '''
        Funtion used for the parallelization.
        Gets input values for Csi and Csi/Csd and returns relevant information regarding the specified system.
    '''
    Ccsi, rel, i, j = input_val # Extract data from input val
    
    # Perform the computation and save it
    name = 'parameterscan_Csi_{:1.0f}'.format(Ccsi)+ '_rel_{:1.1f}'.format(rel)
    
    # Define what parameters we want to consider
    para = makeParaVec(Di = 0.01, Ccsd = Ccsi/rel, Ccsi = Ccsi, CgammaI = 1)  # Define parameters to consider (Cc_1_sd, Cc_1_si, I1, Ca)
   
    # Calculate everything and get interesting data
    end_D, max_I, end_S, end_R, time_array, sum_arrayS, sum_arrayI, \
        sum_arrayR, ges_pop = solve_sweep(name, para, dt, t_max = 0) # 't = 0' makes the simulation run until I is very small

    # Return the data that should be stored in the grid
    return np.array((np.int(i), np.int(j), max_I, end_S, end_R))


###############################################################################


'''
    Time evolution
'''

def solve_sweep(name, parameters, dt_local, t_max = 0, dgl = dgl_SIR_DDFT2D):
    '''
        This function gets a DE and solves it from t=0 to t=t_max with an
        initial step size dt (an adaptive time stepping is used). The model is
        solved using the parameters in the parameter vector 'parameters'. The 
        solution and the data computed here are associated with the name 'name',
        which is also the directory where files are saved.
    '''
    
    # Create the directories for the data and images
    create_dir(name)
    create_dir(str(name + '//S'))
    create_dir(str(name + '//I'))
    create_dir(str(name + '//R'))
    create_dir(str(name + '//bf'))
    
    state = initializeState() # Initialize a state
    ges_pop = np.sum(state[0:2]) # Total population of the system (used for normalization later)

    # Arrays for the flatten-the-curve plots
    time_array = np.array(())
    sum_arrayS = np.array(())
    sum_arrayI = np.array(())
    sum_arrayR = np.array(())
    
    if t_max == 0: # Compute until I is small enough
        time = 0
        while True:
            current_time = time * plot_interval
            if current_time %  plot_interval2 == 0:
                # Plot the current state every time interval     
                
                plt.imshow(state[0], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$S(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/S/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()        
                plt.close()
                
                plt.imshow(state[1], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$I(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/I/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()
                plt.close()
                
                plt.imshow(state[2], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$R(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/R/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()
                plt.close()
                
                if data_save:
                    np.savetxt(name + '/bf/S' +str(int(time*plot_interval))+'.txt', state[0], delimiter='\t', fmt='%1.6f')
                    np.savetxt(name + '/bf/I' +str(int(time*plot_interval))+'.txt', state[1], delimiter='\t', fmt='%1.6f')
                    np.savetxt(name + '/bf/R' +str(int(time*plot_interval))+'.txt', state[2], delimiter='\t', fmt='%1.6f')
                
            # Insert the current data points for the flatten-the-curve plot
            time_array = np.insert(time_array, len(time_array), current_time)
            sum_arrayS = np.insert(sum_arrayS, len(sum_arrayS), np.sum(state[0]))
            sum_arrayI = np.insert(sum_arrayI, len(sum_arrayI), np.sum(state[1]))
            sum_arrayR = np.insert(sum_arrayR, len(sum_arrayR), np.sum(state[2]))
            
            # Calculate the state of the next time we want to take a picture 
            # (plot interval integration)
            state, dt_local  = sol.rkf(dgl, parameters, state, dt_local, plot_interval)
            time += 1
            if np.sum(state[1])/ges_pop < 1E-4 or current_time > 400: # Compute only if I is large enough or a large time is reached.
                break
                
    else:
        # Loop over the time
        for time in range(int(t_max/plot_interval)):
            if time*plot_interval %  plot_interval2 == 0:
                # Plot the current state every time interval     
                
                plt.imshow(state[0], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$S(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/S/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()        
                plt.close()
                
                plt.imshow(state[1], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$I(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/I/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()
                plt.close()
                
                plt.imshow(state[2], vmin = 0, vmax = 1, aspect = 'equal', extent = [0, L, 0, L])
                plt.colorbar()
                plt.title(r'$R(x, y, t = {:1.0f})$'.format(time*plot_interval))
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(name+'/R/{:03d}'.format(int(time*plot_interval))+'.png', format='png', dpi=600) # Save plot
                #plt.show()
                plt.close()
                
                if data_save:
                    np.savetxt(name + '/bf/S' +str(int(time*plot_interval))+'.txt', state[0], delimiter='\t', fmt='%1.6f')
                    np.savetxt(name + '/bf/I' +str(int(time*plot_interval))+'.txt', state[1], delimiter='\t', fmt='%1.6f')
                    np.savetxt(name + '/bf/R' +str(int(time*plot_interval))+'.txt', state[2], delimiter='\t', fmt='%1.6f')
                
            # Insert the current data points for the flatten-the-curve plot
            time_array = np.insert(time_array, len(time_array), time*plot_interval)
            sum_arrayS = np.insert(sum_arrayS, len(sum_arrayS), np.sum(state[0]))
            sum_arrayI = np.insert(sum_arrayI, len(sum_arrayI), np.sum(state[1]))
            sum_arrayR = np.insert(sum_arrayR, len(sum_arrayR), np.sum(state[2]))
            
            # Calculate the state of the next time we want to take a picture 
            # (plot interval integration)
            state, dt_local  = sol.rkf(dgl, parameters, state, dt_local, plot_interval)
    
    
    # Save the data of the flatten-the-curve plots as a txt file
    np.savetxt(name + '/' +name+'.txt', np.transpose(np.array(
            (time_array, sum_arrayS/ges_pop, sum_arrayI/ges_pop,
             sum_arrayR/ges_pop))), delimiter='\t', fmt='%1.6f')
    
    # Plot the data
    plot_flatten_the_curve_plots(name, time_array, sum_arrayS, sum_arrayI, sum_arrayR, ges_pop)
    
    # Get useful information like maximum number of infected individuals in this run, etc.
    max_I = np.max(sum_arrayI/ges_pop)
    end_S = sum_arrayS[-1]/ges_pop
    end_R = sum_arrayR[-1]/ges_pop
    end_D = (1 - (sum_arrayS + sum_arrayI + sum_arrayR) / ges_pop)[-1]
    
    # Return this useful data, etc.
    return end_D, max_I, end_S, end_R, time_array, sum_arrayS, sum_arrayI, \
        sum_arrayR, ges_pop


###############################################################################


'''
    Single-system runs: Example
    
    - Comparison of different values for C_si/C_sd with constant C_si
'''

# Perform a single-system simulation        
name1 = 'example_run_Csi_-20_rel_1' # Define a name for the simulation
para1 =makeParaVec(Di = 0.01, Ccsd = -20, Ccsi = -20, CgammaI = 1) # Define what parameters we want to consider
end_D, max_I, end_S, end_R, time_array1, sum_arrayS1, sum_arrayI1, sum_arrayR1, \
        ges_pop1 = solve_sweep(name1, para1, dt, t_max = 126) # Calculate everything and give some interesting facts for potential futher analysis
    
name2 = 'example_run_Csi_-20_rel_2' 
para2 =makeParaVec(Di = 0.01, Ccsd = -20/2.0, Ccsi = -20, CgammaI = 1)
end_D, max_I, end_S, end_R, time_array2, sum_arrayS2, sum_arrayI2, sum_arrayR2, \
        ges_pop2 = solve_sweep(name2, para2, dt, t_max = 126)

name3 = 'example_run_Csi_-20_rel_3' 
para3 =makeParaVec(Di = 0.01, Ccsd = -20/3.0, Ccsi = -20, CgammaI = 1)
end_D, max_I, end_S, end_R, time_array3, sum_arrayS3, sum_arrayI3, sum_arrayR3, \
        ges_pop2 = solve_sweep(name3, para3, dt, t_max = 126)


'''
    As a second example, we now can consider a circular domain and a source term.
    Note: To make the simulation run smoothly, one has to consider "2 * radius of circle + 1 <= L".
    This is a consequence of the efficient function np.roll used in solver.py. 
    For r=5 we therefore recommend L=11. For L=10 we recommend r<=4.5.
'''

make_circle_bc_array(4) # First, we make the background array to a disk with radius 4 

'''
    Reminder: 
        I_source_amplitude = 0.05
        I_source = I_source_amplitude * np.exp(-1/(L * 2*(L/1000)) * ((L/2 - x_grid)**2+ (L/2 - y_grid)**2)) 
    The desired source term and the source amplitude can be modified at the beginning of the file.
'''
name1 = 'example_source_term_circle' 
para1 =makeParaVec(Di = 0.01, Ccsd = -30/3.0, Ccsi = -30, CgammaI = 1)#makeParaVec(Di = 0.01, Ccsd = -10, Ccsi = -20, CgammaI = 1) # Define what parameters we want to consider
end_D, max_I, end_S, end_R, time_array1, sum_arrayS1, sum_arrayI1, sum_arrayR1, \
        ges_pop1 = solve_sweep(name1, para1, dt, t_max = 126, dgl=dgl_SIR_DDFT2D_bc_w_source) # Note the changed parameter 'dgl=dgl_SIR_DDFT2D_bc'.
                    

###############################################################################


'''
    Parameter scans: Example
    
    - Parameter scans in Csi and Csi/Csd in 2D
      (This can take up to several hours, depending on d_Csi and d_res.)
'''

# Define the grid we want to scan
d_Csi = 1 # Step size of the grid we scan in one direction
d_res = 0.1 # Step size of the grid we scan in the other direction
Csi_array = np.arange(-30, -0 + eps, d_Csi) # Array of the values we want to scan in one direction
res_array = np.arange(1.0, 3.0 + eps, d_res) # Array of the values we want to scan in the other direction

# Initialize arrays for the results obtained later
parameterscanS = np.zeros(shape = (len(Csi_array), len(res_array)))
parameterscanI = np.zeros(shape = (len(Csi_array), len(res_array)))

# Print which scan will be calculated
print ('\n')
print ('C_si parameters: ')
print (Csi_array)
print ('C_si/C_sd parameters: ')
print (res_array)

N_tot = len(Csi_array)*len(res_array) # Number of systems to compute
print ('\nTotal numbers of parameterscans: ' + str(N_tot))

# Loop over all configurations we want to calculate
i = 0 # These indicies i and j give just a more accessible indentification
j = 0
count = 0
input_val = np.array(())
for Ccsi in Csi_array:
    for rel in res_array:
        count += 1  # Number of systems we computed
        input_val = np.append(input_val, np.array((Ccsi, rel, i, j)))
        j += 1
    i += 1
    j = 0
input_val = np.reshape(input_val, (count, 4)) # This defines the workload

# Print combinations
print ('Parameter combinations that are going to be scanned: ')
print (input_val)


total_n_jobs = count

# patch joblib progress callback
class BatchCompletionCallBack(object):
    #  global variables
    global total_n_jobs
    global start_time
    start_time = timer() # Initialize a timer

    def __init__(self, dispatch_timestamp, batch_size, parallel):
        self.dispatch_timestamp = dispatch_timestamp
        self.batch_size = batch_size
        self.parallel = parallel

    def __call__(self, out):
        self.parallel.n_completed_tasks += self.batch_size
        this_batch_duration = time.time() - self.dispatch_timestamp

        self.parallel._backend.batch_completed(self.batch_size,
                                           this_batch_duration)
        self.parallel.print_progress()
        progress = self.parallel.n_completed_tasks / total_n_jobs # Get the progress
        
        # Make a print on the progress of the computation
        print(
            "\nProgress: [{0:50s}] {1:.1f}%".format('#' * int(progress * 50), progress*100) + "  time elapsed [sec]: {:1.1f}".format(timer()-start_time)
            , end="", flush=True)
        if self.parallel.n_completed_tasks == total_n_jobs:
            print('\n')

        if self.parallel._original_iterator is not None:
            self.parallel.dispatch_next()

joblib.parallel.BatchCompletionCallBack = BatchCompletionCallBack # Set the nicer looking callback

input = input_val # Arbitrary list
num_cores = multiprocessing.cpu_count() # Number of cores used = all available. This can be modified if one still has interest in using the machine for other tasks
print("Using {} cores.".format(num_cores))
workload_ = partial(scan_input)
# arg1 is fetched from input list
output = Parallel(n_jobs=num_cores)(delayed(workload_)(i) for i in input) # Compute the workload

# Extract the output data
k = 0
for i in range(len(Csi_array)):
    for j in range(len(res_array)):
        i, j, max_I, end_S, end_R = output[np.int(k)]
        parameterscanI[np.int(i), np.int(j)] = max_I
        parameterscanS[np.int(i), np.int(j)] = end_S
        k += 1

# Swap x and y axes for the plotting (we want res to correspond to the y axis)
parameterscanI = np.transpose(parameterscanI)
parameterscanS = np.transpose(parameterscanS)

# Save the data for later plotting
np.savetxt('parameterscan_example_Imax_2D.txt', parameterscanI, delimiter='\t', fmt='%1.6f')
np.savetxt('parameterscan_example_Sinfty_2D.txt', parameterscanS, delimiter='\t', fmt='%1.6f')

# For a code on how to plot these parameter scans, view 'SIR_DDFT_main_1D'.


###############################################################################
