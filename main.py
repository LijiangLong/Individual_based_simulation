import multiprocessing
from models import population
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pdb
from simulators import *


    
#neutral_simulation(0)
#for each simulator, plot 10 trajectories

#fig, axs = plt.subplots(3,1, figsize=(15, 6))
#axs = axs.ravel()
function_list = [neutral_simulation,positive_selection_simulation,heterozygous_advantage_simulation]
pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count()-1))

for i in range(3):
    figure = plt.figure(figsize=(5, 2))
    results = pool.map(function_list[i], range(10))
    for result in results:
        #axs[i].plot(np.arange(len(result)),result)
        plt.plot(np.arange(len(result)),result)
        
    #axs[i].set_xlim(0, 2000)
    #axs[i].set_ylim(0, 1)
    plt.xlim(0, 2000)
    plt.ylim(0, 1)

    plt.savefig(str(function_list[i]))
# with open('neural.pkl','w') as f:
#     pickle.dump([results], f)