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
#function_list = [neutral_simulation,positive_selection_simulation,heterozygous_advantage_simulation]
function_list = [Neutral_Evolution_with_incompatible_gene_outcrossing0001_large_population_fitness_balancing]
pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count()-1))
for i in range(len(function_list)):
    figure = plt.figure(figsize=(5, 2))
    results = pool.map(function_list[i], range(10))
    for result in results:
        plt.plot(np.arange(len(result)),result)
    
    plt.xlim(0, 2000)
    plt.ylim(0, 1)

    plt.savefig(str(function_list[i]))
