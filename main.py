from multiprocessing import Pool
from models import population
import matplotlib.pyplot as plt
import numpy as np
import pickle


def neutral_simulation(p):
    np.random.seed(p)
    population_simulator = population()
    return population_simulator.evolve_until_fix()
    
def sth(p):
    print(p)
pool = Pool(processes=4)
results = pool.map(neutral_simulation, range(4))


for result in results:
    plt.plot(np.arange(len(result)),result)

plt.show()

with open('neural.pkl','w') as f:
    pickle.dump([results], f)