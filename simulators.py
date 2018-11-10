from models import population
import numpy as np



def neutral_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1,1])
    return population_simulator.evolve_until_fix()
    
    
def positive_selection_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1.01,1.02])
    return population_simulator.evolve_until_fix()
    
    
def heterozygous_advantage_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1.02,1])
    return population_simulator.evolve_until_fix()