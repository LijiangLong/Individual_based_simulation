from models import population
import numpy as np



def neutral_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1,1])
    return population_simulator.evolve_until_fix()
    
    
def positive_selection_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1,1],max_generation=50,outcrossing_rate = 0.01,incompatible = 1)
    return population_simulator.evolve_until_fix()
    
    
def heterozygous_advantage_simulation(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1.02,1])
    return population_simulator.evolve_until_fix()
    

def Neutral_Evolution_with_incompatible_gene_outcrossing1(p):
    np.random.seed(p)
    population_simulator = population(population_size = 10000,growth_speed_array = [1,1,1],incompatible = 1)
    return population_simulator.evolve_until_fix()
    
def Neutral_Evolution_with_incompatible_gene_outcrossing001(p):
    np.random.seed(p)
    population_simulator = population(growth_speed_array = [1,1,1],incompatible = 1,outcrossing_rate = 0.01,max_generation=2000)
    return population_simulator.evolve_until_fix()

def Neutral_Evolution_with_incompatible_gene_outcrossing001_large_population(p):
    np.random.seed(p)
    population_simulator = population(population_size = 10000,growth_speed_array = [1,1,1],incompatible = 1,outcrossing_rate = 0.01)
    return population_simulator.evolve_until_fix()

def Neutral_Evolution_with_incompatible_gene_outcrossing001_large_population_fitness_balancing(p):
    np.random.seed(p)
    population_simulator = population(population_size = 10000,growth_speed_array = [1,1.0035,1.007],incompatible = 1,outcrossing_rate = 0.01)
    return population_simulator.evolve_until_fix()
    
def Neutral_Evolution_with_incompatible_gene_outcrossing0001_large_population_fitness_balancing(p):
    np.random.seed(p)
    population_simulator = population(population_size = 10000,growth_speed_array = [1,1.0035,1.007],incompatible = 1,outcrossing_rate = 0.001)
    return population_simulator.evolve_until_fix()
    
