# write a model to simulate evolution 
import numpy as np
import pdb




class individual:
    def __init__(self, genotype):
        self.genotype = genotype
    


class population:
    def __init__(self, major_allele_freq = 0.5,population_size = 1000,growth_speed_array = [1,1,1],outcrossing_rate = 1,incompatible = 0):
        # create a numpy 2d array, each row is an individual
        # each individual has a allele 0 or 1
        self.population_size = population_size
        self.growth_speed = np.array(growth_speed_array)
        self.outcrossing_rate = outcrossing_rate
        self.incompatible = incompatible
        
        genotypes = np.ndarray((population_size,2))
        
        # 0 or 1 is assigned to the genotypes according the major allele frequency
        for i in range(population_size):
            for j in range(2):
                random_allele = np.random.choice(2,1,p=[major_allele_freq,1-major_allele_freq])
                genotypes[i,j] = random_allele
        self.genotypes = genotypes
        
        

    def calculate_allele_frequency(self):
        minor_allele_frequency = np.sum(self.genotypes)/(2*self.population_size)
        return 1-minor_allele_frequency
        
    def one_generation(self):
        #based on genotypes and growth speed, determine the probability of choosing a particular parent
        
        parents_prob = self.growth_speed[np.sum(self.genotypes,axis=1).astype(int)]
        parents_prob = parents_prob/np.sum(parents_prob)
        #generate offspring to the population size
        new_genotypes = np.ndarray((self.population_size,2))
        for i in range(self.population_size):
            while True:
            #determine if it's from a selfing parent or crossing progeny
                choices = ['selfing','mating']
                choice = np.random.choice(choices,1,p=[1-self.outcrossing_rate,self.outcrossing_rate])
                if choice == 'selfing':
                    ## a random parent is chosen and its child is determined to be live or dead
                    #  if it dies, then resample it. Other wise, add its genotype to the new array.
                    parent = np.random.choice(self.population_size,1,p=parents_prob)
                    parent_genotype = self.genotypes[parent[0]]
                    #randomly generate a child 
                    child_genotype = np.random.choice(parent_genotype,2)
                    #if there's self incompatible site
                    if self.incompatible == 1:
                        if np.sum(parent_genotype) == 1:
                            if np.sum(child_genotype) == 0:
                                continue
                    new_genotypes[i] = child_genotype
                    break
                else:
                    #outcrossing: randomly choose two parents and cross them#
                    #pdb.set_trace()
                    father,mother = np.random.choice(self.population_size,2,p=parents_prob)
                    father_genotype = self.genotypes[father]
                    mother_genotype = self.genotypes[mother]
                    child_genotype = np.zeros(2)
                    child_genotype[0] = np.random.choice(father_genotype,1)
                    child_genotype[1] = np.random.choice(mother_genotype,1)
                    
                    # if the father carries the toxin while the zygote does not carry the antidote, then this child will die
                    if self.incompatible == 1:
                        if np.sum(father_genotype) >= 1 and np.sum(child_genotype) == 0:
                            continue
                    new_genotypes[i] = child_genotype
                    break
        self.genotypes = new_genotypes
                        
    def evolve_until_fix(self,max_generation=2000):
        #pdb.set_trace()
        allele_freq_record = []
        i = 0
        while self.calculate_allele_frequency() != 0 and self.calculate_allele_frequency() != 1 and i < max_generation:
            allele_freq_record.append(self.calculate_allele_frequency())
            self.one_generation()
            i += 1
        return np.array(allele_freq_record)
        
        

