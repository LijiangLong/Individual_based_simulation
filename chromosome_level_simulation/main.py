import numpy as np
import os
from multiprocessing import Pool
import multiprocessing
import pdb


class Chromosome:
    def __init__(self,A,B,recom_rate = 0.5):
        self.locus_A = A # True is functional peel-zeel
        self.locus_B = B # True if big B
        self.recom_rate = recom_rate
        
    def __repr__(self): 
        output = ''
        if self.chrI.locus_A:
            output += 'A'
        else:
            output += 'a'
        if self.chrI.locus_B:
            output += 'B'
        else:
            output += 'b'
        return output

class individual:
    def __init__(self,chrI,chrII):
        self.chrI = chrI
        self.chrII = chrII
        if self.chrI.locus_A or self.chrII.locus_A:
            self.toxin = True
        else:
            self.toxin = False
            
    def __repr__(self): 
        output = ''
        if self.chrI.locus_A:
            output += 'A'
        else:
            output += 'a'
        if self.chrI.locus_B:
            output += 'B'
        else:
            output += 'b'
        if self.chrII.locus_A:
            output += '/A'
        else:
            output += '/a'
        if self.chrII.locus_B:
            output += 'B'
        else:
            output += 'b'
        return output
        
        
        
    def generate_oocyte(self):
        recom_rate = self.chrI.recom_rate
        # first see if there's recombination
        if np.random.random() < recom_rate:
            if np.random.random() < 0.5:
                return Chromosome(self.chrI.locus_A,self.chrII.locus_B,recom_rate)
            else:
                return Chromosome(self.chrII.locus_A,self.chrI.locus_B,recom_rate)
        
        else:
            if np.random.random() < 0.5:
                return Chromosome(self.chrI.locus_A,self.chrI.locus_B,recom_rate)
            else:
                return Chromosome(self.chrII.locus_A,self.chrII.locus_B,recom_rate)
            
    def generate_sperm(self):
        return self.generate_oocyte(),self.toxin
        
    
def one_generation(incoming_groups, output_population_size):
    count = 0
    output_pop = []
    #if in a fully outcrossing population
    while True:
        father,mother = np.random.choice(incoming_groups,2,replace=False)
        chrI,toxin = father.generate_sperm()
        chrII = mother.generate_oocyte()
        new_indi = individual(chrI,chrII)
        if toxin and not new_indi.toxin:
            continue
        else:
            count += 1
            output_pop.append(new_indi)
            if count == output_population_size:
                return output_pop

def report(population):
    locusA_cnt = 0
    locusB_cnt = 0
    for indi in population:
        if indi.chrI.locus_A:
            locusA_cnt += 1
        if indi.chrI.locus_B:
            locusB_cnt += 1
        if indi.chrII.locus_A:
            locusA_cnt += 1
        if indi.chrII.locus_B:
            locusB_cnt += 1
    # print('A freq: {var1}, B freq: {var2}'.format(var1 = locusA_cnt/2/len(population),var2 = locusB_cnt/2/len(population)))
    return locusA_cnt/2/len(population),locusB_cnt/2/len(population)
        

def one_simulation(recombination_rate = 0.5):
    # pdb.set_trace()
    # simulate 1000 genomes, 100 peel-zeel and 900 non-peel-zeel
    np.random.seed()
    group_1 = [individual(Chromosome(True,True,recom_rate = recombination_rate),Chromosome(True,True,recom_rate = recombination_rate)) for i in range(400)]
    group_2 = [individual(Chromosome(False,False,recom_rate = recombination_rate),Chromosome(False,False,recom_rate = recombination_rate)) for i in range(7600)]
    individuals = group_1 + group_2
    population_size = len(individuals)
    A_freqs = []
    B_freqs = []
    freq_A, freq_B = report(individuals)
    A_freqs += [freq_A]*2
    B_freqs += [freq_B]*2
    for i in range(100):
        individuals = one_generation(individuals,population_size)
        freq_A, freq_B = report(individuals)
        A_freqs += [freq_A]
        B_freqs += [freq_B]
    return A_freqs,B_freqs



def main():
    # pdb.set_trace()

    p = Pool(multiprocessing.cpu_count()-1)
    recombination_rates = [0.5]*50+[0.4]*50 + [0.3]*50 + [0.2]*50 + [0.1]*50
    results = p.map(one_simulation,recombination_rates)
    output_file = 'differnt_recom_rate_2.txt'
    with open(output_file,'w') as output:
        for i in range(len(recombination_rates)):
            peel_trajectory,marker_trajectory = results[i]
            output.write(str(recombination_rates[i]))
            output.write('\n')
            output.write(','.join([str(x) for x in peel_trajectory]))
            output.write('\n')
            output.write(','.join([str(x) for x in marker_trajectory]))
            output.write('\n')


if __name__ == "__main__": 
    main()