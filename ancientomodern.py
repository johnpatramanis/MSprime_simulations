import msprime
import numpy as np
import math
import random


N_OUTGROUP=10000
N_ANCIENT=1000
N0_MODERN=1000






r_OUT=0.00
r_ANCIENT=0
r_MODERN=0.001


generation_time = 23


T_split_OUTGROUP= 10000/generation_time
T_split_Ancient_Modern= 2500/generation_time
T_Growth_Begins=700/generation_time
T_Serious_Growth=100/generation_time


T_archaic_sampling=2500/generation_time



########################## Modern population size with 2 different growth rates ??
N_MODERN=N0_MODERN / math.exp(-r_MODERN * T_split_Ancient_Modern)




population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_OUTGROUP),
    msprime.PopulationConfiguration(initial_size=N_ANCIENT),
    msprime.PopulationConfiguration(initial_size=N_MODERN)
]

#migrations
m_=0.00001
#m_T1=131



migration_matrix = [
####DE,NE1,NE2,AF,EU,AS##
    [0,0,0],
    [0,0,0],
    [0,0,0]
]

samples=[msprime.Sample(0,0)]*100 + [msprime.Sample(1,T_archaic_sampling)]*100 + [msprime.Sample(2,0)]*100 


demographic_events = [
######################################################################################
#

msprime.MassMigration(time=T_split_Ancient_Modern,source=2,destination=1,proportion = 1.0),
msprime.PopulationParametersChange(time =T_split_Ancient_Modern , initial_size = N0_MODERN , growth_rate = 0, population_id = 1),

#msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0),

########################################################################################
########################################################################################
#

msprime.MassMigration(time=T_split_OUTGROUP,source=0,destination=1,proportion = 1.0),

#msprime.PopulationParametersChange(time =T_split_AFRICA_OUTOFAFRICA , initial_size = N_OG_SAPIENS , growth_rate = 0, population_id = 4)

########################################################################################
]

chrom=1

#recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(chrom))

n_replicates=1
random_seed=random.randint(0,100000)

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1e-8,
    demographic_events=demographic_events,record_migrations=True,length=1000000,random_seed=random_seed,recombination_rate=2e-8 ,num_replicates=n_replicates)
#recombination_map=recomb_map
