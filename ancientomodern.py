import msprime
import numpy as np
import math



N_AF=10000
N_OUT=1000
N_Ancient=1000






r_AF=0.00
r_Ancient=0.001
r_Modern=0.001


generation_time = 23


T_split_AFRICA_OUTOFAFRICA= 80000/generation_time
T_split_Ancient_Modern= 4500/generation_time
T_Growth_Begins=700/generation_time
T_Serious_Growth=100/generation_time


T_archaic_sampling=4000/generation_time



########################## Modern population size with 2 different growth rates ??
N_Modern=N_OUT / math.exp(-r_Modern * T_split_AFRICA_OUTOFAFRICA)




population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_DENI),
    msprime.PopulationConfiguration(initial_size=N_NEAD1),
    msprime.PopulationConfiguration(initial_size=N_NEAD2),
    msprime.PopulationConfiguration(initial_size=N_AFRICA),
    msprime.PopulationConfiguration(initial_size=N_EU),
    msprime.PopulationConfiguration(initial_size=N_ASIA)
]

#migrations
m_=0.00001
#m_T1=131



migration_matrix = [
####DE,NE1,NE2,AF,EU,AS##
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0.01,0.001],
    [0,0,0,0.01,0,0.01],
    [0,0,0,0.001,0.01,0]
]

samples=[msprime.Sample(0,T_archaic_sampling)]*6 + [msprime.Sample(1,T_archaic_sampling)]*6 + [msprime.Sample(2,T_archaic_sampling)] *6 + [msprime.Sample(3,0)] *6 + [msprime.Sample(4,0)] *6 + [msprime.Sample(5,0)] *6


demographic_events = [
######################################################################################
#Split EU-AS ,migrations set to 0, growth set to 0, population 4 size set to out_of_africa

msprime.MassMigration(time=T_split_EU_ASIA,source=5,destination=4,proportion = 1.0),

msprime.PopulationParametersChange(time =T_split_EU_ASIA , initial_size = N_OUT_OF_AFRICA , growth_rate = 0, population_id = 4),
msprime.PopulationParametersChange(time =T_split_EU_ASIA , growth_rate = 0, population_id = 5),

msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0),

########################################################################################

#Introgration

msprime.MigrationRateChange(time=T_introgration1,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1,rate=0.02, matrix_index=(4, 1)),

#End of Introgration

msprime.MigrationRateChange(time=T_introgration1+5,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1+5,rate=0, matrix_index=(4, 1)),




########################################################################################
#Split Africa-Eurasian, population 3 size set to homo sapiens

msprime.MassMigration(time=T_split_AFRICA_OUTOFAFRICA,source=4,destination=3,proportion = 1.0),

msprime.PopulationParametersChange(time =T_split_AFRICA_OUTOFAFRICA , initial_size = N_OG_SAPIENS , growth_rate = 0, population_id = 4),

########################################################################################
#Split of 2 neanderthal lineages

msprime.MassMigration(time=T_split_NEAD1_NEAD2,source=1,destination=2,proportion = 1.0),
########################################################################################
#Split of Neanderthal and Denisova, population 0 size set to Neanderthal-Denisova

msprime.MassMigration(time=T_split_NEAD_DENI,source=2,destination=0,proportion = 1.0),
msprime.PopulationParametersChange(time =T_split_NEAD_DENI , initial_size = N_OG_NEAD_DENI , population_id = 0),

########################################################################################
#Split of sapiens with archaic

msprime.MassMigration(time=T_split_SAPIENS_NEAD_DENI,source=3,destination=0,proportion = 1.0)

]

chrom=1

#recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(chrom))

n_replicates=1000
random_seed=random.randint(0,100000)

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1e-8,
    demographic_events=demographic_events,record_migrations=True,length=1000000,random_seed=random_seed,recombination_rate=2e-8 ,num_replicates=n_replicates)
#recombination_map=recomb_map
