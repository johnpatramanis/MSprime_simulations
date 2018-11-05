import msprime
import numpy as np
import math




N_OG_SAPIENS=10000
N_OG_NEAD_DENI=10000
N_NEAD1=5000
N_NEAD2=5000
N_DENI=5000
N_OUT_OF_AFRICA=500
N_EU0=250
N_ASIA0=250

r_EU=0.004
r_ASIA=0.004
r_AFRICA=0.001


generation_time = 29


T_split_SAPIENS_NEAD_DENI= 700000/generation_time
T_split_NEAD_DENI= 350000/generation_time
T_split_NEAD1_NEAD2= 130000/generation_time
T_split_AFRICA_OUTOFAFRICA= 80000/generation_time
T_split_EU_ASIA= 40000/generation_time

T_introgration1=50000/generation_time
T_archaic_sampling=50000/generation_time




N_AFRICA=N_OG_SAPIENS / math.exp(-r_AFRICA * T_split_AFRICA_OUTOFAFRICA)
N_EU=N_EU0 / math.exp(-r_EU * T_split_EU_ASIA)
N_ASIA=N_ASIA0 / math.exp(-r_ASIA * T_split_EU_ASIA)

n_replicates=1000


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




migration_matrix = [
####DE,NE1,NE2,AF,EU,AS##
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0.01,0.001],
    [0,0,0,0.01,0,0.01],
    [0,0,0,0.001,0.01,0]
]

samples=[msprime.Sample(0,T_archaic_sampling)]*2 + [msprime.Sample(1,T_archaic_sampling)]*2 + [msprime.Sample(2,T_archaic_sampling)] *2 + [msprime.Sample(3,0)] *2 + [msprime.Sample(4,0)] *2 + [msprime.Sample(5,0)] *2


demographic_events = [
######################################################################################
#Split EU-AS ,migrations set to 0, growth set to 0, population 4 size set to out_of_africa

msprime.MassMigration(time=T_split_EU_ASIA,source=5,destination=4,proportion = 1.0),

msprime.PopulationParametersChange(time =T_split_EU_ASIA , initial_size = N_OUT_OF_AFRICA , growth_rate = 0, population_id = 4),
msprime.PopulationParametersChange(time =T_split_EU_ASIA , growth_rate = 0, population_id = 5),

msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0),

########################################################################################
#Introgration

msprime.MassMigration(time=T_introgration1,source=4,destination=1,proportion = 0.2),

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

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1.25e-8,
    demographic_events=demographic_events,length=1.0,record_migrations=True, recombination_rate=2e-8,num_replicates=n_replicates)

#for i in dd.trees():
#    print(i.draw(format="unicode"))

events=0
for j in dd:
    for i in j.migrations():
        if i.source==4 and i.dest==1:
            events+=1
            break
print(events/n_replicates)