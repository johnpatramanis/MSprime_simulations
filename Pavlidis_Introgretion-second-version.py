import msprime
import numpy as np
import math
import random as random



N_OG_SAPIENS=10000
N_OG_NEAD_DENI=10000
N_NEAD1=5000
N_NEAD2=5000
N_DENI=5000
N_OUT_OF_AFRICA=500
N_EU0=250
N_ASIA0=250


#bottleneck



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

events=0
chrom=1
out=open('genotypes{}.gen'.format(chrom),'w')
for j in dd:
    
    introgressed=[]
    trueintrogressed=[]
    leftchunks=[]
    rightchunks=[]
    who=[]
    when=[]
    
    for i in j.migrations():
        
        if i.source==4 and i.dest==1:
            if ( str(i.left) not in leftchunks ) and ( str(i.right) not in rightchunks ):
                leftchunks.append(str(i.left))
                rightchunks.append(str(i.right))
                who.append(str(i.node))
                when.append(str(i.time))
            if int(i.node) not in introgressed:
                introgressed.append(int(i.node))
                events+=1
                #break
    for tree in j.trees():
        #print(tree.draw(format="unicode"))
        for k in introgressed:
            if tree.is_leaf(k)!=True:
                for my in tree.get_leaves(k):
                    if int(my) not in trueintrogressed:
                        trueintrogressed.append(int(my))
            else:
                if int(k) not in trueintrogressed:
                    trueintrogressed.append(int(k))
                
                
                
                
          

        #print('recombinant: ',tree.interval)
    #print(j.genotype_matrix())
    print('#######################','\n')
    out.write('#######################\n')
    
    modernsamples=[i for i in range(0,j.num_samples)]
    trueintrogressedfinal=[x for x in trueintrogressed if x in modernsamples]
    print('original introgressed nodes: ',introgressed)
    out.write('original introgressed nodes: {}\n'.format(introgressed))
    
    print('number of trees (recombinations) : ',j.num_trees)
    out.write('number of trees (recombinations) : {}\n'.format(j.num_trees))
    
    numberoftrees=j.num_trees
    print('modern introgressed nodes: ',trueintrogressedfinal)
    out.write('modern introgressed nodes: {} \n'.format(trueintrogressedfinal))
    
    for r in range(0,len(leftchunks)):
        print('From {} to {} was an introgression of the individual {} at time = {}'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
        out.write('From {} to {} was an introgression of the individual {} at time = {}\n'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
    print('The total naumber of introgressed segments is {}'.format(len(leftchunks)))
    out.write('The total naumber of introgressed segments is {} \n'.format(len(leftchunks)))
    out.write('#######################\n')
    
    
    
print(events/n_replicates)
