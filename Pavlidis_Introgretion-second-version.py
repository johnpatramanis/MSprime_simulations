import msprime
import numpy as np
import math
import random as random

#Population sizes

N_OG_SAPIENS=10000
N_OG_NEAD_DENI=10000
N_NEAD1=5000
N_NEAD2=5000
N_DENI=5000
N_OUT_OF_AFRICA=500
N_EU0=500
N_ASIA0=500


#time of bottlenecks


#Growth rate
r_EU=0.01
r_ASIA=0.01
r_AFRICA=0.001


#Generation time
generation_time = 29



#Time is calculated in generations so Years/Generation time = Generations
#Population Splits time
T_split_SAPIENS_NEAD_DENI= 700000/generation_time
T_split_NEAD_DENI= 350000/generation_time
T_split_NEAD1_NEAD2= 130000/generation_time
T_split_AFRICA_OUTOFAFRICA= 80000/generation_time
T_split_EU_ASIA= 40000/generation_time

#Time of introgression,time of archaic sampling
T_introgration1=50000/generation_time
T_archaic_sampling=50000/generation_time



#Modern populations size
N_AFRICA=N_OG_SAPIENS / math.exp(-r_AFRICA * T_split_AFRICA_OUTOFAFRICA)
N_EU=N_EU0 / math.exp(-r_EU * T_split_EU_ASIA)
N_ASIA=N_ASIA0 / math.exp(-r_ASIA * T_split_EU_ASIA)



#Set up our populations
population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_DENI),
    msprime.PopulationConfiguration(initial_size=N_NEAD1),
    msprime.PopulationConfiguration(initial_size=N_NEAD2),
    msprime.PopulationConfiguration(initial_size=N_AFRICA),
    msprime.PopulationConfiguration(initial_size=N_EU),
    msprime.PopulationConfiguration(initial_size=N_ASIA)
]



#Our migration matrix
migration_matrix = [
####DE,NE1,NE2,AF,EU,AS##
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0]
]

#How many samples from each population, when msprime.Sample(whichpop, when)
samples=[msprime.Sample(0,T_archaic_sampling)]*6 + [msprime.Sample(1,T_archaic_sampling)]*6 + [msprime.Sample(2,T_archaic_sampling)] *6 + [msprime.Sample(3,0)] *6 + [msprime.Sample(4,0)] *20 + [msprime.Sample(5,0)] *20


demographic_events = [
######################################################################################
#Split EU-AS ,migrations set to 0, growth set to 0, population 4 size set to out_of_africa
#msprime.MigrationRateChange(time=T_split_EU_ASIA-10,rate=0.01, matrix_index=(4, 3)),
#msprime.MigrationRateChange(time=T_split_EU_ASIA-10,rate=0.01, matrix_index=(5, 3)),


msprime.MassMigration(time=T_split_EU_ASIA,source=5,destination=4,proportion = 1.0),

msprime.PopulationParametersChange(time =T_split_EU_ASIA , initial_size = N_OUT_OF_AFRICA , growth_rate = 0, population_id = 4),
msprime.PopulationParametersChange(time =T_split_EU_ASIA , growth_rate = 0, population_id = 5),

msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0),

########################################################################################

#Introgration

msprime.MigrationRateChange(time=T_introgration1,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1,rate=0.01, matrix_index=(4, 1)),

#End of Introgration

msprime.MigrationRateChange(time=T_introgration1+1,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1+1,rate=0, matrix_index=(4, 1)),




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
# I have successfully used the following to load the recombination map of chromosome 1 so the length and recombination rates parameters
# of the Simulation match chromosome (can be one with the others as well0 but had 'memory error' problems

recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(chrom))

#Number of replications, not all reps have an introgression
n_replicates=1000
LENGTH=10e+5
random_seed=random.randint(0,100000)

#The actual simulation begins, all info is stored in the dd object
dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1e-8,
    demographic_events=demographic_events,record_migrations=True,random_seed=random_seed,length=LENGTH, recombination_rate=2e-8 ,num_replicates=n_replicates)
#recombination_map=recomb_map

events=0
chrom=1
#Output file 1 
out=open('Info_{}.gen'.format(chrom),'w')
out2=open('genotypes_{}.gen'.format(chrom),'w')
for j in dd:
    
    introgressed=[]
    trueintrogressed=[]
    leftchunks=[]
    rightchunks=[]
    who=[]
    when=[]
    positions=[]
    person_introgretions={}
    person_chunks={}
    modernsamples=[i for i in range(0,j.num_samples)]
    
    #Cycle through all migrations that happened
    for i in j.migrations():
        
        #Only interested of migrations from Pop1 (NEAN1) to Pop4 (Asia but also Eurasia ancest pop )
        if i.source==4 and i.dest==1:
            if ( str(i.left) not in leftchunks ) and ( str(i.right) not in rightchunks ): #if not yet in our collection, add  it
                leftchunks.append(str(i.left))   #original introgressed chunk
                rightchunks.append(str(i.right)) #original introgressed chunk
                who.append(str(i.node))          #original introgressed individual
                when.append(str(i.time))         # time of introgretion
            if int(i.node) not in introgressed: #all original introgressed individuals
                introgressed.append(int(i.node))
                events+=1
                #break
    for p in j.variants():
        positions.append([p.position,p.index])
        #all variants position and their number, we will need it later to see wich variants belong to neand
    
    
    
    
    for tree in j.trees():
        #Now we cycle through all recombinating segments of our sequence
        #we can print the tree/history of each sement
        #print(tree.draw(format="unicode"))
        
        #We wan to see for each segment where the children of the original introgressed segment went (eg through recombination some segments can migrate to other individuals)
        for k in introgressed:
            if tree.is_leaf(k)!=True:
                for my in tree.get_leaves(k):
                    if int(my) not in trueintrogressed:
                        trueintrogressed.append(int(my))
            else:
                if int(k) not in trueintrogressed:
                    trueintrogressed.append(int(k))

                    
        #keep only the modern people that have introgressedsegments       
        trueintrogressedfinal=[x for x in trueintrogressed if x in modernsamples]
    

        # cycle through modern people that have an introgressed segment and see if this particular we are now (for tree in j.trees)
        # is introgressed in each individual
        if trueintrogressed!=[]:
            for ind in trueintrogressedfinal:
                
                #print('ind',tree.tmrca(ind,8))
                #print('african',tree.tmrca(20,8))

                if tree.tmrca(ind,8) == tree.tmrca(20,8):  # 50000 is the year we have set for the introgretion event, note that for some reason here tmrca return years not generations   
                    pass    
                
                if tree.tmrca(ind,8) < tree.tmrca(20,8):
                    try:
                        person_introgretions[ind]+=1                                        
                    except KeyError:
                        person_introgretions[ind]=1
                    try:
                        person_chunks[ind].append(tree.interval)
                    except KeyError:
                        person_chunks[ind]=[ tree.interval ]
            
    for i in person_chunks:
        temporary_person_chunks=[]
        temporary_person_chunks.append(list(person_chunks[i][0]))
        for k in person_chunks[i][1:]:
            if k[0]==temporary_person_chunks[-1][1]:
                temporary_person_chunks[-1][1]=k[1]
            else:
                temporary_person_chunks.append(list(k))
                
        person_chunks[i]=temporary_person_chunks
    
            






            
            
        #print('recombinant: ',tree.interval)
    #print(j.genotype_matrix())
    print('#######################','\n')
    out.write('#######################\n')
    print('There are ',len(positions),' variants')
    
    
    
    
    
    print('original introgressed nodes: ',introgressed)
    out.write('original introgressed nodes: {}\n'.format(introgressed))
    
    print('number of trees (recombinations) : ',j.num_trees)
    out.write('number of trees (recombinations) : {}\n'.format(j.num_trees))
   
    for k in person_introgretions:
        print('The individual {} has {} introgressed trees in him/her'.format(k,person_introgretions[k]/j.num_trees))
    
    
    
    numberoftrees=j.num_trees
    print('modern introgressed nodes: ',trueintrogressedfinal)
    out.write('modern introgressed nodes: {} \n'.format(trueintrogressedfinal))
    
    for r in range(0,len(leftchunks)):
        print('From {} to {} was an introgression of the individual {} at time = {}'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
        out.write('From {} to {} was an introgression of the individual {} at time = {}\n'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
    
    for chunk in  person_chunks:
        introgressionsum=0
        for i in person_chunks[chunk]:
            add=i[1]-i[0]
            introgressionsum+=add
        print('Person {} has {} chunks of introgretion in him/her '.format(chunk,introgressionsum/LENGTH))    
    

    
    out.write('#######################\n')
    
    
    
print(events/n_replicates)
