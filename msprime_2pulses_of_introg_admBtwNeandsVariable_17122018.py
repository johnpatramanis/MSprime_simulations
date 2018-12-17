import msprime
import numpy as np
import math
import random as random
import sys

adm_prop_btw_neands = sys.argv[1]


####


#N_OUT_OF_AFRICA=1800
##I would recommend adding these based on the estimates given in Gravel et al 2011:
# N_AFR_atSplit = 14474 
# N_EURASIAN_atSplit = 1861 
# N_EUR_atSplit_fromEASN = 1032
# N_EASN_atSplit_fromEUR = 550
#or based on the model given in Vernot and Akey 2014 - which is a combination of Tennessen et al 2012 and Gravel et al 2011:
#N_OUT_OF_AFRICA = 1861 
#N_EASN_atSplit_fromEUR = 4257 
#N_EUR_atSplit_fromEASN = 5632


# #assuming generation time = 25
# r_EU= 0.38 #from Gravel et al 2011 
# r_ASIA= 0.48 #from Gravel et al 2011 
# r_AFRICA= 0.0172 #based on the parameters given in Vernot and Akey 2014 - exponential growth for Africans at 5115 years ago from Ne=14474 to Ne=424000 of today's Africa
# 
# generation_time = 25

# N_AFRICA=N_OG_SAPIENS / math.exp(-r_AFRICA * T_split_AFRICA_OUTOFAFRICA)
# N_EU=N_EU0 / math.exp(-r_EU * T_split_EU_ASIA)
# N_ASIA=N_ASIA0 / math.exp(-r_ASIA * T_split_EU_ASIA)


##assuming generation time = 29 - multiplying above values by  29/25 
r_EU = 0.0226197294141683734 #corrected exponential growth rates
r_ASIA= 0.028572397366705403
r_AFRICA= 0.019148313407347897

generation_time = 29.0

T_split_EU_ASIA_sample = np.array(np.random.uniform(low = 36000, high = 55000, size=1))
#N_OUT_OF_AFRICA_sample = np.array(np.random.randint(low = 5000, high = 15000, size=1))
T_introgration1_sample = np.array(np.random.uniform(low = T_split_EU_ASIA_sample, high = 65000, size=1))
adm_prop2_N_EA_sample = np.array(np.random.uniform(low = -0.000015, high = 0.00025, size=1))
T_introgression_btw_neands_sample = np.array(np.random.uniform(low = int(T_introgration1_sample), high = 130000, size=1))

T_split_SAPIENS_NEAD_DENI= 700000/generation_time
T_split_NEAD_DENI= 400000/generation_time
T_split_NEAD1_NEAD2= 130000/generation_time
T_split_AFRICA_OUTOFAFRICA= 70000/generation_time
T_split_EU_ASIA= T_split_EU_ASIA_sample/generation_time

T_introgration1=T_introgration1_sample/generation_time
T_introgration1_end=(T_introgration1_sample+500.0)/generation_time

T_introgration2= (T_split_EU_ASIA_sample+500.0)/generation_time
T_introgration2_end = (T_split_EU_ASIA_sample+1000.0)/generation_time

T_archaic_sampling=50000/generation_time
T_archaic_sampling_neand2=50000/generation_time

T_introgression_btw_neands=T_introgression_btw_neands_sample/generation_time

#If we try to use the parameters given in Gravel et al 2011, I would recommend using this time point where European and East Asian populations start growing exponentially.
T_exp_growth = 5115/generation_time
###


N_OUT_OF_AFRICA = 1861 
N_EUR_atSplit_EU_ASIA = 1032
N_EASN_atSplit_EU_ASIA = 550

N_OG_SAPIENS=14474
N_OG_NEAD_DENI=10000
N_NEAD1=1000
N_NEAD2=1000
N_DENI=1000
N_EUR_atSplit_EU_ASIA = 1032
N_EASN_atSplit_EU_ASIA = 550
N_AFR_beforeExpGrowth = 14474
N_EUR_beforeExpGrowth = 9475
N_EASN_beforeExpGrowth = 8879


#I would recommend the following: based on the model given in Vernot and Akey 2014 (time at which growth starts taken from Gravel et al 2011)
N_AFRICA= N_AFR_beforeExpGrowth / math.exp(-r_AFRICA * T_exp_growth) #at the time = 5115 years ago
N_EU= N_EUR_beforeExpGrowth / math.exp(-r_EU * T_exp_growth) #at the time = 5115 years ago
N_ASIA= N_EASN_beforeExpGrowth / math.exp(-r_ASIA * T_exp_growth) #at the time = 5115 years ago

r_EU_btw_ExpGrowth_EuEasplit = float(math.log(N_EUR_beforeExpGrowth) - math.log(N_EUR_atSplit_EU_ASIA)) / (T_split_EU_ASIA-T_exp_growth)
r_EA_btw_ExpGrowth_EuEasplit = float(math.log(N_EASN_beforeExpGrowth) - math.log(N_EASN_atSplit_EU_ASIA)) / (T_split_EU_ASIA-T_exp_growth)



population_configurations = [
    msprime.PopulationConfiguration(initial_size=N_DENI),
    msprime.PopulationConfiguration(initial_size=N_NEAD1),
    msprime.PopulationConfiguration(initial_size=N_NEAD2),
    msprime.PopulationConfiguration(initial_size=N_AFRICA,growth_rate = r_AFRICA),
    msprime.PopulationConfiguration(initial_size=N_EU, growth_rate = r_EU),
    msprime.PopulationConfiguration(initial_size=N_ASIA, growth_rate = r_ASIA)
]

#migrations
m_=0.00001
#m_T1=131

m_eu_ea = 0.000025 * 23000.0/float(T_split_EU_ASIA_sample)
m_afr_ea = 0.0000078 * 23000.0/float(T_split_EU_ASIA_sample)
m_afr_eu =  0.0000311 * 23000.0/float(T_split_EU_ASIA_sample)



migration_matrix = [
####DE,NE1,NE2,AF,EU,AS##
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,m_afr_eu,m_afr_ea],
    [0,0,0,m_afr_eu,0,m_eu_ea],
    [0,0,0,m_afr_ea,m_eu_ea,0]
]

samples=[msprime.Sample(0,T_archaic_sampling)]*6 + [msprime.Sample(1,T_archaic_sampling)]*6 + [msprime.Sample(2,T_archaic_sampling_neand2)] *6 + [msprime.Sample(3,0)] *200 + [msprime.Sample(4,0)] *550 + [msprime.Sample(5,0)] *550


demographic_events = [
######################################################################################
#Exponential growth begins


msprime.PopulationParametersChange(time =T_exp_growth , growth_rate = 0 , population_id = 3),

msprime.PopulationParametersChange(time =T_exp_growth , growth_rate = r_EU_btw_ExpGrowth_EuEasplit , population_id = 4),
msprime.PopulationParametersChange(time =T_exp_growth ,growth_rate=r_EA_btw_ExpGrowth_EuEasplit , population_id = 5),

#######################################################################################
#Split of Eu and Asia ,bottleneck,steady growth 


msprime.PopulationParametersChange(time =T_split_EU_ASIA , growth_rate = 0 , population_id = 4),
msprime.PopulationParametersChange(time =T_split_EU_ASIA ,growth_rate=0 , population_id = 5),

msprime.PopulationParametersChange(time =T_split_EU_ASIA , initial_size = N_EUR_atSplit_EU_ASIA , population_id = 4),
msprime.PopulationParametersChange(time =T_split_EU_ASIA , initial_size = N_EASN_atSplit_EU_ASIA, population_id = 5),



msprime.MassMigration(time=T_split_EU_ASIA,source=5,destination=4,proportion = 1.0),

#this seems to be redundant - or perhaps wrong - am I missing something?
msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0), ##didn't know how to change here!

##does this look correct here? - wanted to set a new migration rate between the ancestors of Eurasians and Africans
msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0.0001498975, matrix_index=(4,3)), 
msprime.MigrationRateChange(time=T_split_EU_ASIA,rate=0.0001498975, matrix_index=(3,4)),
    
########################################################################################

#Introgration
msprime.MigrationRateChange(time=T_introgration1,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1,rate=0, matrix_index=(4, 1)),


msprime.MigrationRateChange(time=T_introgration1,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1,rate=0.00075, matrix_index=(4, 1)),

#End of Introgration

msprime.MigrationRateChange(time=T_introgration1_end,rate=0, matrix_index=(1, 4)),
msprime.MigrationRateChange(time=T_introgration1_end,rate=0, matrix_index=(4, 1)),


#msprime.MassMigration(time=T_introgration1,source=1,destination=4,proportion = 0.03),
########################################################################################

#Second Introgration into EA
msprime.MigrationRateChange(time=T_introgration2,rate=0, matrix_index=(2, 5)),
msprime.MigrationRateChange(time=T_introgration2,rate=0, matrix_index=(5, 2)),


msprime.MigrationRateChange(time=T_introgration2,rate=0, matrix_index=(2, 5)),
msprime.MigrationRateChange(time=T_introgration2,rate=adm_prop2_N_EA_sample, matrix_index=(5, 2)),

#End of Introgration

msprime.MigrationRateChange(time=T_introgration2_end,rate=0, matrix_index=(2, 5)),
msprime.MigrationRateChange(time=T_introgration2_end,rate=0, matrix_index=(5, 2)),


########################################################################################
#Split Africa-Eurasian, population 3 size set to homo sapiens

msprime.PopulationParametersChange(time =T_split_AFRICA_OUTOFAFRICA , initial_size = N_OUT_OF_AFRICA , growth_rate = 0, population_id = 4),

msprime.MassMigration(time=T_split_AFRICA_OUTOFAFRICA,source=4,destination=3,proportion = 1.0),
msprime.MigrationRateChange(time=T_split_AFRICA_OUTOFAFRICA,rate=0),
    
msprime.PopulationParametersChange(time =T_split_AFRICA_OUTOFAFRICA , initial_size = N_OG_SAPIENS , growth_rate = 0, population_id = 4),

########################################################################################
#Admixture between 2 neanderthal lineages
#msprime.MassMigration(time=T_introgression_btw_neands,source=2,destination=1,proportion = float(adm_prop_btw_neands)/100.0),
#msprime.MassMigration(time=T_introgression_btw_neands,source=1,destination=2,proportion = float(adm_prop_btw_neands)/100.0),

msprime.MigrationRateChange(time=T_introgression_btw_neands,rate=0, matrix_index=(1, 2)),
msprime.MigrationRateChange(time=T_introgression_btw_neands,rate=0, matrix_index=(2, 1)),


msprime.MigrationRateChange(time=T_introgression_btw_neands,rate=float(adm_prop_btw_neands)/100.0, matrix_index=(1, 2)),
msprime.MigrationRateChange(time=T_introgression_btw_neands,rate=float(adm_prop_btw_neands)/100.0, matrix_index=(2, 1)),

#End of Introgration

msprime.MigrationRateChange(time=T_introgression_btw_neands+1,rate=0, matrix_index=(1, 2)),
msprime.MigrationRateChange(time=T_introgression_btw_neands+1,rate=0, matrix_index=(2, 1)),


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

#recomb_map=msprime.RecombinationMap.read_hapmap('genetic_map_GRCh37_chr{}.txt'.format(chrom))

#Number of replications, not all reps have an introgression
n_replicates=1000
LENGTH=10e+5
random_seed=random.randint(0,100000)

#DEMOGRAPHY DEBUGGER -> un-hash this part to print the model
#dd = msprime.DemographyDebugger(
#    population_configurations=population_configurations,
#    migration_matrix=migration_matrix,
#    demographic_events=demographic_events)
#dd.print_history()
    
    
# Hash it if you want to see only the demography debugger
#The actual simulation begins, all info is stored in the dd object
dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=1e-8,
    demographic_events=demographic_events,record_migrations=True,random_seed=random_seed,length=LENGTH, recombination_rate=2e-8 ,num_replicates=n_replicates)
#recombination_map=recomb_map


chrom=1
#Output file 1 
out=open('Info_{}_2pulse_AdmPropBtwNs_{}perc.gen'.format(chrom,adm_prop_btw_neands),'w')
out2=open('genotypes_{}_2pulse_AdmPropBtwNs_{}perc.gen'.format(chrom,adm_prop_btw_neands),'w')
out3=open('Patterns_{}_2pulse_AdmPropBtwNs_{}perc.gen'.format(chrom,adm_prop_btw_neands),'w')
out4=open('Introgressed_variants_{}_2pulse_AdmPropBtwNs_{}perc.gen'.format(chrom,adm_prop_btw_neands),'w')


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
    random_seed=random.randint(0,100000)

    
    #Cycle through all migrations that happened
    for i in j.migrations():
        
        #Only interested of migrations from Pop1 (NEAN1) to Pop4 (Asia but also Eurasia ancest pop )
        if i.source==4 and i.dest==1:
            if ( str(i.left) not in leftchunks ) and ( str(i.right) not in rightchunks ): #if not yet in our collection, add  it
                leftchunks.append(float(i.left))   #original introgressed chunk
                rightchunks.append(float(i.right)) #original introgressed chunk
                who.append(int(i.node))          #original introgressed individual
                when.append(str(i.time))         # time of introgretion
            if int(i.node) not in introgressed: #all original introgressed individuals
                introgressed.append(int(i.node))
                #break
    for p in j.variants():
        positions.append([p.position,p.index])
        #all variants position and their number, we will need it later to see wich variants belong to neand
    
    
    
    
    RecombinatingSegments={}
    RC=0
    for tree in j.trees():
        INTS=0
        WHO=[]
        #Now we cycle through all recombinating segments of our sequence
        #we can print the tree/history of each sement
        #print(tree.draw(format="unicode"))
        
        #We wan to see for each segment where the children of the original introgressed segment went (eg through recombination some segments can migrate to other individuals)
        for k in range(0,len(who)):
            if leftchunks[k]<=tree.interval[0] and tree.interval[1]<=rightchunks[k]:
                if tree.is_leaf(who[k])!=True:
                    for my in tree.get_leaves(who[k]):
                        if int(my) not in trueintrogressed:
                            trueintrogressed.append(int(my))
                else:
                    if int(who[k]) not in trueintrogressed:
                        trueintrogressed.append(int(who[k]))
    
                        
            #keep only the modern people that have introgressedsegments       
        trueintrogressedfinal=[x for x in trueintrogressed if x in modernsamples]
    

        # cycle through modern people that have an introgressed segment and see if this particular we are now (for tree in j.trees)
        # is introgressed in each individual
        if trueintrogressed!=[]:
            for ind in trueintrogressedfinal:
                
                #print('ind-nead',np.mean([tree.tmrca(ind,x) for x in range(6,11)]))
                #print('ind-afric',np.mean([tree.tmrca(ind,x) for x in range(18,23)]))


                if tree.population(tree.mrca(ind,8))==0 :  #if the mrca of th individual and a Neanderthal (8) belongs to population 0 then its not an introgressed segment   
                    pass    
                
                if tree.population(tree.mrca(ind,8))!=0 and tree.population(tree.mrca(ind,8))!=1 and tree.population(tree.mrca(ind,8))!=2 :
                    print('anomaly')
                
                if tree.population(tree.mrca(ind,8))==1 or tree.population(tree.mrca(ind,8))==2:#if it belongs to population 1 or 2 (Neaderthal pops) then its an introgressed seg
                    INTS+=1
                    WHO.append(int(ind))
                    try:
                        person_introgretions[ind]+=1
                    except KeyError:
                        person_introgretions[ind]=1
                    try:
                        person_chunks[ind].append(tree.interval)
                    except KeyError:
                        person_chunks[ind]=[ tree.interval ]
        RecombinatingSegments[RC]=[INTS,float(tree.interval[0]),float(tree.interval[1]),float(tree.interval[1]-tree.interval[0]),WHO]
        RC+=1
    
    
    
    
    
    
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
        out.write('The individual {} has {} introgressed trees in him/her \n'.format(k,person_introgretions[k]/j.num_trees))
    
    
    numberoftrees=j.num_trees
    print('modern descendants of introgressed nodes: ',trueintrogressedfinal)
    out.write('modern descendants of introgressed nodes: {} \n'.format(trueintrogressedfinal))
    
    for r in range(0,len(leftchunks)):
        print('From {} to {} was an introgression of the individual {} at time = {}'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
        out.write('From {} to {} was an introgression of the individual {} at time = {}\n'.format(leftchunks[r],rightchunks[r],who[r],when[r]))
    
    for chunk in  person_chunks:
        introgressionsum=0
        for i in person_chunks[chunk]:
            add=i[1]-i[0]
            introgressionsum+=add
        print('Person {} has {} chunks of introgretion in him/her '.format(chunk,introgressionsum/LENGTH))    
        out.write('Person {} has {} chunks of introgretion in him/her \n'.format(chunk,introgressionsum/LENGTH))
    
    out.write('#######################\n')
    if trueintrogressedfinal!=[] and len(trueintrogressedfinal)>=3:
        for var in j.variants():
            for indvar in var.genotypes:
                out2.write(str(indvar))
            out2.write('\n')
        out2.write('#################################################################################################################################\n')
        
        for var in j.variants():
            for R in RecombinatingSegments:
                #print(RecombinatingSegments[R][1],RecombinatingSegments[R][2],var.site.position)
                if RecombinatingSegments[R][1]<=float(var.site.position) and float(var.site.position)<=RecombinatingSegments[R][2] and RecombinatingSegments[R][0]>0:
                    for variant in var.genotypes:
                        out4.write(str(variant))
                    out4.write('\t{}'.format(var.site.position))
                    out4.write('\t')
                    out4.write(','.join([str(x) for x in RecombinatingSegments[R][4]]))
                    out4.write('\n')
        out4.write('#####################################################################\n')
        for RC in RecombinatingSegments:
            out3.write('{} {} {} {} \n'.format(RecombinatingSegments[RC][0],RecombinatingSegments[RC][1],RecombinatingSegments[RC][2],RecombinatingSegments[RC][3]))
            
        
        out3.write('#####################################################################\n')
        
    print(len(person_chunks),len(trueintrogressedfinal))
