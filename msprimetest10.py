import msprime
import numpy as np
import math

generation=29

r_AMBRACIA=0

T_PROTOEUROPE=2970/generation #3000
T_PROTOGREEK=1170/generation #1200 BC
T_FOUNDING=600/generation #650 BC
T_SAMPLEING1=440/generation #479BC
T_SAMPLING2=0 #30BC - TODAY

N_PROTOEUROPE=500
N_PROTOGREEK=500
N_CORINTH=300
N0_AMBRACIA=50
N_AMBRACIA= N0_AMBRACIA/ math.exp(-r_AMBRACIA * T_FOUNDING)


population_configurations =[
msprime.PopulationConfiguration(initial_size=N_PROTOEUROPE),
msprime.PopulationConfiguration(initial_size=N_PROTOGREEK),
msprime.PopulationConfiguration(initial_size=N_CORINTH),
msprime.PopulationConfiguration(initial_size=N_AMBRACIA,growth_rate=r_AMBRACIA)]


migration_matrix= [
    [0,0,0,0],
    [0,0,0,0.3],
    [0,0,0,0],
    [0,0.3,0,0]]



samples=[msprime.Sample(1,0)]*100 + [msprime.Sample(2,0)]*100 + [msprime.Sample(3,0)] *100 + [msprime.Sample(3,440/generation)] *100

demographic_events = [
    #FOUNDING OF AMBRACIA
    msprime.MigrationRateChange(time=T_FOUNDING,rate=0),
    msprime.PopulationParametersChange(time=T_FOUNDING,initial_size=N0_AMBRACIA,growth_rate=0,population_id=3),
    msprime.MassMigration(time=T_FOUNDING, source=3, destination=2, proportion=1.0),
    #Dorian Invasion
    msprime.MassMigration(time=T_PROTOGREEK, source=2, destination=1, proportion=1.0),
    #Split wth other proto-europeans
    msprime.MassMigration(time=T_PROTOEUROPE, source=1, destination=0, proportion=1.0),
]

dd = msprime.DemographyDebugger(
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,
    demographic_events=demographic_events)
dd.print_history()

#dd = msprime.simulate(samples=samples,
#    population_configurations=population_configurations,
#    migration_matrix=migration_matrix,mutation_rate=0.01,
#    demographic_events=demographic_events,length=1.0, recombination_rate=2e-8)
#for wow in dd.trees():
#    print(wow.draw(format="unicode"))

dd = msprime.simulate(samples=samples,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,mutation_rate=0.00000004,
    demographic_events=demographic_events,length=54545, recombination_rate=2e-8,num_replicates=22)

#54545

    
j=0
for i in dd:
    wow=open('mynewvcf{}.vcf'.format(j),'w')
    i.write_vcf(wow,2,str(j))
    wow.close()
#    print(dd.samples())
    population_labels= ["protogreek"]*50 + ["corinth"]*50 + ["amvrakia_now"]*50 + ["amvrakia_initial"]*50
    d=0
    newlabels=[]
    for i in range(0,len(population_labels)):
        newlabels.append(population_labels[i]+str(d))
        d+=1
        if i==len(population_labels)-2:
            newlabels.append(population_labels[i]+str(d))
            break
        if population_labels[i]!=population_labels[i+1]:
            d=0
    population_labels=newlabels
    wow=open('mynewvcf{}.vcf'.format(j))
    wowzers=open('myvcf{}.vcf'.format(j),'w')
    for line in wow:
        line=line.strip().split()
        if line[0]=='#CHROM':
            line[9:]=population_labels
        wowzers.write("\t".join(line))
        wowzers.write("\n")
    wow.close()
    wowzers.close()
    j+=1
