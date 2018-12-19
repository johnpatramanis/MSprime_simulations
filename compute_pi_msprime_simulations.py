#inf = open('/Users/Eaaswarkhanth/Desktop/double_introgression/MSprime_simulations_master/Introgressed_variants_1.gen', 'r')

def get_pi_introg_introg(introg_der_allele_count, len_introg_allele_l): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	#introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	#len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
	#print('len_introg_allele_l: ', len_introg_allele_l)
	if int(len_introg_allele_l) >= 2:
		pi = (2.0 * float(introg_der_allele_count) * (len_introg_allele_l - introg_der_allele_count) )/ (len_introg_allele_l * ( len_introg_allele_l - 1 ))
	else:
		pi = 'NA'
	#print ('introg_der_allele_count: ', introg_der_allele_count, 'len_introg_allele_l: ', len_introg_allele_l, 'pi: ', pi )
	return pi




def get_pi_hominin_nonintrog(nonintrog_der_allele_count, nonintrog_der_allele_freq, len_nonintrog_allele_l, hominin_der_allele_count, hominin_der_allele_freq): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	#hominin_der_allele_count = get_hominin_der_al_count_freq(hominin, line2)[0]
	#hominin_der_allele_freq = get_hominin_der_al_count_freq(hominin, line2)[1]
	#nonintrog_der_allele_count = get_non_introg_der_al_count_freq(population, line, line2)[0]
	#nonintrog_der_allele_freq = get_non_introg_der_al_count_freq(population, line, line2)[2]
	#len_nonintrog_allele_l = get_non_introg_der_al_count_freq(population, line, line2)[1]
	if hominin_der_allele_freq == 0:
		pi = float(nonintrog_der_allele_freq)
	elif hominin_der_allele_freq == 1:
		pi = 1.0 - float(nonintrog_der_allele_freq)
	if nonintrog_der_allele_freq == 0:
		pi = float(hominin_der_allele_freq)
	elif nonintrog_der_allele_freq == 1:
		pi = 1.0 - float(hominin_der_allele_freq)
	if (nonintrog_der_allele_freq == 0 and hominin_der_allele_freq == 0) or (nonintrog_der_allele_freq == 1 and hominin_der_allele_freq == 1):
		pi = 0.0
	if nonintrog_der_allele_freq != 0 and hominin_der_allele_freq != 0 and nonintrog_der_allele_freq != 1 and hominin_der_allele_freq != 1:
		pi =  ( (float(hominin_der_allele_count) * (len_nonintrog_allele_l - nonintrog_der_allele_count)) + (float(nonintrog_der_allele_count) * (6.0 - hominin_der_allele_count))   )/ ( len_nonintrog_allele_l * 6.0 )
	#print ('hominin_der_allele_freq: ', hominin_der_allele_freq, 'nonintrog_der_allele_freq: ', nonintrog_der_allele_freq, 'pi: ', pi )
	return 	pi


	
def get_pi_introg_nonintrog(introg_der_allele_count, introg_der_allele_freq, len_introg_allele_l, nonintrog_der_allele_count, nonintrog_der_allele_freq, len_non_introg_allele_l): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	#nonintrog_der_allele_count = get_non_introg_der_al_count_freq(population, line, line2)[0]
	#nonintrog_der_allele_freq = get_non_introg_der_al_count_freq(population, line, line2)[2]
	#len_non_introg_allele_l  = get_non_introg_der_al_count_freq(population, line, line2)[1]
	#introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	#introg_der_allele_freq = get_introg_der_al_count_freq(line, line2)[2]
	#len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
	if introg_der_allele_freq == 0:
		pi = float(nonintrog_der_allele_freq)
	elif introg_der_allele_freq == 1:
		pi = 1.0 - float(nonintrog_der_allele_freq)
	if nonintrog_der_allele_freq == 0:
		pi = float(introg_der_allele_freq)
	elif nonintrog_der_allele_freq == 1:
		pi = 1.0 - float(introg_der_allele_freq)
	if (nonintrog_der_allele_freq == 0 and introg_der_allele_freq == 0) or (nonintrog_der_allele_freq == 1 and introg_der_allele_freq == 1):
		pi = 0.0
	if introg_der_allele_freq != 0 and nonintrog_der_allele_freq != 0 and introg_der_allele_freq != 1 and nonintrog_der_allele_freq != 1:
		pi = (  (float(nonintrog_der_allele_count) * (len_introg_allele_l - introg_der_allele_count)) + (float(introg_der_allele_count) * (len_non_introg_allele_l - nonintrog_der_allele_count)) )/(len_non_introg_allele_l * len_introg_allele_l)
	#print ('introg_der_allele_freq: ', introg_der_allele_freq, 'nonintrog_der_allele_freq: ', nonintrog_der_allele_freq, 'pi: ', pi)
	return pi


	
def get_pi_hominin_hominin(hominin1_der_allele_count, hominin1_der_allele_freq, hominin2_der_allele_count, hominin2_der_allele_freq): #line from the genotypes input file (inf2)	
	#hominin1_der_allele_freq = get_hominin_der_al_count_freq(hominin1, line)[1]
	#hominin2_der_allele_freq = get_hominin_der_al_count_freq(hominin2, line)[1]
	if hominin1_der_allele_freq == 0:
		pi = float(hominin2_der_allele_freq)
	if hominin2_der_allele_freq == 0:
		pi = float(hominin1_der_allele_freq)
	if hominin1_der_allele_freq == 1 and hominin2_der_allele_freq != 1: 
		pi = 1.0 - float(hominin2_der_allele_freq)
	elif hominin2_der_allele_freq == 1 and hominin1_der_allele_freq != 1:
		pi = 1.0 - float(hominin1_der_allele_freq)
	elif hominin1_der_allele_freq == 1 and hominin2_der_allele_freq == 1:
		pi = 0.0
	elif hominin1_der_allele_freq == 0 and hominin2_der_allele_freq == 0:
		pi = 0.0
	if hominin1_der_allele_freq != 1 and hominin1_der_allele_freq != 0 and hominin2_der_allele_freq != 0 and hominin2_der_allele_freq != 1:
		pi =  (    (float(hominin1_der_allele_count) * (6.0 - hominin2_der_allele_count))   + (float(hominin2_der_allele_count) * (6.0 - hominin1_der_allele_count))    ) / ( 6.0 * 6.0 )
	#print('hominin1_der_allele_freq: ', hominin1_der_allele_freq, 'hominin2_der_allele_freq: ', hominin2_der_allele_freq, 'pi: ', pi)
	return pi






def get_pi_introg_hominin(introg_der_allele_count, introg_der_allele_freq, len_introg_allele_l, hominin_der_allele_count, hominin_der_allele_freq): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	#hominin_der_allele_count = get_hominin_der_al_count_freq(hominin, line2)[0]
	#hominin_der_allele_freq = get_hominin_der_al_count_freq(hominin, line2)[1]
	#introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	#introg_der_allele_freq = get_introg_der_al_count_freq(line, line2)[2]
	#len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
	if ((introg_der_allele_freq != 0) and (hominin_der_allele_freq != 0) and (introg_der_allele_freq != 1) and (hominin_der_allele_freq != 1)):		
		pi =  (    (float(hominin_der_allele_count) * (len_introg_allele_l - introg_der_allele_count))   + (float(introg_der_allele_count) * (6.0 - hominin_der_allele_count))    ) / ( len_introg_allele_l * 6.0 )
	if hominin_der_allele_freq == 0:
		pi = float(introg_der_allele_freq)
	elif hominin_der_allele_freq == 1:
		pi = 1.0 - float(introg_der_allele_freq)
	if introg_der_allele_freq == 0:
		pi = float(hominin_der_allele_freq)
	elif introg_der_allele_freq == 1:
		pi = 1.0 - float(hominin_der_allele_freq)
	if ((introg_der_allele_freq == 0 and hominin_der_allele_freq == 0) or (introg_der_allele_freq == 1 and hominin_der_allele_freq == 1)):
		pi = 0.0   
	#print ('hominin_der_allele_freq: ', hominin_der_allele_freq, 'introg_der_allele_freq: ', introg_der_allele_freq, 'pi: ', pi )
	return pi


total_no_variants = 0

count_IndexError =0

count = 0
total_no_sites_pi_50kb = 0
total_pi_introg_N1 = 0.0
total_pi_introg_N2 = 0.0
total_pi_N1_N2 = 0.0
total_pi_introg_introg = 0.0

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

inf = open(input_file, 'r')
outf = open(output_file, 'w')

check = 0
window_size = 50000
for line in inf:
	if line.startswith('#') or line == '\n':
		check = 1
		continue		
	count += 1
	total_no_variants += 1
	l = line.split()
	try: 
		pos = float(l[1])
	except IndexError:
		count_IndexError += 1
		print(line)
		continue
	#print int(pos)
	#print int(first_pos)
	#if line.startswith('#') or line == '\n':
	#	continue
	if count == 1 or check == 1:
		first_pos = pos
		check = 0 
	if (int(pos)>= int(first_pos)) and (int(pos) < (int(first_pos) + window_size)): 
		genos = l[0]
		inds = l[2].split(',')
		#print (pos)
		#print (len(geno))
		neand1_geno = genos[6:12]
		neand2_geno = genos[12:18]
		afr_geno = genos[18:218]
		eur_geno = genos[218:768]
		easn_geno = genos[768:1318]
		#print ('neand1_geno: ', neand1_geno)
		#print ('neand2_geno: ', neand2_geno)
		afr_geno_no_derived = afr_geno.count('1')
		eur_geno_no_derived = eur_geno.count('1')
		easn_geno_no_derived = easn_geno.count('1')
		neand1_geno_no_derived = neand1_geno.count('1')
		neand2_geno_no_derived = neand2_geno.count('1')
		introg_geno_l = []
		nonintrog_geno_l = []
		for ind in inds:
			introg_geno_l.append(genos[int(ind)])
		for i,e in enumerate(genos):
			if str(i) not in inds:
				nonintrog_geno_l.append(e)
		introg_der_allele_c = introg_geno_l.count('1')
		size_introg_allele_l = len(introg_geno_l)
		introg_der_allele_f = float(introg_der_allele_c)/size_introg_allele_l
	
		nonintrog_der_allele_c = nonintrog_geno_l.count('1')
		size_nonintrog_allele_l = len(nonintrog_geno_l)
		nonintrog_der_allele_f = float(nonintrog_der_allele_c)/size_nonintrog_allele_l
	
		N1_der_allele_f = neand1_geno_no_derived/6.0
		N2_der_allele_f = neand2_geno_no_derived/6.0
	
		pi_introg_N1 = get_pi_introg_hominin(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, neand1_geno_no_derived, N1_der_allele_f)
		pi_introg_N2 = get_pi_introg_hominin(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, neand2_geno_no_derived, N2_der_allele_f)
		#pi_nonintrog_N1 = get_pi_hominin_nonintrog(nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l, neand1_geno_no_derived, N1_der_allele_f)
		#pi_nonintrog_N2 = get_pi_hominin_nonintrog(nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l, neand2_geno_no_derived, N2_der_allele_f)
		pi_N1_N2 = get_pi_hominin_hominin(neand1_geno_no_derived, N1_der_allele_f, neand2_geno_no_derived, N2_der_allele_f) #line from the genotypes input file (inf2)	
		#pi_introg_nonintrog = get_pi_introg_nonintrog(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l)
		pi_introg_introg = get_pi_introg_introg(introg_der_allele_c, size_introg_allele_l)	
		total_pi_introg_N1 += pi_introg_N1
		total_pi_introg_N2 += pi_introg_N2
		total_pi_N1_N2 += pi_N1_N2
		if pi_introg_introg != 'NA':
			total_pi_introg_introg += pi_introg_introg
	elif (int(pos) >= (int(first_pos) + window_size)):
		first_pos = pos 
		#print str(total_pi_introg_N1/window_size)
		pi_introg_neand1 = [ str(total_pi_introg_N1/float(window_size)) ]
		pi_introg_neand2 = [ str(total_pi_introg_N2/float(window_size)) ]
		pi_neand1_neand2 = [ str(total_pi_N1_N2/float(window_size)) ]
		pi_introg_introgresed = [ str(total_pi_introg_introg/float(window_size)) ] 
		
		total_no_sites_pi_50kb = 0
		total_pi_introg_N1 = 0.0
		total_pi_introg_N2 = 0.0
		total_pi_N1_N2 = 0.0
		total_pi_introg_introg = 0.0
		#print ('\t'.join(pi_introg_neand1, pi_introg_neand2, pi_neand1_neand2, pi_introg_introgresed), file = outf  )	
		print >>outf, '\t'.join(pi_introg_neand1 + pi_introg_neand2 + pi_neand1_neand2 + pi_introg_introgresed)	
		genos = l[0]
		inds = l[2].split(',')
		neand1_geno = genos[6:12]
		neand2_geno = genos[12:18]
		afr_geno = genos[18:218]
		eur_geno = genos[218:768]
		easn_geno = genos[768:1318]
		afr_geno_no_derived = afr_geno.count('1')
		eur_geno_no_derived = eur_geno.count('1')
		easn_geno_no_derived = easn_geno.count('1')
		neand1_geno_no_derived = neand1_geno.count('1')
		neand2_geno_no_derived = neand2_geno.count('1')
		introg_geno_l = []
		nonintrog_geno_l = []
		for ind in inds:
			introg_geno_l.append(genos[int(ind)])
		for i,e in enumerate(genos):
			if str(i) not in inds:
				nonintrog_geno_l.append(e)
		introg_der_allele_c = introg_geno_l.count('1')
		size_introg_allele_l = len(introg_geno_l)
		introg_der_allele_f = float(introg_der_allele_c)/size_introg_allele_l
	
		nonintrog_der_allele_c = nonintrog_geno_l.count('1')
		size_nonintrog_allele_l = len(nonintrog_geno_l)
		nonintrog_der_allele_f = float(nonintrog_der_allele_c)/size_nonintrog_allele_l
	
		N1_der_allele_f = neand1_geno_no_derived/6.0
		N2_der_allele_f = neand2_geno_no_derived/6.0
	
		pi_introg_N1 = get_pi_introg_hominin(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, neand1_geno_no_derived, N1_der_allele_f)
		pi_introg_N2 = get_pi_introg_hominin(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, neand2_geno_no_derived, N2_der_allele_f)
		#pi_nonintrog_N1 = get_pi_hominin_nonintrog(nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l, neand1_geno_no_derived, N1_der_allele_f)
		#pi_nonintrog_N2 = get_pi_hominin_nonintrog(nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l, neand2_geno_no_derived, N2_der_allele_f)
		pi_N1_N2 = get_pi_hominin_hominin(neand1_geno_no_derived, N1_der_allele_f, neand2_geno_no_derived, N2_der_allele_f) #line from the genotypes input file (inf2)	
		#pi_introg_nonintrog = get_pi_introg_nonintrog(introg_der_allele_c, introg_der_allele_f, size_introg_allele_l, nonintrog_der_allele_c, nonintrog_der_allele_f, size_nonintrog_allele_l)
		pi_introg_introg = get_pi_introg_introg(introg_der_allele_c, size_introg_allele_l)	
		total_pi_introg_N1 += pi_introg_N1
		total_pi_introg_N2 += pi_introg_N2
		total_pi_N1_N2 += pi_N1_N2
		if pi_introg_introg != 'NA':
			total_pi_introg_introg += pi_introg_introg

		
	
	

	 

inf.close()
outf.close()
		



