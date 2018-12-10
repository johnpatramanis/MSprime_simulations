
import sys
 
# file_name = sys.argv[1]
# out_name = sys.argv[2]
# out_name2 = sys.argv[3]
# out_name3 = sys.argv[4]
# 
# chr_no = sys.argv[5] 

#file_name = 'Neand_haplos_inds_newIDs_fineThres_EASN_phased_run1_10.sorted.dupl2.chr.21.bed'
#out_name =  'Neand_haplos_inds_newIDs_fineThres_EASN_phased_run1_10.sorted.dupl2.pi.snpfreq.chr.21.bed'
#out_name2 =  'Neand_haplos_inds_newIDs_fineThres_EASN_phased_run1_10.sorted.dupl2.pi.chr.21.bed'

# inf = open(file_name, 'r')
# outf = open(out_name, 'w')
# outf2 = open(out_name2, 'w')
# outf3 = open(out_name3, 'w')


inf = open('/projects/academic/omergokc/ozgur_Sstats/double_introgression/new_Sstar_runs_easn/Neand_haplos_inds_newIDs_fineThres_EASN_phased_run1_10.sorted.dupl2.chr.21.bed', 'r')
#outf = open('/Users/Eaaswarkhanth/Desktop/double_introgression/haplotypes_Sstar/a.bed', 'w')
outf2 = open('/projects/academic/omergokc/ozgur_Sstats/double_introgression/new_Sstar_runs_easn/output_a_oocode_chr21.bed', 'w')
outf3 = open('/projects/academic/omergokc/ozgur_Sstats/double_introgression/new_Sstar_runs_easn/output_b_oocode_chr21.bed', 'w')


	
	
	
#Get individuals, prop of derived alleles for Sstar significant SNPs per phased haplotype, number of introgressed haplotypes 	
def get_inds_percDerived_noHaplos(line): #line from the phased haplotypes input file (inf1)
	print (line[-1].split('\n')[0])
	if line[-1].split('\n')[0] == '7':
		individuals = [ str(line[4])]
		percent_derived_avr = [ str(line[5])]
		noHaplos = ['1']		
	elif line[-1].split('\n')[0] == '9':
		individuals = [ line[4] + '_' + line[5] ] 
		percent_derived_avr = [ str( ( float(line[6]) + float(line[7]) )/2.0 ) ]
		noHaplos = ['2']		
	elif int(line[-1].split('\n')[0]) > 9:
		ind_and_hapl_l = []
		percent_shared_l = []
		no_ind_perc = len(line[4:-2])
		#print (no_ind_perc)
		stop_index =  4 + (int(no_ind_perc)/2)
		#print(stop_index)
		for x in line[4 : int(stop_index)]:
			ind_and_hapl_l.append(x)
		for x in line[int(stop_index):-2]:
			percent_shared_l.append(float(x))				
		individuals = '_'.join(ind_and_hapl_l) 
		percent_derived_avr = [ str( sum(percent_shared_l) / float(len(percent_shared_l)) )]
		noHaplos = [str( len(ind_and_hapl_l) ) ]
		#print percent_derived_avg
	return individuals, percent_derived_avr, noHaplos



##Get Sstar significant SNP list
def get_Sstar_snps(line): #line from the phased haplotypes input file (inf1)
	snps = line[-2]
	snps = snps.split(', ')
	SNP1 = snps[0].split('[')[1]
	SNP_LAST = snps[-1].split(']')[0]
	count2 = 0
	new_snps = []
	for SNP in snps:
		count2 += 1
		if count2 == 1 :
			new_snps.append(SNP1[1:-1])
		elif count2 > 1 and count2 < len(snps):
			new_snps.append(SNP[1:-1])
		elif count2 == len(snps):
			new_snps.append(SNP_LAST[1:-1])	
	newSNPs = new_snps
	len_haplo = int(SNP_LAST[1:-1]) - int(SNP1[1:-1])
	return (newSNPs, len_haplo)
	

			

##Get introgressed haplotype info
def get_ind_col_id_list(line): #line from the phased haplotypes input file (inf1) 
	kg1_header = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00104', 'HG00106', 'HG00108', 'HG00109', 'HG00110', 'HG00111', 'HG00112', 'HG00113', 'HG00114', 'HG00116', 'HG00117', 'HG00118', 'HG00119', 'HG00120', 'HG00121', 'HG00122', 'HG00123', 'HG00124', 'HG00125', 'HG00126', 'HG00127', 'HG00128', 'HG00129', 'HG00130', 'HG00131', 'HG00133', 'HG00134', 'HG00135', 'HG00136', 'HG00137', 'HG00138', 'HG00139', 'HG00140', 'HG00141', 'HG00142', 'HG00143', 'HG00146', 'HG00148', 'HG00149', 'HG00150', 'HG00151', 'HG00152', 'HG00154', 'HG00155', 'HG00156', 'HG00158', 'HG00159', 'HG00160', 'HG00171', 'HG00173', 'HG00174', 'HG00176', 'HG00177', 'HG00178', 'HG00179', 'HG00180', 'HG00182', 'HG00183', 'HG00185', 'HG00186', 'HG00187', 'HG00188', 'HG00189', 'HG00190', 'HG00231', 'HG00232', 'HG00233', 'HG00234', 'HG00235', 'HG00236', 'HG00237', 'HG00238', 'HG00239', 'HG00240', 'HG00242', 'HG00243', 'HG00244', 'HG00245', 'HG00246', 'HG00247', 'HG00249', 'HG00250', 'HG00251', 'HG00252', 'HG00253', 'HG00254', 'HG00255', 'HG00256', 'HG00257', 'HG00258', 'HG00259', 'HG00260', 'HG00261', 'HG00262', 'HG00263', 'HG00264', 'HG00265', 'HG00266', 'HG00267', 'HG00268', 'HG00269', 'HG00270', 'HG00271', 'HG00272', 'HG00273', 'HG00274', 'HG00275', 'HG00276', 'HG00277', 'HG00278', 'HG00280', 'HG00281', 'HG00282', 'HG00284', 'HG00285', 'HG00306', 'HG00309', 'HG00310', 'HG00311', 'HG00312', 'HG00313', 'HG00315', 'HG00318', 'HG00319', 'HG00320', 'HG00321', 'HG00323', 'HG00324', 'HG00325', 'HG00326', 'HG00327', 'HG00328', 'HG00329', 'HG00330', 'HG00331', 'HG00332', 'HG00334', 'HG00335', 'HG00336', 'HG00337', 'HG00338', 'HG00339', 'HG00341', 'HG00342', 'HG00343', 'HG00344', 'HG00345', 'HG00346', 'HG00349', 'HG00350', 'HG00351', 'HG00353', 'HG00355', 'HG00356', 'HG00357', 'HG00358', 'HG00359', 'HG00360', 'HG00361', 'HG00362', 'HG00364', 'HG00366', 'HG00367', 'HG00369', 'HG00372', 'HG00373', 'HG00375', 'HG00376', 'HG00377', 'HG00378', 'HG00381', 'HG00382', 'HG00383', 'HG00384', 'HG00403', 'HG00404', 'HG00406', 'HG00407', 'HG00418', 'HG00419', 'HG00421', 'HG00422', 'HG00427', 'HG00428', 'HG00436', 'HG00437', 'HG00442', 'HG00443', 'HG00445', 'HG00446', 'HG00448', 'HG00449', 'HG00451', 'HG00452', 'HG00457', 'HG00458', 'HG00463', 'HG00464', 'HG00472', 'HG00473', 'HG00475', 'HG00476', 'HG00478', 'HG00479', 'HG00500', 'HG00501', 'HG00512', 'HG00513', 'HG00524', 'HG00525', 'HG00530', 'HG00531', 'HG00533', 'HG00534', 'HG00536', 'HG00537', 'HG00542', 'HG00543', 'HG00553', 'HG00554', 'HG00556', 'HG00557', 'HG00559', 'HG00560', 'HG00565', 'HG00566', 'HG00577', 'HG00578', 'HG00580', 'HG00581', 'HG00583', 'HG00584', 'HG00589', 'HG00590', 'HG00592', 'HG00593', 'HG00595', 'HG00596', 'HG00607', 'HG00608', 'HG00610', 'HG00611', 'HG00613', 'HG00614', 'HG00619', 'HG00620', 'HG00625', 'HG00626', 'HG00628', 'HG00629', 'HG00634', 'HG00635', 'HG00637', 'HG00638', 'HG00640', 'HG00641', 'HG00650', 'HG00651', 'HG00653', 'HG00654', 'HG00656', 'HG00657', 'HG00662', 'HG00663', 'HG00671', 'HG00672', 'HG00683', 'HG00684', 'HG00689', 'HG00690', 'HG00692', 'HG00693', 'HG00698', 'HG00699', 'HG00701', 'HG00702', 'HG00704', 'HG00705', 'HG00707', 'HG00708', 'HG00731', 'HG00732', 'HG00734', 'HG00736', 'HG00737', 'HG00740', 'HG01047', 'HG01048', 'HG01051', 'HG01052', 'HG01055', 'HG01060', 'HG01061', 'HG01066', 'HG01067', 'HG01069', 'HG01070', 'HG01072', 'HG01073', 'HG01075', 'HG01079', 'HG01080', 'HG01082', 'HG01083', 'HG01085', 'HG01095', 'HG01097', 'HG01098', 'HG01101', 'HG01102', 'HG01104', 'HG01105', 'HG01107', 'HG01108', 'HG01112', 'HG01113', 'HG01124', 'HG01125', 'HG01133', 'HG01134', 'HG01136', 'HG01137', 'HG01140', 'HG01148', 'HG01149', 'HG01167', 'HG01168', 'HG01170', 'HG01171', 'HG01173', 'HG01174', 'HG01176', 'HG01183', 'HG01187', 'HG01188', 'HG01190', 'HG01191', 'HG01197', 'HG01198', 'HG01204', 'HG01250', 'HG01251', 'HG01257', 'HG01259', 'HG01271', 'HG01272', 'HG01274', 'HG01275', 'HG01277', 'HG01278', 'HG01334', 'HG01342', 'HG01344', 'HG01345', 'HG01350', 'HG01351', 'HG01353', 'HG01354', 'HG01356', 'HG01357', 'HG01359', 'HG01360', 'HG01365', 'HG01366', 'HG01374', 'HG01375', 'HG01377', 'HG01378', 'HG01383', 'HG01384', 'HG01389', 'HG01390', 'HG01437', 'HG01440', 'HG01441', 'HG01455', 'HG01456', 'HG01461', 'HG01462', 'HG01465', 'HG01488', 'HG01489', 'HG01491', 'HG01492', 'HG01494', 'HG01495', 'HG01497', 'HG01498', 'HG01515', 'HG01516', 'HG01518', 'HG01519', 'HG01521', 'HG01522', 'HG01550', 'HG01551', 'HG01617', 'HG01618', 'HG01619', 'HG01620', 'HG01623', 'HG01624', 'HG01625', 'HG01626', 'NA06984', 'NA06986', 'NA06989', 'NA06994', 'NA07000', 'NA07037', 'NA07048', 'NA07051', 'NA07056', 'NA07347', 'NA07357', 'NA10847', 'NA10851', 'NA11829', 'NA11830', 'NA11831', 'NA11843', 'NA11892', 'NA11893', 'NA11894', 'NA11919', 'NA11920', 'NA11930', 'NA11931', 'NA11932', 'NA11933', 'NA11992', 'NA11993', 'NA11994', 'NA11995', 'NA12003', 'NA12004', 'NA12006', 'NA12043', 'NA12044', 'NA12045', 'NA12046', 'NA12058', 'NA12144', 'NA12154', 'NA12155', 'NA12249', 'NA12272', 'NA12273', 'NA12275', 'NA12282', 'NA12283', 'NA12286', 'NA12287', 'NA12340', 'NA12341', 'NA12342', 'NA12347', 'NA12348', 'NA12383', 'NA12399', 'NA12400', 'NA12413', 'NA12489', 'NA12546', 'NA12716', 'NA12717', 'NA12718', 'NA12748', 'NA12749', 'NA12750', 'NA12751', 'NA12761', 'NA12763', 'NA12775', 'NA12777', 'NA12778', 'NA12812', 'NA12814', 'NA12815', 'NA12827', 'NA12829', 'NA12830', 'NA12842', 'NA12843', 'NA12872', 'NA12873', 'NA12874', 'NA12889', 'NA12890', 'NA18486', 'NA18487', 'NA18489', 'NA18498', 'NA18499', 'NA18501', 'NA18502', 'NA18504', 'NA18505', 'NA18507', 'NA18508', 'NA18510', 'NA18511', 'NA18516', 'NA18517', 'NA18519', 'NA18520', 'NA18522', 'NA18523', 'NA18525', 'NA18526', 'NA18527', 'NA18528', 'NA18530', 'NA18532', 'NA18534', 'NA18535', 'NA18536', 'NA18537', 'NA18538', 'NA18539', 'NA18541', 'NA18542', 'NA18543', 'NA18544', 'NA18545', 'NA18546', 'NA18547', 'NA18548', 'NA18549', 'NA18550', 'NA18552', 'NA18553', 'NA18555', 'NA18557', 'NA18558', 'NA18559', 'NA18560', 'NA18561', 'NA18562', 'NA18563', 'NA18564', 'NA18565', 'NA18566', 'NA18567', 'NA18570', 'NA18571', 'NA18572', 'NA18573', 'NA18574', 'NA18576', 'NA18577', 'NA18579', 'NA18582', 'NA18592', 'NA18593', 'NA18595', 'NA18596', 'NA18597', 'NA18599', 'NA18602', 'NA18603', 'NA18605', 'NA18606', 'NA18608', 'NA18609', 'NA18610', 'NA18611', 'NA18612', 'NA18613', 'NA18614', 'NA18615', 'NA18616', 'NA18617', 'NA18618', 'NA18619', 'NA18620', 'NA18621', 'NA18622', 'NA18623', 'NA18624', 'NA18626', 'NA18627', 'NA18628', 'NA18630', 'NA18631', 'NA18632', 'NA18633', 'NA18634', 'NA18635', 'NA18636', 'NA18637', 'NA18638', 'NA18639', 'NA18640', 'NA18641', 'NA18642', 'NA18643', 'NA18645', 'NA18647', 'NA18740', 'NA18745', 'NA18747', 'NA18748', 'NA18749', 'NA18757', 'NA18853', 'NA18856', 'NA18858', 'NA18861', 'NA18867', 'NA18868', 'NA18870', 'NA18871', 'NA18873', 'NA18874', 'NA18907', 'NA18908', 'NA18909', 'NA18910', 'NA18912', 'NA18916', 'NA18917', 'NA18923', 'NA18924', 'NA18933', 'NA18934', 'NA18939', 'NA18940', 'NA18941', 'NA18942', 'NA18943', 'NA18944', 'NA18945', 'NA18946', 'NA18947', 'NA18948', 'NA18949', 'NA18950', 'NA18951', 'NA18952', 'NA18953', 'NA18954', 'NA18956', 'NA18957', 'NA18959', 'NA18960', 'NA18961', 'NA18962', 'NA18963', 'NA18964', 'NA18965', 'NA18966', 'NA18968', 'NA18971', 'NA18973', 'NA18974', 'NA18975', 'NA18976', 'NA18977', 'NA18978', 'NA18980', 'NA18981', 'NA18982', 'NA18983', 'NA18984', 'NA18985', 'NA18986', 'NA18987', 'NA18988', 'NA18989', 'NA18990', 'NA18992', 'NA18994', 'NA18995', 'NA18998', 'NA18999', 'NA19000', 'NA19002', 'NA19003', 'NA19004', 'NA19005', 'NA19007', 'NA19009', 'NA19010', 'NA19012', 'NA19020', 'NA19028', 'NA19035', 'NA19036', 'NA19038', 'NA19041', 'NA19044', 'NA19046', 'NA19054', 'NA19055', 'NA19056', 'NA19057', 'NA19058', 'NA19059', 'NA19060', 'NA19062', 'NA19063', 'NA19064', 'NA19065', 'NA19066', 'NA19067', 'NA19068', 'NA19070', 'NA19072', 'NA19074', 'NA19075', 'NA19076', 'NA19077', 'NA19078', 'NA19079', 'NA19080', 'NA19081', 'NA19082', 'NA19083', 'NA19084', 'NA19085', 'NA19087', 'NA19088', 'NA19093', 'NA19095', 'NA19096', 'NA19098', 'NA19099', 'NA19102', 'NA19107', 'NA19108', 'NA19113', 'NA19114', 'NA19116', 'NA19117', 'NA19118', 'NA19119', 'NA19121', 'NA19129', 'NA19130', 'NA19131', 'NA19137', 'NA19138', 'NA19146', 'NA19147', 'NA19149', 'NA19150', 'NA19152', 'NA19160', 'NA19171', 'NA19172', 'NA19175', 'NA19185', 'NA19189', 'NA19190', 'NA19197', 'NA19198', 'NA19200', 'NA19204', 'NA19207', 'NA19209', 'NA19213', 'NA19222', 'NA19223', 'NA19225', 'NA19235', 'NA19236', 'NA19247', 'NA19248', 'NA19256', 'NA19257', 'NA19307', 'NA19308', 'NA19309', 'NA19310', 'NA19311', 'NA19312', 'NA19313', 'NA19315', 'NA19316', 'NA19317', 'NA19318', 'NA19319', 'NA19321', 'NA19324', 'NA19327', 'NA19328', 'NA19331', 'NA19332', 'NA19334', 'NA19338', 'NA19346', 'NA19347', 'NA19350', 'NA19351', 'NA19352', 'NA19355', 'NA19359', 'NA19360', 'NA19371', 'NA19372', 'NA19373', 'NA19374', 'NA19375', 'NA19376', 'NA19377', 'NA19379', 'NA19380', 'NA19381', 'NA19382', 'NA19383', 'NA19384', 'NA19385', 'NA19390', 'NA19391', 'NA19393', 'NA19394', 'NA19395', 'NA19396', 'NA19397', 'NA19398', 'NA19399', 'NA19401', 'NA19403', 'NA19404', 'NA19428', 'NA19429', 'NA19430', 'NA19431', 'NA19434', 'NA19435', 'NA19436', 'NA19437', 'NA19438', 'NA19439', 'NA19440', 'NA19443', 'NA19444', 'NA19445', 'NA19446', 'NA19448', 'NA19449', 'NA19451', 'NA19452', 'NA19453', 'NA19455', 'NA19456', 'NA19457', 'NA19461', 'NA19462', 'NA19463', 'NA19466', 'NA19467', 'NA19468', 'NA19469', 'NA19470', 'NA19471', 'NA19472', 'NA19473', 'NA19474', 'NA19625', 'NA19648', 'NA19651', 'NA19652', 'NA19654', 'NA19655', 'NA19657', 'NA19660', 'NA19661', 'NA19663', 'NA19664', 'NA19672', 'NA19675', 'NA19676', 'NA19678', 'NA19679', 'NA19681', 'NA19682', 'NA19684', 'NA19685', 'NA19700', 'NA19701', 'NA19703', 'NA19704', 'NA19707', 'NA19711', 'NA19712', 'NA19713', 'NA19716', 'NA19717', 'NA19719', 'NA19720', 'NA19722', 'NA19723', 'NA19725', 'NA19726', 'NA19728', 'NA19729', 'NA19731', 'NA19732', 'NA19734', 'NA19735', 'NA19737', 'NA19738', 'NA19740', 'NA19741', 'NA19746', 'NA19747', 'NA19749', 'NA19750', 'NA19752', 'NA19753', 'NA19755', 'NA19756', 'NA19758', 'NA19759', 'NA19761', 'NA19762', 'NA19764', 'NA19770', 'NA19771', 'NA19773', 'NA19774', 'NA19776', 'NA19777', 'NA19779', 'NA19780', 'NA19782', 'NA19783', 'NA19785', 'NA19786', 'NA19788', 'NA19789', 'NA19794', 'NA19795', 'NA19818', 'NA19819', 'NA19834', 'NA19835', 'NA19900', 'NA19901', 'NA19904', 'NA19908', 'NA19909', 'NA19914', 'NA19916', 'NA19917', 'NA19920', 'NA19921', 'NA19922', 'NA19923', 'NA19982', 'NA19984', 'NA19985', 'NA20126', 'NA20127', 'NA20276', 'NA20278', 'NA20281', 'NA20282', 'NA20287', 'NA20289', 'NA20291', 'NA20294', 'NA20296', 'NA20298', 'NA20299', 'NA20314', 'NA20317', 'NA20322', 'NA20332', 'NA20334', 'NA20336', 'NA20339', 'NA20340', 'NA20341', 'NA20342', 'NA20344', 'NA20346', 'NA20348', 'NA20351', 'NA20356', 'NA20357', 'NA20359', 'NA20363', 'NA20412', 'NA20414', 'NA20502', 'NA20503', 'NA20504', 'NA20505', 'NA20506', 'NA20507', 'NA20508', 'NA20509', 'NA20510', 'NA20512', 'NA20513', 'NA20515', 'NA20516', 'NA20517', 'NA20518', 'NA20519', 'NA20520', 'NA20521', 'NA20522', 'NA20524', 'NA20525', 'NA20527', 'NA20528', 'NA20529', 'NA20530', 'NA20531', 'NA20532', 'NA20533', 'NA20534', 'NA20535', 'NA20536', 'NA20537', 'NA20538', 'NA20539', 'NA20540', 'NA20541', 'NA20542', 'NA20543', 'NA20544', 'NA20581', 'NA20582', 'NA20585', 'NA20586', 'NA20588', 'NA20589', 'NA20752', 'NA20753', 'NA20754', 'NA20755', 'NA20756', 'NA20757', 'NA20758', 'NA20759', 'NA20760', 'NA20761', 'NA20765', 'NA20766', 'NA20768', 'NA20769', 'NA20770', 'NA20771', 'NA20772', 'NA20773', 'NA20774', 'NA20775', 'NA20778', 'NA20783', 'NA20785', 'NA20786', 'NA20787', 'NA20790', 'NA20792', 'NA20795', 'NA20796', 'NA20797', 'NA20798', 'NA20799', 'NA20800', 'NA20801', 'NA20802', 'NA20803', 'NA20804', 'NA20805', 'NA20806', 'NA20807', 'NA20808', 'NA20809', 'NA20810', 'NA20811', 'NA20812', 'NA20813', 'NA20814', 'NA20815', 'NA20816', 'NA20818', 'NA20819', 'NA20826', 'NA20828']
	ind_column_id_l = []
	if line[-1].split('\n')[0] == '7':
		ind = line[4].split('_')[0]
		for i,e in enumerate(kg1_header):
			if ind == e:
				ind_column_id = i + 5 
				ind_column_id_l.append(ind_column_id)
	elif line[-1].split('\n')[0] == '9':
		ind1 = line[4].split('_')[0]
		ind2 = line[5].split('_')[0]
		for i,e in enumerate(kg1_header):
			if ind1 == e:
				ind_column_id_1 = i + 5
				ind_column_id_l.append(ind_column_id_1)
			if ind2 == e:
				ind_column_id_2 = i + 5
				ind_column_id_l.append(ind_column_id_2)
	elif int(line[-1].split('\n')[0]) > 9:
		ind_l = []
		no_ind_perc = len(line[4:-2])
		for x in line[4: int(4 + (no_ind_perc/2))]:
			ind = x.split('_')[0]
			ind_l.append(ind)
		for k in ind_l:
			for i,e in enumerate(kg1_header):
				if e == k:
					ind_column_id = i + 5
					ind_column_id_l.append(ind_column_id)
	return ind_column_id_l



def get_ind_hapl_list(line): #line from the phased haplotypes input file (inf1) 
	ind_hapl_l = []
	if int(line[-1].split('\n')[0]) == 7:		
		haplo_no = line[4].split('_')[1][-1]
		ind_hapl_l.append(haplo_no)
	elif int(line[-1].split('\n')[0]) == 9:
		haplo_no1 = line[4].split('_')[1][-1]
		haplo_no2 = line[5].split('_')[1][-1]
		ind_hapl_l.append(haplo_no1)
		ind_hapl_l.append(haplo_no2)
	elif int(line[-1].split('\n')[0]) > 9:
		no_ind_perc = len(line[4:-2])
		for x in line[4: int(4 + (no_ind_perc/2))]:
			haplo_no = x.split('_')[1][-1]
			ind_hapl_l.append(haplo_no)
	return ind_hapl_l



#Get Derived Allele Frequencies
def get_eurasian_der_al_freq(population, line): #line from the genotypes input file (inf2)
	population = population.lower()
	if population == 'easn':
		index_list = [181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729]
		count = 0.0
		for i in index_list:
			genos = line[i+5].split('|')
			geno1 = genos[0]
			geno2 = genos[1]
			if geno1 == '1': 
				count += 1
			if geno2 == '1': 
				count += 1
		pi = float(count)/ (len(index_list) * 2.0)
	elif population == 'nweur':
		index_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 357, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495]
		count = 0.0
		for i in index_list:
			genos = line[i+5].split('|')
			geno1 = genos[0]
			geno2 = genos[1]
			if geno1 == '1': 
				count += 1
			if geno2 == '1': 
				count += 1
		pi = float(count)/ (len(index_list) * 2.0)
	elif population == 'yri':
		index_list = [497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778]
		count = 0.0
		for i in index_list:
			genos = line[i].split('|')
			geno1 = genos[0]
			geno2 = genos[1]
			if geno1 == '1': 
				count += 1
			if geno2 == '1': 
				count += 1
		pi = float(count)/ (len(index_list) * 2.0)
	return pi
	
	
def get_hominin_der_al_count_freq(hominin, line): #line from the genotypes input file (inf2)
	hominin = hominin.lower()
	if hominin == 'altai':
		index = line[-19]
	elif hominin == 'vindija':
		index = line[-13]
	elif hominin == 'denisova':
		index = line[-7]
	elif hominin == 'chagyrskaya':
		index = line[-1]
	genos = index.split('/')
	return (genos.count('1'),  genos.count('1')/2.0)



def get_introg_der_al_count_freq(line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	ind_column_id_l = get_ind_col_id_list(line)
	ind_hapl_l = get_ind_hapl_list(line)
	#print (ind_column_id_l)
	#print (ind_hapl_l)
	introg_allele_l = []
	for i,e in enumerate(ind_column_id_l):
		alleles = line2[e].split('|')
		introg_allele = alleles[int(ind_hapl_l[i])-1]
		introg_allele_l.append(introg_allele)
	#print (introg_allele_l)
	return (introg_allele_l.count('1'), len(introg_allele_l), float(introg_allele_l.count('1'))/len(introg_allele_l)) 


def get_non_introg_der_al_count_freq(population, line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	ind_column_id_l = get_ind_col_id_list(line)
	non_introg_allele_l = []
	population = population.lower()
	if population == 'easn':
		index_list = [181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729]
	elif population == 'nweur':
		index_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 357, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495]
	for j in index_list:
		if int(j+5) not in ind_column_id_l:
			non_introg_geno = line2[j+5].split('|')
			non_introg_allele1 = non_introg_geno[0]   
			non_introg_allele2 = non_introg_geno[1]   
			non_introg_allele_l.append(non_introg_allele1)	
			non_introg_allele_l.append(non_introg_allele2)
	return (non_introg_allele_l.count('1'), len(non_introg_allele_l), float(non_introg_allele_l.count('1'))/len(non_introg_allele_l))





##Get pi 
def get_pi_nonintrog_nonintrog(population, line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	non_introg_der_allele_count = get_non_introg_der_al_count_freq(population, line, line2)[0] 
	len_non_introg_allele_l = get_non_introg_der_al_count_freq(population, line, line2)[1]
	pi = (2.0 * float(non_introg_der_allele_count) * (len_non_introg_allele_l - non_introg_der_allele_count) )/ (len_non_introg_allele_l * ( len_non_introg_allele_l - 1 ))
	print ('non_introg_der_allele_count: ', non_introg_der_allele_count, 'len_non_introg_allele_l: ', len_non_introg_allele_l, 'pi: ', pi )
	return pi


def get_pi_introg_introg(line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
	print('len_introg_allele_l: ', len_introg_allele_l)
	if int(len_introg_allele_l) >= 2:
		pi = (2.0 * float(introg_der_allele_count) * (len_introg_allele_l - introg_der_allele_count) )/ (len_introg_allele_l * ( len_introg_allele_l - 1 ))
	else:
		pi = 'NA'
	print ('introg_der_allele_count: ', introg_der_allele_count, 'len_introg_allele_l: ', len_introg_allele_l, 'pi: ', pi )
	return pi


def get_pi_introg_hominin(hominin, line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	hominin_der_allele_count = get_hominin_der_al_count_freq(hominin, line2)[0]
	hominin_der_allele_freq = get_hominin_der_al_count_freq(hominin, line2)[1]
	introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	introg_der_allele_freq = get_introg_der_al_count_freq(line, line2)[2]
	len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
	if ((introg_der_allele_freq != 0) and (hominin_der_allele_freq != 0) and (introg_der_allele_freq != 1) and (hominin_der_allele_freq != 1)):		
		pi =  (    (float(hominin_der_allele_count) * (len_introg_allele_l - introg_der_allele_count))   + (float(introg_der_allele_count) * (2.0 - hominin_der_allele_count))    ) / ( len_introg_allele_l * 2.0 )
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
	print ('hominin_der_allele_freq: ', hominin_der_allele_freq, 'introg_der_allele_freq: ', introg_der_allele_freq, 'pi: ', pi )
	return pi



def get_pi_hominin_nonintrog(hominin, population, line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	hominin_der_allele_count = get_hominin_der_al_count_freq(hominin, line2)[0]
	hominin_der_allele_freq = get_hominin_der_al_count_freq(hominin, line2)[1]
	nonintrog_der_allele_count = get_non_introg_der_al_count_freq(population, line, line2)[0]
	nonintrog_der_allele_freq = get_non_introg_der_al_count_freq(population, line, line2)[2]
	len_nonintrog_allele_l = get_non_introg_der_al_count_freq(population, line, line2)[1]
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
		pi =  ( (float(hominin_der_allele_count) * (len_nonintrog_allele_l - nonintrog_der_allele_count)) + (float(nonintrog_der_allele_count) * (2.0 - hominin_der_allele_count))   )/ ( len_nonintrog_allele_l * 2.0 )
	print ('hominin_der_allele_freq: ', hominin_der_allele_freq, 'nonintrog_der_allele_freq: ', nonintrog_der_allele_freq, 'pi: ', pi )
	return 	pi


	
def get_pi_introg_nonintrog(population, line, line2): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	nonintrog_der_allele_count = get_non_introg_der_al_count_freq(population, line, line2)[0]
	nonintrog_der_allele_freq = get_non_introg_der_al_count_freq(population, line, line2)[2]
	len_non_introg_allele_l  = get_non_introg_der_al_count_freq(population, line, line2)[1]
	introg_der_allele_count = get_introg_der_al_count_freq(line, line2)[0]
	introg_der_allele_freq = get_introg_der_al_count_freq(line, line2)[2]
	len_introg_allele_l = get_introg_der_al_count_freq(line, line2)[1]
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
	print ('introg_der_allele_freq: ', introg_der_allele_freq, 'nonintrog_der_allele_freq: ', nonintrog_der_allele_freq, 'pi: ', pi)
	return pi


	
def get_pi_hominin_hominin(hominin1, hominin2, line): #line from the genotypes input file (inf2)
	hominin1_der_allele_freq = get_hominin_der_al_count_freq(hominin1, line)[1]
	hominin2_der_allele_freq = get_hominin_der_al_count_freq(hominin2, line)[1]
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
	if hominin1_der_allele_freq != 1 and hominin1_der_allele_freq != 0 and hominin2_der_allele_freq != 0 and hominin2_der_allele_freq != 1:
		pi = 2.0 * hominin1_der_allele_freq *  (1.0 - hominin2_der_allele_freq) 
	print('hominin1_der_allele_freq: ', hominin1_der_allele_freq, 'hominin2_der_allele_freq: ', hominin2_der_allele_freq, 'pi: ', pi)
	return pi




##Save SNP sharing data
def Sstar_snp_sharing(population, line, line2, snp_pos_shared_dict, snp_pos_population_freq_dict, notSstar_snp_pos_shared_dict, notSstar_snp_pos_population_freq_dict, notSstar_snp_pos_noYri_shared_dict, notSstar_snp_pos_noYri_population_freq_dict, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_dict, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corrected, per_haplo_noSNPs_nonSstar_noYri_count): #line: from the phased haplotypes input file (inf1), line2: from the genotypes input file (inf2)
	snps = get_Sstar_snps(line)[0]
	pos = line2[0] + '_' + line2[2]
	snp = line2[2]
	yri_freq = get_eurasian_der_al_freq('Yri', line2)
	population_freq = get_eurasian_der_al_freq(population, line2)	
	altai_freq = get_hominin_der_al_count_freq('Altai', line2)[1]
	vindija_freq = get_hominin_der_al_count_freq('Vindija', line2)[1]
	chagyrskaya_freq = get_hominin_der_al_count_freq('Chagyrskaya', line2)[1]
	denisova_freq = get_hominin_der_al_count_freq('Denisova', line2)[1]
	sharing = ''
	if altai_freq != 0:
		sharing += 'Altai' + '_'
	if vindija_freq != 0:
		sharing += 'Vindija' + '_'
	if chagyrskaya_freq != 0:
		sharing += 'Chagyrsakaya' + '_'
	if denisova_freq != 0:
		sharing += 'Denisova' + '_'
	sharing += 'shared'
	if sharing == 'shared':
		sharing = 'human_specific'
	if snp in snps:	
		snp_pos_shared_dict[pos] = sharing
		snp_pos_population_freq_dict[pos] = population_freq
	elif snp not in snps:
		notSstar_snp_pos_shared_dict[pos] = sharing
		notSstar_snp_pos_population_freq_dict[pos] = population_freq
		if yri_freq == 0:
			notSstar_snp_pos_noYri_shared_dict[pos] = sharing
			notSstar_snp_pos_noYri_population_freq_dict[pos] = population_freq
			if sharing == 'human_specific':
				count_no_derived_in_introg_haplos_per_haplo = get_introg_der_al_count_freq(line, line2)[0]/float(get_introg_der_al_count_freq(line, line2)[1])
				notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_dict[pos] = count_no_derived_in_introg_haplos_per_haplo
				per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corrected += count_no_derived_in_introg_haplos_per_haplo 			
				per_haplo_noSNPs_nonSstar_noYri_count += 1
	print ('pos: ', pos, 'sharing: ', sharing, 'population_freq: ', population_freq)
	return (snp_pos_shared_dict, snp_pos_population_freq_dict, notSstar_snp_pos_shared_dict, notSstar_snp_pos_population_freq_dict, notSstar_snp_pos_noYri_shared_dict, notSstar_snp_pos_noYri_population_freq_dict, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_dict, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corrected, per_haplo_noSNPs_nonSstar_noYri_count)

	



def main():
	count_index_error = 0
	count_value_error = 0 
	count_name_error = 0 
	all_easn_der_allele_freq_l = []
	all_easn_der_allele_freq_d = {}
	all_easn_der_allele_freq_Sstar_snp_l = []
	all_easn_der_allele_freq_notSstar_snp_l = []
	all_easn_der_allele_freq_Sstar_snp_d = {}
	all_snp_pos_shared_d = {}
	notStar_all_snp_pos_shared_d = {}
	for Line in inf:
		l = Line.split('\t')
		try:
			chr = l[0]
			chr_no = chr[3:]
			start = l[1]
			end = l[2]
			wind = l[3]
			len_l = int(l[-1].split('\n')[0])
			newSNPs = get_Sstar_snps(l)[0]
			len_haplo = get_Sstar_snps(l)[1]
			individuals = get_inds_percDerived_noHaplos(l)[0]
			percent_derived_avr = get_inds_percDerived_noHaplos(l)[1]
			no_haplos = get_inds_percDerived_noHaplos(l)[2]
		except IndexError:
			count_index_error += 1
			continue
		except ValueError:
			count_value_error += 1
			continue	
		altai_der_allele_freq_l = []
		vindija_der_allele_freq_l = []
		denisova_der_allele_freq_l = []
		chagyrskaya_der_allele_freq_l = []

		altai_der_allele_freq_Sstar_snp_l = []
		altai_der_allele_freq_notSstar_snp_l = []
		vindija_der_allele_freq_Sstar_snp_l = []
		vindija_der_allele_freq_notSstar_snp_l = []
		denisova_der_allele_freq_Sstar_snp_l = []
		denisova_der_allele_freq_notSstar_snp_l = []
		chagyrskaya_der_allele_freq_Sstar_snp_l = []
		chagyrskaya_der_allele_freq_notSstar_snp_l = []

		easn_der_allele_freq_l = []
		easn_der_allele_freq_Sstar_snp_l = []
		easn_der_allele_freq_notSstar_snp_l = [] 

		yri_der_allele_freq_l = []
		yri_der_allele_freq_Sstar_snp_l = []
		yri_der_allele_freq_notSstar_snp_l = [] 

		snp_pos_shared_d = {}
		snp_pos_population_freq_d = {}
	
		notSstar_snp_pos_shared_d = {}
		notSstar_snp_pos_population_freq_d = {}
		notSstar_all_snp_pos_shared_d = {}
	
		notSstar_snp_pos_noYri_shared_d = {} 
		notSstar_snp_pos_noYri_population_freq_d = {}
		notSstar_all_snp_pos_noYri_shared_d = {}
		notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d = {}

		per_haplo_noSNPs_nonSstar = 0.0
		per_haplo_noSNPs_nonSstar_noYri = 0.0
		per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr = 0.0

		total_pi_introg_vindija = 0.0
		total_pi_introg_altai = 0.0
		total_pi_introg_chagyrskaya = 0.0
		total_pi_introg_denisova = 0.0
		total_pi_introg_nonintrog = 0.0
		total_pi_nonintrog_nonintrog = 0.0
		total_pi_altai_nonintrogressed = 0.0
		total_pi_vindija_nonintrogressed = 0.0
		total_pi_chagyrskaya_nonintrogressed = 0.0
		total_pi_denisova_nonintrogressed = 0.0
		total_pi_altai_vindija = 0.0
		total_pi_altai_chagyrskaya = 0.0
		total_pi_altai_denisova = 0.0
		total_pi_vindija_denisova = 0.0
		total_pi_vindija_chagyrskaya = 0.0
		total_pi_chagyrskaya_denisova = 0.0
		total_pi_introg_introg = 0.0

		inf2 = open('/projects/academic/omergokc/ozgur_Sstats/double_introgression/human_neand_denis/1kgp1_altai_vindija_denisova_chagyrskaya.chr.' + chr_no + '.uniq.sorted.bed','r')	
		
		count_w_snps = 0
		for line2 in inf2:
			l2 = line2.split()
			chr2 = l2[0]
			start2 = l2[1]
			end2 = l2[2]
			ref = l2[3]
			alt = l2[4]
			if  int(end2) > int(newSNPs[-1]):
				break
			elif (int(end2) >= int(newSNPs[0]) and int(end2) <= int(newSNPs[-1])):
				count_w_snps += 1 
				pos = chr2 + '_' + end2
				try: 
					print ('snp in this haplotype %s %s %s is: %s' % (chr, start, end, end2))
				except NameError:
					count_name_error += 1
					continue
				print ('genotype line: ', l2)
				print ('get_ind_col_id_list: ', get_ind_col_id_list(l) )
				print ('get_ind_hapl_list: ', get_ind_hapl_list(l) ) 
				altai_dr_al_freq = get_hominin_der_al_count_freq('altai', l2)				
				vindija_dr_al_freq = get_hominin_der_al_count_freq('vindija', l2)				
				chagyrskaya_dr_al_freq = get_hominin_der_al_count_freq('chagyrskaya', l2)				
				denisova_dr_al_freq = get_hominin_der_al_count_freq('denisova', l2)	
				easn_dr_al_freq = get_eurasian_der_al_freq('easn', l2)			
				yri_dr_al_freq = get_eurasian_der_al_freq('yri', l2)	
				#print(snp_pos_population_freq_d)		
				Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)
# 				snp_pos_shared_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[0]
# 				snp_pos_population_freq_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[1]
# 				notSstar_snp_pos_shared_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[2]
# 				notSstar_snp_pos_population_freq_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[3]
# 				notSstar_snp_pos_noYri_shared_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[4]
# 				notSstar_snp_pos_noYri_population_freq_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[5]
# 				notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[6]
# 				per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[7]
# 				per_haplo_noSNPs_nonSstar_noYri = Sstar_snp_sharing('easn', l, l2, snp_pos_shared_d, snp_pos_population_freq_d, notSstar_snp_pos_shared_d, notSstar_snp_pos_population_freq_d, notSstar_snp_pos_noYri_shared_d, notSstar_snp_pos_noYri_population_freq_d, notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d, per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr, per_haplo_noSNPs_nonSstar_noYri)[8]																												
				
				pi_introg_altai = get_pi_introg_hominin('altai', l, l2)
				total_pi_introg_altai += pi_introg_altai
				pi_introg_vindija = get_pi_introg_hominin('vindija', l, l2)
				total_pi_introg_vindija += pi_introg_vindija
				pi_introg_chagyrskaya = get_pi_introg_hominin('chagyrskaya', l, l2)
				total_pi_introg_chagyrskaya += pi_introg_chagyrskaya
				pi_introg_denisova = get_pi_introg_hominin('denisova', l, l2)
				total_pi_introg_denisova += pi_introg_denisova

				pi_altai_vindija =  get_pi_hominin_hominin('altai', 'vindija', l2)
				total_pi_altai_vindija += pi_altai_vindija		

				pi_altai_chagyrskaya = get_pi_hominin_hominin('altai', 'chagyrskaya', l2) 
				total_pi_altai_chagyrskaya += pi_altai_chagyrskaya		

				pi_altai_denisova = get_pi_hominin_hominin('altai', 'denisova', l2)
				total_pi_altai_denisova += pi_altai_denisova	

				pi_chagyrskaya_denisova = get_pi_hominin_hominin('chagyrskaya', 'denisova', l2)
				total_pi_chagyrskaya_denisova += pi_chagyrskaya_denisova	

				pi_vindija_denisova =  get_pi_hominin_hominin('vindija', 'denisova', l2)
				total_pi_vindija_denisova += pi_vindija_denisova	

				pi_vindija_chagyrskaya = get_pi_hominin_hominin('vindija', 'chagyrskaya', l2)
				total_pi_vindija_chagyrskaya += pi_vindija_chagyrskaya	

				pi_nonintrog_nonintrog = get_pi_nonintrog_nonintrog('easn', l, l2)
				total_pi_nonintrog_nonintrog += pi_nonintrog_nonintrog 
				pi_altai_nonintrogressed = get_pi_hominin_nonintrog('altai', 'easn', l, l2)
				total_pi_altai_nonintrogressed += pi_altai_nonintrogressed 
				pi_vindija_nonintrogressed = get_pi_hominin_nonintrog('vindija', 'easn', l, l2)
				total_pi_vindija_nonintrogressed += pi_vindija_nonintrogressed
				pi_chagyrskaya_nonintrogressed = get_pi_hominin_nonintrog('chagyrskaya', 'easn', l, l2)
				total_pi_chagyrskaya_nonintrogressed += pi_chagyrskaya_nonintrogressed			 	
				pi_denisova_nonintrogressed = get_pi_hominin_nonintrog('denisova', 'easn', l, l2)
				total_pi_denisova_nonintrogressed += pi_denisova_nonintrogressed
				pi_introg_nonintrogressed = get_pi_introg_nonintrog('easn', l, l2)
				total_pi_introg_nonintrog += pi_introg_nonintrogressed
				print('len_l: ', len_l)
				
				if len_l != 7:
					pi_introg_introg = get_pi_introg_introg(l, l2)
					total_pi_introg_introg += pi_introg_introg 
		if count_w_snps > 0:
			pi_introg_altai = [ str(total_pi_introg_altai / float(len_haplo)) ]
			pi_introg_vindija = [ str(total_pi_introg_vindija / float(len_haplo)) ]
			pi_introg_chagyrskaya = [ str(total_pi_introg_chagyrskaya / float(len_haplo)) ]
			pi_introg_introg = [ str(total_pi_introg_introg / float(len_haplo)) ]
			pi_introg_nonintrog = [ str(total_pi_introg_nonintrog / float(len_haplo)) ]
			pi_nonintrog_nonintrog = [ str(total_pi_nonintrog_nonintrog/ float(len_haplo)) ]
			pi_altai_nonintrogressed = [ str(total_pi_altai_nonintrogressed/ float(len_haplo)) ]
			pi_vindija_nonintrogressed = [ str(total_pi_vindija_nonintrogressed/ float(len_haplo)) ]
			pi_chagyrskaya_nonintrogressed = [ str(total_pi_chagyrskaya_nonintrogressed/float(len_haplo)) ]
			pi_altai_vindija = [ str(total_pi_altai_vindija/ float(len_haplo)) ]
			pi_altai_chagyrskaya = [ str(total_pi_altai_chagyrskaya/float(len_haplo)) ]
			pi_vindija_chagyrskaya = [ str(total_pi_vindija_chagyrskaya/float(len_haplo)) ]
			pi_chagyrskaya_denisova = [ str(total_pi_chagyrskaya_denisova/float(len_haplo)) ]			
			pi_altai_denisova = [ str(total_pi_altai_denisova/float(len_haplo)) ]
			pi_vindija_denisova = [ str(total_pi_vindija_denisova/float(len_haplo)) ]
			pi_introg_denisova = [ str(total_pi_introg_denisova / float(len_haplo)) ]
			pi_denisova_nonintrogressed = [ str(total_pi_denisova_nonintrogressed/float(len_haplo)) ]

			#altai_der_allele_freq_l = [ str(altai_der_allele_freq_l) ]
			#vindija_der_allele_freq_l = [ str(vindija_der_allele_freq_l)]
			#denisova_der_allele_freq_l = [ str(denisova_der_allele_freq_l)]
			#chagyrskaya_der_allele_freq_l = [ str(chagyrskaya_der_allele_freq_l) ]

			#altai_der_allele_freq_Sstar_snp_l = [ str(altai_der_allele_freq_Sstar_snp_l)]
			#vindija_der_allele_freq_Sstar_snp_l = [ str(vindija_der_allele_freq_Sstar_snp_l)]
			#denisova_der_allele_freq_Sstar_snp_l = [ str(denisova_der_allele_freq_Sstar_snp_l)]
			#chagyrskaya_der_allele_freq_Sstar_snp_l = [ str(chagyrskaya_der_allele_freq_Sstar_snp_l)]

			#altai_der_allele_freq_notSstar_snp_l = [ str(altai_der_allele_freq_notSstar_snp_l)]
			#vindija_der_allele_freq_notSstar_snp_l = [ str(vindija_der_allele_freq_notSstar_snp_l)]
			#denisova_der_allele_freq_notSstar_snp_l = [ str(denisova_der_allele_freq_notSstar_snp_l)]
			#chagyrskaya_der_allele_freq_notSstar_snp_l = [ str(chagyrskaya_der_allele_freq_notSstar_snp_l)]

			#easn_der_allele_freq_l = [ str(easn_der_allele_freq_l)]
			#easn_der_allele_freq_Sstar_snp_l = [ str(easn_der_allele_freq_Sstar_snp_l)]
			#easn_der_allele_freq_notSstar_snp_l = [ str(easn_der_allele_freq_notSstar_snp_l)]

			snp_pos_shared_d = [ str(snp_pos_shared_d)]
			snp_pos_population_freq_d = [ str(snp_pos_population_freq_d)]

			notSstar_snp_pos_shared_d = [ str(notSstar_snp_pos_shared_d)]
			notSstar_snp_pos_population_freq_d = [ str(notSstar_snp_pos_population_freq_d)]
			
			notSstar_snp_pos_noYri_shared_d = [ str(notSstar_snp_pos_noYri_shared_d) ]
			notSstar_snp_pos_noYri_population_freq_d = [ str(notSstar_snp_pos_noYri_population_freq_d)]
			notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d = [ str(notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d) ]
			
			per_haplo_noSNPs_nonSstar = [str(per_haplo_noSNPs_nonSstar)]			
			per_haplo_noSNPs_nonSstar_noYri = [ str(per_haplo_noSNPs_nonSstar_noYri)]
			per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr = [ str(per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr) ] 
			
			individuals = [ str(individuals)]
			percent_derived_avr = [ str(percent_derived_avr)]
			no_haplos = [str(no_haplos)]
			
			#print>>outf, '\t'.join(l[:4] + individuals + percent_shared_avg + no_haplos + pi_introg_altai + pi_introg_vindija + pi_introg_nonintrog + pi_nonintrog_nonintrog + pi_altai_nonintrogressed + pi_vindija_nonintrogressed + pi_altai_vindija + pi_introg_introg + pi_altai_denisova + pi_vindija_denisova + pi_introg_denisova + pi_denisova_nonintrogressed + pi_altai_chagyrskaya + pi_vindija_chagyrskaya + pi_chagyrskaya_denisova + pi_introg_chagyrskaya + pi_chagyrskaya_nonintrogressed  + snp_pos_shared_d + snp_pos_easn_freq_d) 		 	
			print ('\t'.join(l[:4] + individuals + percent_derived_avr + no_haplos + pi_introg_altai + pi_introg_vindija + pi_introg_nonintrog + pi_nonintrog_nonintrog + pi_altai_nonintrogressed + pi_vindija_nonintrogressed + pi_altai_vindija + pi_introg_introg + pi_altai_denisova + pi_vindija_denisova + pi_introg_denisova + pi_denisova_nonintrogressed + pi_altai_chagyrskaya + pi_vindija_chagyrskaya + pi_chagyrskaya_denisova + pi_introg_chagyrskaya + pi_chagyrskaya_nonintrogressed), file = outf2  )	
			print ('\t'.join(l[:4] + individuals + percent_derived_avr + no_haplos + notSstar_snp_pos_shared_d + notSstar_snp_pos_population_freq_d + notSstar_snp_pos_noYri_shared_d + notSstar_snp_pos_noYri_population_freq_d + notSstar_snp_pos_noYri_count_no_der_in_intr_hapl_d + per_haplo_noSNPs_nonSstar + per_haplo_noSNPs_nonSstar_noYri + per_haplo_noSNPs_nonSstar_noYri_total_haplo_no_corr), file = outf3  )
			#print l[:4], individuals, percent_shared_avg, no_haPlos , pi_introg_altai , pi_introg_vindija , pi_introg_introg , pi_introg_nonintrog , pi_nonintrog_nonintrog , pi_altai_nonintrogressed , pi_vindija_nonintrogressed , pi_altai_vindija , pi_altai_denisova, pi_vindija_denisova, pi_introg_denisova, pi_denisova_nonintrogressed, easn_der_allele_freq_l, easn_der_allele_freq_Sstar_snp_l, snp_pos_shared_d, snp_pos_easn_freq_d, altai_der_allele_freq_Sstar_snp_l, vindija_der_allele_freq_Sstar_snp_l, denisova_der_allele_freq_Sstar_snp_l, newSNPs
	inf.close()
	outf2.close()
	outf3.close()


	


if __name__ == '__main__':
	main()

