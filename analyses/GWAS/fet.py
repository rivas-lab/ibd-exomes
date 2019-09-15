import os
import subprocess
import hail as hl

#hailctl dataproc start guhancluster --zone us-west1-b --num-workers 120 --worker-machine-type n1-standard-16 --master-machine-type n1-highmem-32 --worker-boot-disk-size 100

hl.init(log='/home/fet.log')

pop_strings = ['AJ', 'FIN', 'LIT', 'NFE']

print("Reading in QC'ed MT...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36+ccdg_qc.mt/')
print(mt.count())

print("Reading in diagnoses...")
diagnosis_info = hl.import_table('gs://ibd-exomes/v36meta/v36+ccdg_pop+diagnosis.tsv')
diagnosis_info = diagnosis_info.key_by('SAMPLE_ID')

print("Annotating MT with diagnoses...")
# WHY LOSE 1000 HERE???
mt = mt.annotate_cols(DIAGNOSIS = diagnosis_info[mt.col_key].DIAGNOSIS)
mt = mt.filter_cols(hl.is_defined(mt.DIAGNOSIS))
print(mt.count())

#FILTER OUT UNKNOWNS/BAD DIAGNOSES
print("Filtering out undesirable diagnoses...")
filter_out = hl.literal({'UNKNOWN', 'PSC', 'IBS', 'ID', 'T1D', 'AIH', 'CVID'})
mt = mt.annotate_globals(x = filter_out)
mt = mt.filter_cols(mt.x.contains(mt.DIAGNOSIS), keep=False)
print(mt.count())

shared_controls = hl.literal({'ASC', 'CONTROL'})
ibd_case = hl.literal({'CD', 'UC', 'IBD', 'IBD-U', 'IBDU', 'IC', 'ENTEROPATHY', 'ESOPHAGITIS, DUODENITIS'})

for pop_string in pop_strings:
	print("Running " + pop_string + " FET...")
	pop_table = hl.import_table('gs://ibd-exomes/v36meta/' + pop_string + '_sample_list.tsv')
	pop_table = pop_table.key_by('s')

	print("Subsetting by population...")
	pop_mt = mt.filter_cols(hl.is_defined(pop_table[mt.s]), keep=True)
	print(pop_mt.count())

	print("Splitting samples into diagnosis sets for:")
	print("CD...")
	pop_mt = pop_mt.annotate_globals(cd = hl.literal({'CD'}))
	pop_cd_mt = pop_mt.filter_cols(pop_mt.cd.contains(pop_mt.DIAGNOSIS), keep=True)
	print("IBD...")
	pop_mt = pop_mt.annotate_globals(ibd = ibd_case)
	pop_ibd_mt = pop_mt.filter_cols(pop_mt.ibd.contains(pop_mt.DIAGNOSIS), keep=True)
	print("UC...")
	pop_mt = pop_mt.annotate_globals(uc = hl.literal({'UC'}))
	pop_uc_mt = pop_mt.filter_cols(pop_mt.uc.contains(pop_mt.DIAGNOSIS), keep=True)
	print("Controls...")
	pop_mt = pop_mt.annotate_globals(control = shared_controls)
	pop_control_mt = pop_mt.filter_cols(pop_mt.control.contains(pop_mt.DIAGNOSIS), keep=True)

	print("Performing variant QC on whole table and shared controls...")
	pop_mt = hl.variant_qc(pop_mt, name='variant_qc')
	pop_control_mt = hl.variant_qc(pop_control_mt, name='variant_qc')

	for disease_string, pop_case_mt in zip(['cd', 'ibd', 'uc'], [pop_cd_mt, pop_ibd_mt, pop_uc_mt]):
		print("Performing variant QC for " + disease_string + " cases...")
		pop_case_mt = hl.variant_qc(pop_case_mt, name='variant_qc')

		print('Annotating case and control pHWE...')
		pop_mt = pop_mt.annotate_rows(case_phwe = pop_case_mt.index_rows(pop_mt.row_key).variant_qc.p_value_hwe)
		pop_mt = pop_mt.annotate_rows(control_phwe = pop_control_mt.index_rows(pop_mt.row_key).variant_qc.p_value_hwe)

		print('Annotating CaAC, CoAC, CaNAC, CoNAC...')
		pop_mt = pop_mt.annotate_rows(caac = hl.agg.filter(pop_mt[disease_string].contains(pop_mt.DIAGNOSIS), hl.agg.sum(pop_mt.GT.n_alt_alleles())))
		pop_mt = pop_mt.annotate_rows(coac = hl.agg.filter(pop_mt['control'].contains(pop_mt.DIAGNOSIS), hl.agg.sum(pop_mt.GT.n_alt_alleles())))
		pop_mt = pop_mt.annotate_rows(canac = hl.agg.filter(pop_mt[disease_string].contains(pop_mt.DIAGNOSIS), hl.agg.sum(pop_mt.GT.is_hom_ref() * 2 + pop_mt.GT.is_het_ref())))
		pop_mt = pop_mt.annotate_rows(conac = hl.agg.filter(pop_mt['control'].contains(pop_mt.DIAGNOSIS), hl.agg.sum(pop_mt.GT.is_hom_ref() * 2 + pop_mt.GT.is_het_ref())))

		print('Filtering on CaAC and CoAC...')
		pop_mt = pop_mt.filter_rows((pop_mt.caac > 0) | (pop_mt.coac > 0))

		print('Performing Fisher Exact Test...')
		pop_mt = pop_mt.annotate_rows(fet = hl.fisher_exact_test(hl.int32(pop_mt.caac), hl.int32(pop_mt.canac), hl.int32(pop_mt.coac), hl.int32(pop_mt.conac)))

		print('Annotating Standard Error...')
		pop_mt = pop_mt.annotate_rows(OR = (((pop_mt.caac + 0.5) * (pop_mt.conac + 0.5)) / ((pop_mt.coac + 0.5) * (pop_mt.canac + 0.5))))
		pop_mt = pop_mt.annotate_rows(se = 1.0/(pop_mt.caac + 0.5) + 1.0/(pop_mt.canac + 0.5) + 1.0/(pop_mt.coac + 0.5) + 1.0/(pop_mt.conac + 0.5))

		print('Exporting population table for other analyses...')
		pop_mt = pop_mt.annotate_rows(fet_p_value = pop_mt.fet.p_value)
		pop_mt = pop_mt.annotate_rows(CHROM = pop_mt.locus.contig.replace('chr', ''), POS = pop_mt.locus.position, REF = pop_mt.alleles[0], ALT = pop_mt.alleles[1])
		pop_mt = pop_mt.annotate_rows(V = hl.delimit(hl.array([pop_mt.CHROM, hl.str(pop_mt.POS), pop_mt.REF, pop_mt.ALT]), ":"))

		pop_mt.write('gs://ibd-exomes/v36meta/' + pop_string + '_' + disease_string + '.mt', overwrite=True)

		print('Exporting FET results to table...')
		rows = pop_mt.rows()
		rows.select(V = rows.V, P = rows.fet_p_value, OR = rows.OR, SE = rows.se, CaAC = rows.caac, CaNAC = rows.canac, CoAC = rows.coac, CoNAC = rows.conac, maf = rows.variant_qc.AF, call_rate = rows.variant_qc.call_rate, mean_dp = rows.variant_qc.dp_stats.mean, case_phwe = rows.case_phwe, control_phwe = rows.control_phwe).export('gs://ibd-exomes/' + pop_string + '_' + disease_string + '_FET_results.tsv.bgz')
