import hail as hl




#ADAM17: 2-9661450-A-G
#RELA: 11:65658293:C:T
#SDF2L1: 22-21998280-G-A

chrs = ["chr2", "chr11", "chr22"]
positions = [9521321, 65658293, 21643991]
refs = ["A", "C", "G"]
alts = ["G", "T", "A"]


for pop in ["AJ", "FIN", "LIT", "NFE"]:
	for disease in ["cd", "ibd", "uc"]:
		print("Reading in " + pop + " " + disease + " MT...")
		mt = hl.read_matrix_table("gs://ibd-exomes/v36meta/" + pop + "_" + disease + ".mt/")

		for chrom, position, ref, alt in zip(chrs, positions, refs, alts):
			print(chrom, position, ref, alt)
			var_mt = mt.filter_rows(
				(mt.locus == hl.locus(chrom, position, reference_genome="GRCh38"))
				& (mt.alleles == hl.array([ref, alt]))
			)
			
			print(var_mt.count())