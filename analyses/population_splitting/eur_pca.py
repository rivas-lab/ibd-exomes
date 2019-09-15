import hail as hl

hl.init(log='/home/eur_pca.log')
# mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36+ccdg_082119.mt')

# #Pull out EUR PCA samples
# print("Grabbing train/test European samples for PCA...")
# eur_ht = hl.import_table('gs://ibd-exomes/v36meta/eur_samples_to_pca.tsv')
# eur_ht = eur_ht.key_by('s')
# mt = mt.filter_cols(hl.is_defined(eur_ht[mt.s]), keep=True)

# #Get rid of duplicate ID samples
# problematic_samples = hl.import_table('gs://ibd-exomes/v36meta/problematic_samples.tsv')
# problematic_samples = problematic_samples.key_by('s')
# mt = mt.filter_cols(hl.is_defined(problematic_samples[mt.s]), keep=False)
# print(mt.count())

# #Filter on call rate
# print("Performing variant QC and filtering on call rate...")
# mt = hl.variant_qc(mt, name='variant_qc')
# mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_eur_call_rate.mt', overwrite=True)
# print(mt.count())

# #Filter on gnomAD exome NFE MAF, Konrad's table contains NFE AF
# #Lifted over version
# nfe_ht = hl.read_table('gs://ibd-exomes/v36meta/fin_enriched_exomes_38.ht')
# mt = mt.annotate_rows(nfe_ht = nfe_ht[mt.row_key])

# print("Filtering gnomAD NFE AF > 0.01...")
# mt = mt.filter_rows(mt.nfe_ht.nfe.AF > 0.01)
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_eur_nfe.mt', overwrite=True)
# print(mt.count())

# print("Performing VEP...")
# mt = hl.vep(mt, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")
# print("Filtering for synonymous variants...")
# mt = mt.filter_rows(mt.vep.most_severe_consequence == "synonymous_variant")
# mt = mt.naive_coalesce(500)
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_eur_filtered.mt', overwrite=True)
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36+ccdg_eur_filtered.mt')
print(mt.count())

#LD PRUNE
pruned_variants = hl.ld_prune(mt.GT)
print("Pruning...")
mt = mt.filter_rows(hl.is_defined(pruned_variants[mt.row_key]))
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_eur_pruned.mt', overwrite=True)

#PCA on the EUR samples
print("Performing PCA...")
eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=40)

print("Exporting results...")
scores.export('gs://ibd-exomes/v36meta/eur_40_pca_scores.tsv.bgz')
print(eigenvalues)