import hail as hl

# Adding liftover capabilities
rg37 = hl.get_reference("GRCh37")
rg38 = hl.get_reference("GRCh38")
rg37.add_liftover("gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38)

# Loop through summary statistics and assign new loci
# for disease in ["cd", "ibd", "uc"]:
#     for pop in ["AFR", "HISZ12", "HISZ3"]:
#         tbl = hl.import_table(
#             "gs://ibd-exomes/v36meta/curated."
#             + pop
#             + "_"
#             + disease
#             + "_FET_results.tsv.gz",
#             force=True,
#         )
#         print(tbl.count())
#         tbl = tbl.annotate(variant=hl.parse_variant(tbl.V, reference_genome="GRCh37"))
#         tbl = tbl.annotate(locus=tbl.variant.locus, alleles=tbl.variant.alleles)
#         tbl = tbl.annotate(new_locus=hl.liftover(tbl.locus, "GRCh38"))
#         tbl = tbl.filter(hl.is_defined(tbl.new_locus))
#         tbl = tbl.annotate(
#             new_V=hl.variant_str(tbl.new_locus, tbl.alleles).replace("chr", "")
#         )
#         tbl.select(
#             V=tbl.new_V,
#             P=tbl.P,
#             OR=tbl.OR,
#             SE=tbl.SE,
#             CaAC=tbl.CaAC,
#             CaNAC=tbl.CaNAC,
#             CoAC=tbl.CoAC,
#             CoNAC=tbl.CoNAC,
#             maf=tbl.maf,
#             call_rate=tbl.call_rate,
#             mean_dp=tbl.mean_dp,
#             case_phwe=tbl.case_phwe,
#             control_phwe=tbl.control_phwe,
#         ).export(
#             "gs://ibd-exomes/v36meta/curated."
#             + pop
#             + "_"
#             + disease
#             + "_FET_results_liftover.tsv.gz"
#         )
tbl = hl.import_table("gs://ukbb-exome/ld_indep_array.tsv")
tbl = tbl.annotate(variant=hl.parse_variant(tbl.V, reference_genome="GRCh37"))
tbl = tbl.annotate(locus=tbl.variant.locus, alleles=tbl.variant.alleles)
tbl = tbl.annotate(new_locus=hl.liftover(tbl.locus, "GRCh38"))
tbl = tbl.filter(hl.is_defined(tbl.new_locus))
tbl = tbl.annotate(
    new_V=hl.variant_str(tbl.new_locus, tbl.alleles).replace("chr", "")
)
tbl.select(
    V=tbl.new_V,
).export("gs://ukbb-exome/ld_indep_exome.tsv")