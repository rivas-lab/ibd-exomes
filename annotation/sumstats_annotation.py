import hail as hl

# vep GRCh38

# Import mt file
print("Reading in sumstats...")

hl.init(log='/home/sumstats_annotation.log')

# pops = ["curated.AFR", "AJ", "curated.HISZ12", "curated.HISZ3", "FIN", "LIT", "NFE"]
pops = ["curated.HISZ12"]

# diseases = ["cd", "ibd", "uc"]
diseases = ["ibd", "uc"]

for pop in pops:
    for disease in diseases:
        filename = (
            "gs://ibd-exomes/v36meta/"
            + pop
            + "_"
            + disease
            + "_FET_results.tsv.gz"
        )
        print("Importing " + filename + "...")
        tbl = hl.import_table(filename, force_bgz=True, min_partitions=100)
        tbl = tbl.annotate(variant=hl.parse_variant("chr" + tbl.V, reference_genome="GRCh38"))
        tbl = tbl.annotate(locus=tbl.variant.locus, alleles=tbl.variant.alleles)
        tbl = tbl.key_by(tbl.locus, tbl.alleles)

        # mt = hl.read_matrix_table(
        #     "gs://ibd-exomes/v36meta/"
        #     + pop
        #     + "_"
        #     + disease
        #     + ".mt"
        # )

        # tbl = tbl.annotate(variant=hl.parse_variant("chr" + tbl.V, reference_genome="GRCh38"))
        # tbl = tbl.annotate(locus=tbl.variant.locus, alleles=tbl.variant.alleles)
        # tbl = tbl.key_by(tbl.locus, tbl.alleles)
        # print("Annotating with VEP....")
        # tbl = tbl.annotate(vep=mt.index_rows(tbl.key).vep)

        tbl = tbl.checkpoint('gs://ibd-exomes/v36meta/tbl.ht', overwrite=True)
        
        print("Annotating with VEP....")
        tbl = hl.vep(tbl, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")
        print("Annotating specific fields...")
        tbl = tbl.annotate(
            HGVSp=tbl.vep.transcript_consequences.map(lambda x: x.hgvsp),
            HGVSc=tbl.vep.transcript_consequences.map(lambda x: x.hgvsc),
            consequence=tbl.vep.most_severe_consequence,
            gene_symbol=tbl.vep.transcript_consequences.map(lambda x: x.gene_symbol),
        )
        tbl = tbl.key_by()
        print("Exporting...")
        tbl.select(
            V=tbl.V,
            P=tbl.P,
            OR=tbl.OR,
            SE=tbl.SE,
            CaAC=tbl.CaAC,
            CaNAC=tbl.CaNAC,
            CoAC=tbl.CoAC,
            CoNAC=tbl.CoNAC,
            maf=tbl.maf,
            call_rate=tbl.call_rate,
            mean_dp=tbl.mean_dp,
            case_phwe=tbl.case_phwe,
            control_phwe=tbl.control_phwe,
            HGVSp=tbl.HGVSp,
            HGVSc=tbl.HGVSc,
            consequence=tbl.consequence,
            gene_symbol=tbl.gene_symbol,
        ).export(
            "gs://ibd-exomes/v36meta/"
            + pop
            + "_"
            + disease
            + "_FET_results_annotated.tsv.gz"
        )
