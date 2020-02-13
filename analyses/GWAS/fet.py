import hail as hl

# hailctl dataproc start guhancluster --zone us-west1-b --num-preemptible-workers 148 --worker-machine-type n1-highcpu-16 --master-machine-type n1-highmem-32

hl.init(log="/home/fet.log")

pop_strings = ["AJ", "FIN", "LIT", "NFE"]
# pop_strings = ["NFE"]

for pop_string in pop_strings:
    # for disease_string in ["cd", "ibd", "uc"]:
    for disease_string in ["cd"]:
        print("Running " + pop_string + " " + disease_string + " FET...")
        pop_ca_co_mt = hl.read_matrix_table(
            "gs://ibd-exomes/v36meta/" + pop_string + "_" + disease_string + ".mt"
        )

        pop_ca_co_mt.naive_coalesce(500)

        print("Performing Fisher Exact Test...")
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
            fet=hl.fisher_exact_test(
                hl.int32(pop_ca_co_mt.caac),
                hl.int32(pop_ca_co_mt.canac),
                hl.int32(pop_ca_co_mt.coac),
                hl.int32(pop_ca_co_mt.conac),
            )
        )

        print("Annotating Standard Error and p-value...")
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
            OR=(
                ((pop_ca_co_mt.caac + 0.5) * (pop_ca_co_mt.conac + 0.5))
                / ((pop_ca_co_mt.coac + 0.5) * (pop_ca_co_mt.canac + 0.5))
            )
        )
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
            se=hl.sqrt(1.0 / (pop_ca_co_mt.caac + 0.5)
            + 1.0 / (pop_ca_co_mt.canac + 0.5)
            + 1.0 / (pop_ca_co_mt.coac + 0.5)
            + 1.0 / (pop_ca_co_mt.conac + 0.5))
        )
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(fet_p_value=pop_ca_co_mt.fet.p_value)

        pop_ca_co_mt = pop_ca_co_mt.annotate(
            HGVSp=pop_ca_co_mt.vep.transcript_consequences.map(lambda x: x.hgvsp),
            HGVSc=pop_ca_co_mt.vep.transcript_consequences.map(lambda x: x.hgvsc),
            consequence=pop_ca_co_mt.vep.most_severe_consequence,
            gene_symbol=pop_ca_co_mt.vep.transcript_consequences.map(lambda x: x.gene_symbol),
        )

        print("Exporting FET results to table...")
        rows = pop_ca_co_mt.rows()
        rows = rows.key_by()
        rows.select(
            V=rows.V,
            P=rows.P,
            OR=rows.OR,
            SE=rows.SE,
            CaAC=rows.CaAC,
            CaNAC=rows.CaNAC,
            CoAC=rows.CoAC,
            CoNAC=rows.CoNAC,
            maf=rows.maf,
            call_rate=rows.call_rate,
            mean_dp=rows.mean_dp,
            case_phwe=rows.case_phwe,
            control_phwe=rows.control_phwe,
            HGVSp=rows.HGVSp,
            HGVSc=rows.HGVSc,
            consequence=rows.consequence,
            gene_symbol=rows.gene_symbol,
        ).export(
            "gs://ibd-exomes/v36meta/"
            + pop
            + "_"
            + disease
            + "_FET_results.tsv.gz"
        )