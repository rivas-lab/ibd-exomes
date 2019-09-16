import hail as hl

# hailctl dataproc start guhancluster --zone us-west1-b --num-preemptible-workers 148 --worker-machine-type n1-highcpu-16 --master-machine-type n1-highmem-32

hl.init(log="/home/fet.log")

# pop_strings = ["AJ", "FIN", "LIT", "NFE"]
pop_strings = ["NFE"]

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
            se=1.0 / (pop_ca_co_mt.caac + 0.5)
            + 1.0 / (pop_ca_co_mt.canac + 0.5)
            + 1.0 / (pop_ca_co_mt.coac + 0.5)
            + 1.0 / (pop_ca_co_mt.conac + 0.5)
        )
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(fet_p_value=pop_ca_co_mt.fet.p_value)

        print("Exporting FET results to table...")
        rows = pop_ca_co_mt.rows()
        rows.select(
            V=rows.V,
            P=rows.fet_p_value,
            OR=rows.OR,
            SE=rows.se,
            CaAC=rows.caac,
            CaNAC=rows.canac,
            CoAC=rows.coac,
            CoNAC=rows.conac,
            maf=rows.variant_qc.AF,
            call_rate=rows.variant_qc.call_rate,
            mean_dp=rows.variant_qc.dp_stats.mean,
            case_phwe=rows.case_phwe,
            control_phwe=rows.control_phwe,
        ).export(
            "gs://ibd-exomes/v36meta/"
            + pop_string
            + "_"
            + disease_string
            + "_FET_results.tsv.gz"
        )
