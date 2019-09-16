import hail as hl

# hailctl dataproc start guhancluster --zone us-west1-b --num-preemptible-workers 148 --worker-machine-type n1-highcpu-16 --master-machine-type n1-highmem-32

hl.init(log="/home/wald.log")

pop_strings = ["AJ", "FIN", "LIT", "NFE"]

for pop_string in pop_strings:
    print("Running " + pop_string + " FET...")
    for disease_string in ["cd", "ibd", "uc"]:
        pop_ca_co_mt = hl.read_matrix_table(
            "gs://ibd-exomes/v36meta/" + pop_string + "_" + disease_string + ".mt"
        )

        pop_ca_co_mt.naive_coalesce(500)
        print("Annotating case/control status...")
        pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
            is_case=(~(pop_ca_co_mt.control.contains(pop_ca_co_mt.DIAGNOSIS)))
        )

        print("Computing and annotating PCs...")
        eigenvalues, pcs, _ = hl.hwe_normalized_pca(pop_ca_co_mt.GT, k=10)
        pop_ca_co_mt = pop_ca_co_mt.annotate_cols(scores=pcs[pop_ca_co_mt.s].scores)

        for test in ["wald"]:  # , 'lrt', 'firth']:
            print(
                "Conducting "
                + test
                + " test for "
                + pop_string
                + " population, "
                + disease_string
                + "..."
            )
            result = hl.logistic_regression_rows(
                test=test,
                y=pop_ca_co_mt.is_case,
                x=pop_ca_co_mt.GT.n_alt_alleles(),
                covariates=[
                    1,
                    pop_ca_co_mt.scores[0],
                    pop_ca_co_mt.scores[1],
                    pop_ca_co_mt.scores[2],
                    pop_ca_co_mt.scores[3],
                    pop_ca_co_mt.scores[4],
                    pop_ca_co_mt.scores[5],
                    pop_ca_co_mt.scores[6],
                    pop_ca_co_mt.scores[7],
                    pop_ca_co_mt.scores[8],
                    pop_ca_co_mt.scores[9],
                ],
            )
            result.export(
                "gs://ibd-exomes/"
                + pop_string
                + "_"
                + disease_string
                + "_"
                + test
                + ".tsv.gz"
            )
