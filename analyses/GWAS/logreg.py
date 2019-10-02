import hail as hl

# hailctl dataproc start guhancluster --vep GRCh38 --zone us-west1-b --num-preemptible-workers 100 --worker-machine-type n1-standard-16 --master-machine-type n1-highmem-32

hl.init(log="/home/gwas.log")

def compute_pcs(mt, pc_num):
    # Run PCA on synonymous, NFE <AF > 0.01, LD Pruned variants

    print("Performing variant QC and filtering on call rate...")
    mt = hl.variant_qc(mt, name="variant_qc")
    mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)

    # Lifted over version, generated with above code
    nfe_ht = hl.read_table("gs://ibd-exomes/v36meta/fin_enriched_exomes_38.ht")

    mt = mt.annotate_rows(nfe_ht=nfe_ht[mt.row_key])
    print("Filtering gnomAD NFE AF > 0.01...")
    mt = mt.filter_rows(mt.nfe_ht.nfe.AF > 0.01)
    print(mt.count())

    print("Filtering for synonymous variants...")
    mt = mt.filter_rows(mt.vep.most_severe_consequence == "synonymous_variant")
    print(mt.count())

    # LD PRUNE
    pruned_variants = hl.ld_prune(mt.GT)
    print("Pruning...")
    mt = mt.filter_rows(hl.is_defined(pruned_variants[mt.row_key]))
    print(mt.count())

    print("Computing and annotating PCs...")
    eigenvalues, pcs, _ = hl.hwe_normalized_pca(mt.GT, k=pc_num)
    
    return eigenvalues, pcs

def run_logistic_regression(test, mt, pop, disease):
    result = hl.logistic_regression_rows(
        test=test,
        y=mt.is_case,
        x=mt.GT.n_alt_alleles(),
        covariates=[
            1,
            mt.scores[0],
            mt.scores[1],
            mt.scores[2],
            mt.scores[3],
            mt.scores[4],
            mt.scores[5],
            mt.scores[6],
            mt.scores[7],
            mt.scores[8],
            mt.scores[9],
        ],
    )
    result = result.checkpoint(
        "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_" + test + "_results.ht",
        overwrite=True,
    )
    return result


def export_tsv(mt, pop, disease):
    rows = mt.rows()
    rows = rows.key_by()
    print("Annotating specific fields...")
    rows = rows.annotate(
        maf=rows.variant_qc.AF,
        call_rate=rows.variant_qc.call_rate,
        HGVSp=rows.vep.transcript_consequences.map(lambda x: x.hgvsp),
        HGVSc=rows.vep.transcript_consequences.map(lambda x: x.hgvsc),
        consequence=rows.vep.most_severe_consequence,
        mean_dp=rows.variant_qc.dp_stats.mean,
        gene_symbol=rows.vep.transcript_consequences.map(lambda x: x.gene_symbol),
    )
    rows.select(
        V=rows.V,
        P=hl.cond(rows.wald_fit.converged == True, rows.wald_p, rows.firth_p),
        OR=hl.cond(
            rows.wald_fit.converged == True,
            hl.exp(rows.wald_beta),
            hl.exp(rows.firth_beta),
        ),
        SE=hl.cond(
            rows.wald_fit.converged == True,
            rows.wald_se,
            (rows.firth_beta / (1 - hl.qnorm(rows.firth_p))),
        ),
        CaAC=rows.caac,
        CaNAC=rows.canac,
        CoAC=rows.coac,
        CoNAC=rows.conac,
        maf=rows.maf,
        call_rate=rows.call_rate,
        mean_dp=rows.mean_dp,
        case_phwe=rows.case_phwe,
        control_phwe=rows.control_phwe,
        TEST=hl.cond(rows.wald_fit.converged == True, "WALD", "FIRTH"),
        HGVSp=rows.HGVSp,
        HGVSc=rows.HGVSc,
        consequence=rows.consequence,
        gene_symbol=rows.gene_symbol,
    ).export(
        "gs://ibd-exomes/v36meta/"
        + pop
        + "_"
        + disease
        + "_logreg_results.tsv.gz"
    )



def main():
    # for pop in ["AJ", "FIN", "LIT", "NFE"]:
    for pop in ["NFE"]:
        print("Running " + pop + " Wald...")
        for disease in ["cd", "ibd", "uc"]:
        # for disease in ["cd"]:
            # Read in MT
            pop_ca_co_mt = hl.read_matrix_table(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + ".mt"
            )

            # Coalesce
            pop_ca_co_mt = pop_ca_co_mt.naive_coalesce(500)

            # Prep for GWAS with C/C status and PCs
            print("Annotating case/control status...")
            pop_ca_co_mt = pop_ca_co_mt.annotate_cols(
                is_case=(~(pop_ca_co_mt.control.contains(pop_ca_co_mt.DIAGNOSIS)))
            )

            # Checkpoint
            pop_ca_co_mt = pop_ca_co_mt.checkpoint(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_diagnosis.mt", overwrite=True
            )

            # print("Performing VEP...")
            # pop_ca_co_mt = hl.vep(pop_ca_co_mt, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")    

            eigenvalues, pcs = compute_pcs(pop_ca_co_mt, 10)

            pop_ca_co_mt = pop_ca_co_mt.annotate_cols(scores=pcs[pop_ca_co_mt.s].scores)

            # Checkpoint
            pop_ca_co_mt = pop_ca_co_mt.checkpoint(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_pcs.mt", overwrite=True
            )

            # Read in MT
            pop_ca_co_mt = hl.read_matrix_table(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_pcs.mt"
            )

            # print(pop_ca_co_mt.count())

            print(
                "Conducting Wald test for "
                + pop
                + " population, "
                + disease
                + "..."
            )

            # Run, annotate with results
            wald_result = run_logistic_regression("wald", pop_ca_co_mt, pop, disease)
            print(wald_result.count())

            wald_result = hl.read_table(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_wald_results.ht"
            )

            pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
                wald_beta=wald_result[pop_ca_co_mt.row_key].beta,
                wald_se=wald_result[pop_ca_co_mt.row_key].standard_error,
                wald_z_stat=wald_result[pop_ca_co_mt.row_key].z_stat,
                wald_p=wald_result[pop_ca_co_mt.row_key].p_value,
                wald_fit=wald_result[pop_ca_co_mt.row_key].fit,
            )

            # # Rerun with Firth
            # pop_ca_co_firth_mt = pop_ca_co_mt.filter_rows(
            #     pop_ca_co_mt.wald_fit.converged == False
            # )

            # print(
            #     "Conducting Firth test for "
            #     + pop
            #     + " population, "
            #     + disease
            #     + "..."
            # )

            # print(pop_ca_co_firth_mt.count())
            # firth_result = run_logistic_regression(
            #     "firth", pop_ca_co_firth_mt, pop, disease
            # )

            firth_result = hl.read_table(
                "gs://ibd-exomes/v36meta/" + pop + "_" + disease + "_firth_results.ht"
            )

            # Merge
            pop_ca_co_mt = pop_ca_co_mt.annotate_rows(
                firth_beta=firth_result[pop_ca_co_mt.row_key].beta,
                firth_chi_sq_stat=firth_result[pop_ca_co_mt.row_key].chi_sq_stat,
                firth_p=firth_result[pop_ca_co_mt.row_key].p_value,
                firth_fit=firth_result[pop_ca_co_mt.row_key].fit,
            )

            print("Exporting logreg results to table...")
            export_tsv(pop_ca_co_mt, pop, disease)

if __name__ == "__main__":
    main()
