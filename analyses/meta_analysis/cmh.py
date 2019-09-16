from __future__ import print_function
from __future__ import division
import gzip
import os

import scipy.stats
from scipy.stats.distributions import chi2
import statsmodels.stats.contingency_tables


def read_in_data(pops, pop_tables):
    vardict = {}
    varids = set()
    for pop, pop_table in zip(pops, pop_tables):
        for line in pop_table:
            line = line.rstrip()
            line = line.split("\t")
            if line[0] == "V":
                continue
            vardict[line[0], pop] = [
                [int(line[4]), int(line[5])],
                [int(line[6]), int(line[7])],
            ]
            vardict[line[0], pop + "_P"] = line[1]
            vardict[line[0], pop + "_OR"] = line[2]
            vardict[line[0], pop + "_SE"] = line[3]
            vardict[line[0], pop + "_MAF"] = line[8]
            varids.add(line[0])
    return vardict, varids


def check_if_meta_analyzable(var, vardict, pops):
    # Need at least 2 pops
    count = 0
    valid_keys = []
    for pop in pops:
        if (var, pop + "_P") in vardict:
            count += 1
            valid_keys.append(pop)
    if count > 1:
        return valid_keys, True
    else:
        return valid_keys, False


def return_cmh_results(list_of_arr):
    df = statsmodels.stats.contingency_tables.StratifiedTable(list_of_arr)
    p = chi2.sf(df.test_null_odds().statistic, 1)
    phet = chi2.sf(df.test_equal_odds().statistic, 1)
    oddsr = df.oddsratio_pooled
    se = df.logodds_pooled_se
    return p, phet, oddsr, se


def write_line(var, vardict, valid_keys, meta_analyzable, pvalarr, orarr, mafarr, fout):
    if meta_analyzable:
        list_of_arr = [vardict[var, valid_key] for valid_key in valid_keys]
        p, phet, oddsr, se = return_cmh_results(list_of_arr)
        print(
            var,
            str(p),
            str(phet),
            str(oddsr),
            str(se),
            pvalarr,
            orarr,
            mafarr,
            ",".join(valid_keys),
            sep="\t",
            file=fout,
        )
    else:
        # According to our criteria, this means that there is only 1 population in the dict, so we can safely do this
        valid_key = valid_keys[0]
        print(
            var,
            pvalarr,
            "NA",
            orarr,
            vardict[var, valid_key + "_SE"],
            pvalarr,
            orarr,
            mafarr,
            valid_key,
            sep="\t",
            file=fout,
        )


def write_file(vardict, varids, pops, disease, test):
    print("Writing file for " + test + " meta-analysis of " + disease + "...")
    fout = gzip.open("output/CMH_" + disease + "_" + "_".join(pops) + ".tsv.gz", "wt")
    print("V\tPnull\tPhet\tOR\tSE\tParr\tORarr\tMAFarr\tpoparr", file=fout)
    line_count = 0
    for var in varids:
        if line_count % 10000 == 0:
            print(line_count)
        # Count check for number of pops that have variant
        valid_keys, meta_analyzable = check_if_meta_analyzable(var, vardict, pops)
        pvalarr = ",".join(
            [str(vardict[var, valid_key + "_P"]) for valid_key in valid_keys]
        )
        orarr = ",".join(
            [str(vardict[var, valid_key + "_OR"]) for valid_key in valid_keys]
        )
        mafarr = ",".join(
            [str(vardict[var, valid_key + "_MAF"]) for valid_key in valid_keys]
        )
        write_line(
            var, vardict, valid_keys, meta_analyzable, pvalarr, orarr, mafarr, fout
        )
        line_count += 1
    print(
        "Wrote CMH_" + disease + "_" + "_".join(pops) + ".tsv.gz to output directory."
    )
    fout.close()


def main():
    pops = ["curated.AFR", "AJ", "curated.HISZ12", "curated.HISZ3", "FIN", "LIT", "NFE"]
    diseases = ["cd", "ibd", "uc"]
    # tests = ["FET", "wald"]
    tests = ["FET"]
    for disease in diseases:
        pop_tables = []
        popf = []
        for pop in pops:
            for test in tests:
                if not os.path.exists(
                    "../GWAS/output/fet/"
                    + pop
                    + "_"
                    + disease
                    + "_"
                    + test
                    + "_results.tsv.gz"
                ):
                    continue
                popf.append(pop)
                print("Opening " + pop + " " + disease + " " + test + " file...")
                pop_tables.append(
                    gzip.open(
                        "../GWAS/output/fet/"
                        + pop
                        + "_"
                        + disease
                        + "_"
                        + test
                        + "_results.tsv.gz",
                        "rt",
                    )
                )
        print("Storing data for " + test + " meta-analysis of " + disease + "...")
        vardict, varids = read_in_data(popf, pop_tables)
        write_file(vardict, varids, popf, disease, test)


if __name__ == "__main__":
    main()
