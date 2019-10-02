from __future__ import print_function
from __future__ import division
import gzip
import os
import ast

import scipy.stats
from scipy.stats.distributions import chi2
import statsmodels.stats.contingency_tables


def read_in_data(pops, pop_tables, test):
    vardict = {}
    varids = {}
    for pop, pop_table in zip(pops, pop_tables):
        for line in pop_table:
            line = line.rstrip()
            line = line.split("\t")
            varid = line[0]
            if varid == "V":
                continue
            vardict[varid, pop] = [
                [int(line[4]), int(line[5])],
                [int(line[6]), int(line[7])],
            ]
            vardict[varid, pop + "_P"] = line[1]
            vardict[varid, pop + "_OR"] = line[2]
            vardict[varid, pop + "_SE"] = line[3]
            vardict[varid, pop + "_MAF"] = line[8]
            #HGVSp, HGVSc, consequence, gene_symbol
            if test == 'logreg':
                varids[varid] = [line[14], line[15], line[16], line[17]]
            else:
                varids[varid] = [line[13], line[14], line[15], line[16]]
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


def write_line(var, varids, vardict, valid_keys, meta_analyzable, pvalarr, orarr, mafarr, fout):
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
            varids[var][0],
            varids[var][1],
            varids[var][2],
            varids[var][3],
            min(ast.literal_eval(vardict[var, "curated.AFR_MAF"])) if "curated.AFR" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "AJ_MAF"])) if "AJ" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "curated.HISZ12_MAF"])) if "curated.HISZ12" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "curated.HISZ3_MAF"])) if "curated.HISZ3" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "FIN_MAF"])) if "FIN" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "LIT_MAF"])) if "LIT" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "NFE_MAF"])) if "NFE" in valid_keys else "",
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
            varids[var][0],
            varids[var][1],
            varids[var][2],
            varids[var][3],
            min(ast.literal_eval(vardict[var, "curated.AFR_MAF"])) if "curated.AFR" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "AJ_MAF"])) if "AJ" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "curated.HISZ12_MAF"])) if "curated.HISZ12" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "curated.HISZ3_MAF"])) if "curated.HISZ3" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "FIN_MAF"])) if "FIN" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "LIT_MAF"])) if "LIT" in valid_keys else "",
            min(ast.literal_eval(vardict[var, "NFE_MAF"])) if "NFE" in valid_keys else "",
            sep="\t",
            file=fout,
        )


def write_file(vardict, varids, pops, disease, test):
    print("Writing file for " + test + " meta-analysis of " + disease + "...")
    fout = gzip.open("output/full/CMH_" + disease + "_" + test + "_" + "_".join(pops) + ".tsv.gz", "wt")
    header = "V\tPnull\tPhet\tOR\tSE\tParr\tORarr\tMAFarr\tpoparr\tHGVSp\tHGVSc\tconsequence\tgene_symbol\t" + "\t".join([s + "_MAF" for s in ["curated.AFR", "AJ", "curated.HISZ12", "curated.HISZ3", "FIN", "LIT", "NFE"]])
    print(header, file=fout)
    line_count = 0
    for var in varids.keys():
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
            var, varids, vardict, valid_keys, meta_analyzable, pvalarr, orarr, mafarr, fout
        )
        line_count += 1
    print(
        "Wrote CMH_" + disease + "_" + test + "_" + "_".join(pops) + ".tsv.gz to output/full directory."
    )
    fout.close()


def main():
    pops = ["curated.AFR", "AJ", "curated.HISZ12", "curated.HISZ3", "FIN", "LIT", "NFE"]
    diseases = ["cd", "ibd", "uc"]
    tests = ["FET", "logreg"]
    for disease in diseases:
        for test in tests:
            print("Storing data for " + test + " meta-analysis of " + disease + "...")
            pop_tables = []
            popf = []
            for pop in pops:
                filename = ("../GWAS/output/"
                    + test + "/"
                    + pop
                    + "_"
                    + disease
                    + "_"
                    + test
                    + "_results.tsv.gz"
                )
                if os.path.exists(filename):
                    popf.append(pop)
                else:
                    continue
                print("Appending " + pop + " " + disease + " " + test + " file to valid tables...")
                pop_tables.append(
                    gzip.open(
                        filename,
                        "rt",
                    )
                )
            vardict, varids = read_in_data(popf, pop_tables, test)
            write_file(vardict, varids, popf, disease, test)

if __name__ == "__main__":
    main()
