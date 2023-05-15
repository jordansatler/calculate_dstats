"""Calculate D-statistics (ABBA/BABA) for four taxon trees.
Assumes SNPs are unlinked and in phylip format. Also requires
a phylogeny and a taxon designated as the outgroup. Script uses
1000 bootstrap replicates to test significance of D-stat patterns. 

usage: python dstats.py snps.file tree.file outgroup
"""
import sys
import random
import datetime
import numpy as np
from ete3 import Tree
from scipy import stats
from itertools import combinations

def read_data(file):
    """read snp input file"""
    with open(file, "r") as f:
        return{i.split()[0]:i.split()[1] for i in f if not i.split()[1].isnumeric()}

def get_taxa_combinations(tr, OG, taxa):
    """get all test combinations ordered by phylogeny"""
    tree = Tree(tr)

    # check if outgroup is in tree
    if not OG in tree.get_leaf_names():
        print("OUTGROUP is not in tree. Please provide outgroup species.\n")
        print("SAMPLED TAXA: {0}".format(', '.join(i for i in 
                                                   sorted(tree.get_leaf_names()))))
        sys.exit()
    
    # get all combinations of 3 ingroup taxa
    tree.set_outgroup(OG)
    ingroup = [i for i in taxa if not i == OG]
    sp_sets = list(combinations(ingroup, 3))

    # order 3 ingroup taxa + outgroup by phylogeny
    res = []
    for sp in sp_sets:
        taxa_set = [t for t in sp]
        taxa_set.append(OG)
        tree2 = tree.copy("newick")
        tree2.prune(taxa_set)
        # O,P3,P2,P1 order
        tax = [node.name for node in tree2.traverse("levelorder") if node.name]
        res.append(tax)
    return res

def get_valid_snps(data, taxon_set):
    """retain snps that match ABBA/BABA pattern"""
    snps = []
    t = [data[i].upper() for i in taxon_set]
    allowed = ['A', 'C', 'G', 'T']
    for i in range(len(t[0])):
        snp = [ind[i] for ind in t]
        # retain snp if matches ABBA or BABA pattern
        if len(set(snp).intersection(allowed)) == 2 and 'N' not in snp:
            # ABBA
            if snp[0] == snp[3] and snp[1] == snp[2] and snp[0] != snp[1]:
                snps.append('ABBA')
            # BABA
            elif snp[0] == snp[2] and snp[1] == snp[3] and snp[0] != snp[1]:
               snps.append('BABA')
    return snps

def get_z_score(Dvals):
    """calculate z-score for replicates"""
    mean = np.mean(Dvals)
    std = np.std(Dvals)
    Z = (0 - mean) / std
    pval = stats.norm.sf(abs(Z)) * 2
    return mean, std, Z, pval

def write_to_file(res):
    """write results to file"""
    with open("Dstats_results.txt", "w") as out:
        h = "O\tP3\tP2\tP1\tNsnps\tZ\tp_value\tSig?\tmean\tstdev"
        out.write(h + "\n")
        sig_level = 0.01 / len(res)
        for i in sorted(res, key = lambda x: x[3]):
            out.write("{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4}\t{5:.4f}\t{6:.4f}\n".format('\t'.join(j for j in i[0]),
                                                                                   i[1],
                                                                                   i[2],
                                                                                   i[3],
                                                                                   "NS" if i[3] > sig_level else "SIG",
                                                                                   i[4],
                                                                                   i[5]))

def summarize_significant_results(res):
    """summarize patterns of significant results"""
    sum_res = {}
    for line in res:
        if line[3] < (0.01 / len(res)):
            tax_set = '\t'.join(sorted([line[0][1], line[0][3]]))
            if line[4] < 0:
                if tax_set not in sum_res:
                    sum_res[tax_set] = 1
                else:
                    sum_res[tax_set] += 1
            else:
                if tax_set not in sum_res:
                    sum_res[tax_set] = 1
                else:
                    sum_res[tax_set] += 1
    
    with open("Dstats_results_summarized.txt", "w") as out:
        h = "t1\tt2\tcount\n"
        out.write(h)
        # sort dictionary entries by counts
        sum_res_sort = {k:v for k, v in sorted(sum_res.items(), 
                        key = lambda x: x[1], reverse = True)}
        for k, v in sum_res_sort.items():
            out.write("{0}\t{1}\n".format(k, v))

def main():
    if len(sys.argv) != 4:
        print("python dstats.py snps.file tree.file outgroup")
        sys.exit()
    
    data = read_data(sys.argv[1])
    taxa_set = get_taxa_combinations(sys.argv[2], sys.argv[3], data.keys())
    
    # keep track of time
    begin_time = datetime.datetime.now()
    
    # store results
    res = []
    for idx, test in enumerate(taxa_set):
        filtered_snps = get_valid_snps(data, test)
        
        # if no SNPs match ABBA or BABA pattern
        if not filtered_snps:
            continue

        print("{0}\t{1}/{2}\t{3}".format(datetime.datetime.now(),
                                         idx + 1,
                                         len(taxa_set),
                                         ', '.join(i for i in test)))

        D = []
        # do 1000 replicates
        for i in range(1000):
            ABBA = 0
            BABA = 0
            # bootstrap replicate matching empirical number of valid snps
            for j in range(len(filtered_snps)):
                picked_snp = random.choice(filtered_snps)
                if picked_snp.startswith("A"):
                    ABBA += 1
                else:
                    BABA += 1
            d = (ABBA - BABA) / (ABBA + BABA)
            D.append(d)
        mu, sd, Z, p = get_z_score(D)
        res.append([test, len(filtered_snps), Z, p, mu, sd])

    write_to_file(res)
    summarize_significant_results(res)
    print("Total time: {0}".format(datetime.datetime.now() - begin_time))

if __name__ == '__main__':
    main()
