import sys
import cPickle as pickle
from scipy.stats import chisquare
from copy import deepcopy
from pprint import pprint

MK_RESULTS_FILENAME = 'mk-results-01.p'
PVAL_FILENAME = 'p-val.p'
LOW_PVAL_FILENAME_PREFIX = 'low-p-val-'
LOW_HIGH_PVAL_FILENAME = 'low-high-p-val.p'
P_DIFF_FILENAME = 'p-diff.txt'

NUM_P_DIFFS = 10

NEUTRALITY_IDX = 1
CONTINGENCY_IDX = 2

CEU_IDX = 0
CHB_IDX = 1
YRI_IDX = 2

PVAL_MIN_THRESHOLD = 0.01
PVAL_MAX_THRESHOLD = 0.50


try:
    MK_RESULTS_FILENAME = sys.argv[1]
    print "Using file", MK_RESULTS_FILENAME, "as input."
except:
    print "Using default file", MK_RESULTS_FILENAME, "as input."
# Add the '01.p' onto the filename
LOW_PVAL_FILENAME = LOW_PVAL_FILENAME_PREFIX + MK_RESULTS_FILENAME[
    MK_RESULTS_FILENAME.index('.')-2:]

print "Loading MK results..."
mk_results = pickle.load(open(MK_RESULTS_FILENAME, 'rb'))

p_val = dict()
low_p_val = dict()
low_high_p_val = dict()
p_diff_to_gene = dict()

for gene in mk_results:
    p_val[gene] = []
    for i in xrange(3):
        contingency_table = mk_results[gene][CONTINGENCY_IDX][i]
        PsPn = contingency_table[0:2]
        DsDn = contingency_table[2:4]
        p_val[gene].append(chisquare(PsPn,DsDn)[1])
    min_p_val = min(p_val[gene])
    max_p_val = max(p_val[gene])
    p_diff_to_gene[max_p_val-min_p_val] = gene
    if min_p_val < PVAL_MIN_THRESHOLD:
        low_p_val[gene] = deepcopy(p_val[gene])
        if max_p_val > PVAL_MAX_THRESHOLD:
            low_high_p_val[gene] = deepcopy(p_val[gene])

print "Writing p-values to files..."
p_diff_file = open(P_DIFF_FILENAME, 'w+b')
# Write the top n genes ordered by the difference in their
# max and min p-values across populations
for gene_i in xrange(NUM_P_DIFFS):
    key = max(p_diff_to_gene)
    p_diff_file.write(p_diff_to_gene[key] + "\n")
    del p_diff_to_gene[key]
print "Wrote to file", P_DIFF_FILENAME + "..."
# pickle.dump(p_val, open(PVAL_FILENAME, 'w+b'))
# print "Wrote to file", PVAL_FILENAME + "."
pickle.dump(low_p_val, open(LOW_PVAL_FILENAME, 'w+b'))
print "Wrote to file", LOW_PVAL_FILENAME + "."
# pickle.dump(low_high_p_val, open(LOW_HIGH_PVAL_FILENAME, 'w+b'))
# print "Wrote to file", LOW_HIGH_PVAL_FILENAME + "."
num_low = len(low_p_val)
num_low_hi = len(low_high_p_val)
print "Generated", num_low, "genes with p values below the lower threshold."
print "Generated", num_low_hi, "genes with variation in p value across populations."
