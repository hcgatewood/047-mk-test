import cPickle as pickle
from matplotlib.pyplot import *
import numpy as np
import heapq

GENE_FILENAME = 'low-p-val-01.p'
# GENE_FILENAME = 'p-diff.txt'
MK_RESULTS_FILENAME = 'mk-results.p'

NEUTRALITY_IDX = 1

YRI_COLOR = '#9AAFC3'
CHB_COLOR = '#427097'
CEU_COLOR = '#174F80'


COLORS = [CEU_COLOR, CHB_COLOR, YRI_COLOR]

NUM_POPS = 3

print "Loading MK results..."
mk_results = pickle.load(open(MK_RESULTS_FILENAME, 'rb'))

print "Loading genes of interest..."
pvals = pickle.load(open(GENE_FILENAME, 'rb'))
# gene_file = open(GENE_FILENAME, 'r')

# genes = []
# for line in gene_file:
#     genes.append(line.strip())


avg_pvals = dict()
for gene in pvals:
    avg_pvals[gene] = sum(pvals[gene])/len(pvals[gene])

#average pvals for genes under positive selection
pos_sel_avg_pvals = dict()

# find genes under positive selection
for gene in avg_pvals:
    # print mk_data[gene][ALPHA_IDX]
    neutrality_indices = mk_results[gene][NEUTRALITY_IDX]
    if max(neutrality_indices) > 1:
        pos_sel_avg_pvals[gene] = avg_pvals[gene]

genes = heapq.nsmallest(10, pos_sel_avg_pvals)

NUM_GENES = len(genes)
print "Number of genes: ", NUM_GENES

WIDTH = 0.35

#Formatting plot
ax = subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)


for i in xrange(NUM_GENES):
    gene = genes[i]
    NI_values = []
    for j in xrange(NUM_POPS):
        NI_values.append((j, mk_results[gene][NEUTRALITY_IDX][j]))

    p1 = bar(i, NI_values[0][1], color = COLORS[NI_values[0][0]], edgecolor = "none")
    p2 = bar(i, NI_values[1][1], color = COLORS[NI_values[1][0]], bottom = NI_values[0][1], edgecolor = "none")
    p3 = bar(i, NI_values[2][1], color = COLORS[NI_values[2][0]], bottom = NI_values[0][1] + NI_values[1][1], edgecolor = "none")

#add labels to the x-axis
ind = np.arange(NUM_GENES)
tick_params(axis="both", which="both", bottom="off", top="off",labelbottom="on", left="off", right="off", labelleft="on")
xticks(ind + WIDTH, genes)
ylabel('Sum of Neutrality Index Across Populations')
xlabel('Gene')
title('Neutrality Index Distribution Across Populations Of Genes Under Positive Selection With P-Values < 0.01')
#title('Neutrality Index Distribution Across Populations Of Genes With Maximum Difference in P-Values')

#add a dashed line accross the plot
p4 = plot([0,NUM_GENES], [3, 3], color='k', linestyle='--', linewidth=2)

#add a legend
legend((p1[0], p2[0], p3[0], p4[0]), ('CEU', 'CHB', 'YRI', 'Null Hypothesis'), loc='best', fancybox=True, framealpha=0.5)



show()


