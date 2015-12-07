import sys
import scipy.stats
import cPickle as pickle
import matplotlib.pyplot as plt

PLOTS_DIRECTORY_PREFIX = 'plots/'
MK_RESULTS_01_FILENAME = 'mk-results-01.p'
MK_RESULTS_05_FILENAME = 'mk-results-05.p'
MK_RESULTS_08_FILENAME = 'mk-results-08.p'
LOW_GENES_01_FILENAME = 'low-p-val-01.p'
LOW_GENES_05_FILENAME = 'low-p-val-05.p'
LOW_GENES_08_FILENAME = 'low-p-val-08.p'

DEFAULT_POLYMORPHISM_THRESHOLD = 0.01

POPULATION_NAMES = ['CEU', 'CHB', 'YRI']

LIGHT_BLUE = '#9AAFC3'
NORMAL_BLUE = '#427097'
DARKER_BLUE = '#174F80'

# MK data idxs
ALPHA_IDX = 0
NEUTRALITY_IDX = 1

# Regression idxs
SLOPE_IDX = 0
Y_0_IDX = 1
P_VALUE_IDX = 3

# Selection functions
ALL = lambda neutr_idx: True
NEGATIVE = lambda neutr_idx: neutr_idx < 1
POSITIVE = lambda neutr_idx: neutr_idx > 1


def p_value_plot(
        included_genes, mk_data, figure_name, selection_function=ALL, use_alpha=False):
    # Generate plot values
    x_vals = []
    y_vals = []
    # for gene in mk_data:
    for gene in included_genes:
        # print mk_data[gene][ALPHA_IDX]
        for population_i in xrange(num_populations):
            alpha = mk_data[gene][ALPHA_IDX][population_i]
            neutrality_index = mk_data[gene][NEUTRALITY_IDX][population_i]
            included_value = alpha if use_alpha else neutrality_index
            # Only include these values if they meet the selection criteria
            if selection_function(neutrality_index):
                x_vals.append(population_i)
                y_vals.append(included_value)

    # Generate the regression data
    regression_data = scipy.stats.linregress(x_vals, y_vals)
    slope = regression_data[SLOPE_IDX]
    y_0 = regression_data[Y_0_IDX]
    p_value = round(regression_data[P_VALUE_IDX], 2)
    x_regression = [-1] + range(num_populations+1)
    y_regression = [el*slope + y_0 for el in x_regression]

    # Generate the plot
    # Remove the plot spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.scatter(x_vals, y_vals, color=NORMAL_BLUE)
    plt.plot(
        x_regression, y_regression,
        label='Regression over populations', color=NORMAL_BLUE)
    line1_x = [-1] + range(num_populations+1)
    line1_y = [1] * len(line1_x)
    plt.plot(
        line1_x, line1_y, '--',
        label='MK null hypothesis', color=LIGHT_BLUE)
    plt.xticks(range(num_populations), POPULATION_NAMES)
    plt.xlim(-0.5, num_populations-0.5)
    plt.ylim(ymin=0)
    plt.title(figure_name)
    plt.xlabel("Populations")
    plt.ylabel("Neutrality index")
    plt.suptitle("Regression p-value = " + str(p_value))
    plt.legend(loc='upper center', ncol=2, fancybox=True)
    # Turn off the small axis ticks
    for tic in plt.axes().xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    for tic in plt.axes().yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    plt.savefig(PLOTS_DIRECTORY_PREFIX + figure_name)
    plt.clf()


print "Loading MK results file..."
mk_file_01 = open(MK_RESULTS_01_FILENAME, 'rb')
mk_data_01 = pickle.load(mk_file_01)
mk_file_01.close()
mk_file_05 = open(MK_RESULTS_05_FILENAME, 'rb')
mk_data_05 = pickle.load(mk_file_05)
mk_file_05.close()
mk_file_08 = open(MK_RESULTS_08_FILENAME, 'rb')
mk_data_08 = pickle.load(mk_file_08)
mk_file_08.close()
print "Loading low p-value genes files..."
low_file_01 = open(LOW_GENES_01_FILENAME, 'rb')
low_genes_01_dict = pickle.load(low_file_01)
low_file_01.close()
low_file_05 = open(LOW_GENES_05_FILENAME, 'rb')
low_genes_05_dict = pickle.load(low_file_05)
low_file_05.close()
low_file_08 = open(LOW_GENES_08_FILENAME, 'rb')
low_genes_08_dict = pickle.load(low_file_08)
low_file_08.close()
print "Loaded."

low_genes_01 = low_genes_01_dict.keys()
low_genes_05 = low_genes_05_dict.keys()
low_genes_08 = low_genes_08_dict.keys()

# Get number of populations
random_gene_info = mk_data_01[mk_data_01.keys()[0]]
num_populations = len(random_gene_info[ALPHA_IDX])

print "Producing plots..."
p_value_plot(low_genes_01, mk_data_01,
             "Genes under selection, polymorphism threshold = 0%")
p_value_plot(low_genes_01, mk_data_01,
             "Genes under positive selection, polymorphism threshold = 0%",
             selection_function=POSITIVE)
p_value_plot(low_genes_01, mk_data_01,
             "Genes under negative selection, polymorphism threshold = 0%",
             selection_function=NEGATIVE)
p_value_plot(low_genes_05, mk_data_05,
             "Genes under selection, polymorphism threshold = 5%")
p_value_plot(low_genes_05, mk_data_05,
             "Genes under positive selection, polymorphism threshold = 5%",
             selection_function=POSITIVE)
p_value_plot(low_genes_05, mk_data_05,
             "Genes under negative selection, polymorphism threshold = 5%",
             selection_function=NEGATIVE)
p_value_plot(low_genes_08, mk_data_08,
             "Genes under selection, polymorphism threshold = 8%")
p_value_plot(low_genes_08, mk_data_08,
             "Genes under positive selection, polymorphism threshold = 8%",
             selection_function=POSITIVE)
p_value_plot(low_genes_08, mk_data_08,
             "Genes under negative selection, polymorphism threshold = 8%",
             selection_function=NEGATIVE)
print "Plotted."
