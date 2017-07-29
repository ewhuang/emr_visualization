### Author: Edward Huang

import argparse
import json
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

### Given a file, plot the fraction of clusters enriched with at least one
### biomarker across FDRs.

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--prosnet_clus_fname', required=True, type=str,
        help='File name for ProSNet json-dumped dictionary mapping clusters to enriched biomarkers.')
    parser.add_argument('-b', '--baseline_clus_fname', required=True, type=str,
        help='File name for baseline json-dumped dictionary mapping clusters to enriched biomarkers.')
    args = parser.parse_args()

def compute_frac_enriched(fdr, clus_p_lst):
    '''
    Compute the number of clusters enriched with at least one biomarker for the
    given FDR.
    '''
    num_enriched_clusters = 0.0
    sorted_clus_p_lst = sorted(clus_p_lst)
    m = len(sorted_clus_p_lst)
    for bm_idx, p_value in enumerate(sorted_clus_p_lst):
        # Compute the upper bound.
        upper_bound = (bm_idx + 1.0) / m * fdr
        if p_value <= upper_bound:
            num_enriched_clusters += 1
    return num_enriched_clusters

def plot_method(fname, color, label):
    with open(fname, 'r') as fp:
        clus_enrichment_dct = json.load(fp)
    fp.close()

    # print len(clus_enrichment_dct)

    # Get the lowest p-value for each biomarker.
    bm_best_p_dct = {}
    for clus_id in clus_enrichment_dct:
        biomarker_p_lst = clus_enrichment_dct[clus_id]
        for biomarker, f_table, p_value in biomarker_p_lst:
            if biomarker in bm_best_p_dct:
                p_value = min(p_value, bm_best_p_dct[biomarker])
            bm_best_p_dct[biomarker] = p_value
    print len(clus_enrichment_dct)
    clus_p_lst = []
    for clus_id in clus_enrichment_dct:
        biomarker_p_lst = clus_enrichment_dct[clus_id]
        # print biomarker_p_lst
        # Sort the list of tuples by the p-values.
        curr_p_lst = sorted(biomarker_p_lst, key=lambda x:x[2])
        clus_p_lst += [curr_p_lst[0][2]]
        print curr_p_lst[:50]
        # # Keep going until a tuple matches the best p-value for that biomarker.
        # for biomarker, f_table, p_value in curr_p_lst:
        #     if bm_best_p_dct[biomarker] == p_value:
        #         clus_p_lst += [p_value]
        #         break

    plot_point_lst = []
    for fdr in np.arange(0, 0.06, 0.001):
        frac_enriched_clusters = compute_frac_enriched(fdr, clus_p_lst)
        plot_point_lst += [(fdr, frac_enriched_clusters)]

    plt.plot(*zip(*plot_point_lst), color=color, label=label, lw=2.5)
    return len(clus_p_lst)

def main():
    parse_args()

    plt.figure(figsize=(10, 7.5))

    n_pros_clus = plot_method(args.prosnet_clus_fname, '#7CAE00', 'prosnet')
    n_base_clus = plot_method(args.baseline_clus_fname, '#00BFC4', 'baseline')

    # Plot settings.
    plt.ylim(0, 40)
    plt.xlim(0, 0.05)
    # plt.legend(loc='lower right')
    plt.text(0.025, 41, 'FDR of Each Method\'s Biomarker Enrichments', fontsize=17,
        ha='center')
    plt.xlabel('False Discovery Rate')
    plt.ylabel('Number of Biomarker-Enriched Clusters')

    # Experimenting with plot prettiness.
    # Remove the plot frame lines. They are unnecessary chartjunk.    
    ax = plt.subplot(111)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    # Ensure that the axis ticks only show up on the bottom and left of the plot.    
    # Ticks on the right and top of the plot are generally unnecessary chartjunk.    
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left()    
    for y in range(5, 41, 5):    
        plt.plot(np.arange(0, 2), [y] * len(np.arange(0, 2)), "--", lw=0.5, color="black", alpha=0.3)    
    # Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
    plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")    

    plt.text(0.05, n_base_clus, 'baseline', fontsize=14, color='#00BFC4')
    plt.text(0.05, n_pros_clus, 'ProSNet', fontsize=14, color='#7CAE00')

    plt.show()
    plt.savefig('./results/biomarker_enrichments/fdr_plot.pdf')
    plt.close()

if __name__ == '__main__':
    main()