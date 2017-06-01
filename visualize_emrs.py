### Author: Edward Huang

import argparse
import matplotlib
import numpy as np
import os
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE

matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import pylab

### Visualizes EMR data via dimensionality reduction.

def generate_directories():
    global results_folder
    results_folder = './results'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    results_folder = './results/emr_points'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

def read_feature_matrix(suffix):
    '''
    Reads the feature matrix of the patient data. Takes an optional argument in
    the form number of dimensions.
    '''
    feature_matrix, updrs_lst, master_feature_list = [], [], []

    fname = './data/feature_matrices/feature_matrix_%s.txt' % suffix
    f = open(fname, 'r')
    for i, line in enumerate(f):
        line = line.strip().split('\t')
        feature_list = line[2:]
        if i == 0:
            master_feature_list = feature_list[:]
            continue
        assert len(feature_list) == len(master_feature_list)
        feature_matrix += [map(float, feature_list)]
        updrs_lst += [float(line[1])]
    f.close()

    return np.array(feature_matrix), updrs_lst, master_feature_list

def plot_histogram(updrs_lst):
    bins = range(200)
    matplotlib.pyplot.hist(updrs_lst, bins, alpha=0.25)

    plt.show()
    pylab.savefig('./results/emr_points/updrs_histogram.png')
    plt.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--norm_type', help='normalization method: [l1, l2, max]')
    parser.add_argument('-d', '--num_dim', help='number of ProSNet dimensions')
    parser.add_argument('-s', '--sim_thresh', help='similarity threshold for ProSNet imputation')
    parser.add_argument('-w', '--where_norm', help='where to normalize: before or after (or both) ProSNet imputation')
    parser.add_argument('-p', '--n_pca_comp', help='number of PCA components')
    parser.add_argument('-t', '--tsne_init', help='t-SNE initialization method')
    return parser.parse_args()

def main():
    args = parse_args()
    suffix = args.norm_type
    assert suffix in ['l1', 'l2', 'max']

    if args.num_dim != None:
        suffix += '_%s_%s_%s' % (args.num_dim, args.sim_thresh, args.where_norm)

    generate_directories()

    feature_matrix, updrs_lst, feature_list = read_feature_matrix(suffix)

    # TODO: currently not plotting histogram.
    # plot_histogram(updrs_lst)

    # TODO: converting arbitrary intervals of UPDRS scores to colors.
    for i, e in enumerate(updrs_lst):
        if e < 50:
            updrs_lst[i] = 'black'
        elif e < 80:
            updrs_lst[i] = 'blue'
        else:
            updrs_lst[i] = 'red'

    # PCA for ProSNet-enriched matrices.
    if args.num_dim != None:
        dim_reduc = PCA(n_components=int(args.n_pca_comp))
    else:
        # Truncated SVD on the feature matrix for non-ProSNet. TODO: ncomponents = 50?
        dim_reduc = TruncatedSVD(n_components=50, n_iter=10, random_state=9305)
    feature_matrix = dim_reduc.fit_transform(feature_matrix)

    # TODO: Use different initiation techniques.
    tsne = TSNE(n_components=2, init=args.tsne_init, random_state=9305)
    feature_matrix = tsne.fit_transform(feature_matrix)

    x_points = [point[0] for point in feature_matrix]
    y_points = [point[1] for point in feature_matrix]
    # Plot resulting feature matrix.
    plt.scatter(x=x_points, y=y_points, c=updrs_lst, s=20)
    plt.show()

    pylab.savefig('%s/%s_%s_%s.png' % (results_folder, suffix, args.n_pca_comp,
        args.tsne_init))
    plt.close()

if __name__ == '__main__':
    main()