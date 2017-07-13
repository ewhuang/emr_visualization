### Author: Edward Huang

import argparse
import matplotlib
import numpy as np
import os
from process_loni_parkinsons import get_updrs_dct
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier

matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
import pylab

### Visualizes EMR data via dimensionality reduction.

def generate_directories():
    for results_folder in ('./results', './results/emr_points',
        './results/baseline_accuracy', './results/prosnet_accuracy'):
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)

def read_feature_matrix(suffix):
    '''
    Reads the feature matrix of the patient data. Takes an optional argument in
    the form number of dimensions.
    '''
    # feature_matrix, label_lst, master_feature_list = [], [], []
    feature_matrix, patient_lst = [], []

    fname = './data/feature_matrices/feature_matrix_%s.tsv' % suffix
    f = open(fname, 'r')
    header = f.readline().strip().split('\t')
    for line in f:
        line = line.strip().split('\t')
        assert len(line) == len(header)
        feature_matrix += [map(float, line[1:])]
        patient_lst += [line[0]]
        # TODO: SWEDD/Control/PD won't be floats.
        # label_lst += [float(line[1])]
        # label_lst += [line[1]]
    f.close()

    # return np.array(feature_matrix), label_lst, master_feature_list
    return np.array(feature_matrix), patient_lst

# def plot_histogram(label_lst):
#     bins = range(200)
#     matplotlib.pyplot.hist(label_lst, bins, alpha=0.25)

#     plt.show()
#     pylab.savefig('./results/emr_points/updrs_histogram.png')
#     plt.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--norm_type', help='normalization method: [l1, l2, max]')
    parser.add_argument('-d', '--num_dim', help='number of ProSNet dimensions')
    parser.add_argument('-s', '--sim_thresh', help='similarity threshold for ProSNet imputation')
    parser.add_argument('-w', '--where_norm', help='where to normalize: before or after (or both) ProSNet imputation')
    parser.add_argument('-p', '--n_pca_comp', help='number of PCA components')
    parser.add_argument('-t', '--tsne_init', help='t-SNE initialization method')
    parser.add_argument('-r', '--l_rate', help='learning rate of t-SNE')
    # parser.add_argument('-k', '--knn', help='number of nearest neighbors')
    return parser.parse_args()

def knn_predict(args, fname, X, y):
    '''
    Gets the classification accuracy on a 3-fold CV with a KNN classifier.
    Writes out the score to file if this is a baseline run.
    '''
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,
    #     random_state=9305)
    # TODO: make n_neighbors be a command line option.
    knn = KNeighborsClassifier(n_neighbors=10)
    # knn.fit(X_train, y_train)
    # pred = knn.predict(X_test)

    # score = accuracy_score(y_test, pred)
    # print score
    # TODO: 10-fold CV.
    score = cross_val_score(knn, X, y, cv=10, scoring='accuracy').mean()

    # Write the score out to file if this is a baseline run.
    if args.num_dim == None:
        folder = './results/baseline_accuracy'
    else:
        folder = './results/prosnet_accuracy'
    out = open('%s/%s.txt' % (folder, fname), 'w')
    out.write('%g' % score)
    out.close()
    return score

def label_visualization(args, feature_matrix, label_lst, suffix, test_name):
    # PCA for ProSNet-enriched matrices.
    # if args.num_dim != None:
    # TODO: everyone running Truncated.
    if True:
        # print 'running PCA...'
        dim_reduc = PCA(n_components=int(args.n_pca_comp))
    else:
        print 'running truncated SVD...'
        # Truncated SVD on the feature matrix for non-ProSNet. TODO: ncomponents = 50?
        dim_reduc = TruncatedSVD(n_components=int(args.n_pca_comp), n_iter=10,
            random_state=9305)
    feature_matrix = dim_reduc.fit_transform(feature_matrix)

    tsne = TSNE(n_components=2, init=args.tsne_init, random_state=9305,
        learning_rate=int(args.l_rate))
    feature_matrix = tsne.fit_transform(feature_matrix)

    # Predict the accuracy of the visualization.
    fname = '%s_%s_%s_%s_%s' % (suffix, args.n_pca_comp, args.tsne_init,
        test_name, args.l_rate)
    accuracy = knn_predict(args, fname, feature_matrix, label_lst)
    print accuracy

    # Always plot baseline plots.
    if args.num_dim == None: # Title for baseline is just the accuracy score.
        title = 'baseline accuracy (%g)' % accuracy
        plot_embeddings(feature_matrix, label_lst, fname, title)
    # Plot a ProSNet plot if its accuracy is better than its corresponding baseline.
    else:
        baseline_score = get_baseline_accuracy(args, test_name)
        if accuracy > baseline_score:
            title = '%g vs. baseline (%g)' % (accuracy, baseline_score)
            plot_embeddings(feature_matrix, label_lst, fname, title)

def get_baseline_accuracy(args, test_name):
    '''
    Get the classification accuracy of the baseline method.
    '''
    fname = '%s_%s_%s_%s_%s' % (args.norm_type, args.n_pca_comp,
        args.tsne_init, test_name, args.l_rate)
    f = open('./results/baseline_accuracy/%s.txt' % fname, 'r')
    score = float(f.readline())
    f.close()
    return score

def plot_embeddings(feature_matrix, label_lst, fname, title):
    x_points = [point[0] for point in feature_matrix]
    y_points = [point[1] for point in feature_matrix]
    # Plot resulting feature matrix.
    plt.scatter(x=x_points, y=y_points, c=label_lst, s=20, alpha=0.5)
    plt.title(title)
    plt.show()

    pylab.savefig('./results/emr_points/%s.png' % fname)
    plt.close()

def main():
    generate_directories()
    args = parse_args()

    # Process the filename for each run.
    suffix = args.norm_type # Baseline only contains normalization method.
    assert suffix in ['l1', 'l2', 'max']
    if args.num_dim != None:
        assert args.num_dim.isdigit()
        suffix += '_%s_%s_%s' % (args.num_dim, args.sim_thresh, args.where_norm)

    # feature_matrix, label_lst, feature_list = read_feature_matrix(suffix)
    feature_matrix, patient_lst = read_feature_matrix(suffix)

    label_dct, score_names = get_updrs_dct()
    assert set(label_dct.keys()) == set(patient_lst)

    # Get the label list from the label dictionary.
    # 1. Sum of UPDRS scores.
    label_lst = []
    for patno in patient_lst:
        score = sum(label_dct[patno].values())
        label_lst += [score]
    for i, e in enumerate(label_lst):
        # Convert UPDRS scores to colors based on their range.
        if e < 45:
            label_lst[i] = 'white'
        elif e < 70:
            label_lst[i] = 'yellow'
        else:
            label_lst[i] = 'red'
    label_visualization(args, feature_matrix, label_lst, suffix, 'sum')

    # # 2. Individual UPDRS scores. # TODO: Only doing NP2WALK.
    # for score_name in score_names:
    # # for score_name in ['NP4OFF']:
    #     label_lst = []
    #     for patno in patient_lst:
    #         if score_name not in label_dct[patno]:
    #             score = 0
    #         else:
    #             score = label_dct[patno][score_name]
    #         label_lst += [score]
    #     # Convert labels to colors for each UPDRS test.
    #     for i, e in enumerate(label_lst):
    #         if e < 1:
    #             label_lst[i] = 'black'
    #         elif e < 3:
    #             label_lst[i] = 'blue'
    #         else:
    #             label_lst[i] = 'red'
    #     label_visualization(args, feature_matrix, label_lst, suffix, score_name)

    # # TODO: currently not plotting histogram.
    # plot_histogram(label_lst)

    # for i, e in enumerate(label_lst):
    #     if e == 'PD':
    #         label_lst[i] = 'red'
    #     elif e == 'SWEDD':
    #         label_lst[i] = 'blue'
    #     else:
    #         assert e == 'Control'
    #         label_lst[i] = 'black'

    # for i, e in enumerate(label_lst):
    #     if e < 25:
    #         label_lst[i] = 'red'
    #     elif e < 75:
    #         label_lst[i] = 'blue'
    #     else:
    #         label_lst[i] = 'blue'

if __name__ == '__main__':
    main()