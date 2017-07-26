### Author: Edward Huang

import argparse
from collections import Counter
import json
import matplotlib
import numpy as np
import os
from process_loni_parkinsons import *
from scipy.stats import fisher_exact, ttest_rel
from sklearn.cluster import DBSCAN
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

np.random.seed(seed=930519)
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
         (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
         (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
         (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
         (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i, (r, g, b) in enumerate(tableau20):
    tableau20[i] = (r / 255., g / 255., b / 255.)

# These are the "Tableau 20" colors as RGB.
traffic_light = [(216, 37, 38), (255, 193, 86), (159, 205, 153)]
for i in range(len(traffic_light)):
    r, g, b = traffic_light[i]
    traffic_light[i] = (r / 255., g / 255., b / 255.)

def generate_directories():
    for results_folder in ('./results', './results/emr_points',
        './results/baseline_accuracy', './results/prosnet_accuracy',
        './results/biomarker_enrichments', './results/color_biomarker_mappings'):
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)

def read_feature_matrix(suffix):
    '''
    Reads the feature matrix of the patient data. Takes an optional argument in
    the form number of dimensions.
    '''
    # feature_matrix, label_lst, master_feature_list = [], [], []
    feature_matrix, patient_lst = [], []

    if args.excl_feat == 'None' or suffix == 'unnorm':
        fname = './data/feature_matrices/feature_matrix_%s.tsv' % suffix
    elif args.excl_feat == 'biospecimen':
        fname = './data/feature_matrices/feature_matrix_%s_no_biospecimen.tsv' % suffix

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
    return np.array(feature_matrix), header[1:], patient_lst

# def plot_histogram(label_lst):
#     bins = range(200)
#     matplotlib.pyplot.hist(label_lst, bins, alpha=0.25)

#     plt.show()
#     pylab.savefig('./results/emr_points/updrs_histogram.png')
#     plt.close()

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--norm_type', required=True, choices=['l1', 'l2', 'max'],
        help='Normalization method.')
    parser.add_argument('-d', '--num_dim', help='number of ProSNet dimensions')
    parser.add_argument('-s', '--sim_thresh', help='similarity threshold for ProSNet imputation')
    parser.add_argument('-w', '--where_norm', help='where to normalize: before or after (or both) ProSNet imputation')
    parser.add_argument('-p', '--n_pca_comp', help='number of PCA components',
        required=True, type=int)
    parser.add_argument('-t', '--tsne_init', help='t-SNE initialization method',
        required=True, choices=['pca', 'random'])
    parser.add_argument('-r', '--l_rate', help='learning rate of t-SNE',
        required=True, type=int)
    parser.add_argument('-e', '--excl_feat', choices=['biospecimen'],
        help='Feature types to exclude.')
    # parser.add_argument('-k', '--knn', help='number of nearest neighbors')
    args = parser.parse_args()

# def knn_predict(args, fname, X, y):
#     '''
#     Gets the classification accuracy on a 3-fold CV with a KNN classifier.
#     Writes out the score to file if this is a baseline run.
#     '''
#     # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1,
#     #     random_state=9305)
#     # TODO: make n_neighbors be a command line option.
#     knn = KNeighborsClassifier(n_neighbors=30)
#     # knn.fit(X_train, y_train)
#     # pred = knn.predict(X_test)

#     # score = accuracy_score(y_test, pred)
#     # print score
#     # TODO: 10-fold CV.
#     # score = cross_val_score(knn, X, y, cv=10, scoring='accuracy').mean()
#     score_lst = cross_val_score(knn, X, y, cv=10, scoring='accuracy')

#     # Write the score out to file if this is a baseline run.
#     if args.num_dim == None:
#         folder = './results/baseline_accuracy'
#     else:
#         folder = './results/prosnet_accuracy'
#     out = open('%s/%s.txt' % (folder, fname), 'w')
#     out.write('%s' % '\n'.join(map(str,score_lst)))
#     out.close()
#     return score_lst

def get_symptom_drug_set():
    '''
    Gets the set of symptom and drug names in the dataset.
    '''
    symptom_drug_set = set([])
    code_dct = read_code_file()
    symptom_drug_set = symptom_drug_set.union(read_clinical_diagnosis(code_dct)[1])
    symptom_drug_set = symptom_drug_set.union(read_medical_conditions()[1])
    symptom_drug_set = symptom_drug_set.union(read_binary_tests('rem_disorder')[1])
    symptom_drug_set = symptom_drug_set.union(read_test_analysis('concom_medications')[1])
    symptom_drug_set = symptom_drug_set.union(read_binary_tests('medication')[1])
    return symptom_drug_set

def get_biospecimen_set():
    '''
    Gets the set of biospecimen names from the PPMI dataset.
    '''
    biospecimen_set = set([])
    biospecimen_set = biospecimen_set.union(read_test_analysis('biospecimen')[1])
    return biospecimen_set

def compute_cluster_enrichment(clus_feat_matrix, non_clus_feat_matrix, feat_idx):
    '''
    Given a cluster of patients and a feature, determine whether the cluster is
    enriched in the feature.
    '''
    # Find the number of patients in the cluster with the feature.

    clus_and_feat = np.count_nonzero(clus_feat_matrix[:,feat_idx])
    clus_no_feat = len(clus_feat_matrix) - clus_and_feat
    non_clus_and_feat = np.count_nonzero(non_clus_feat_matrix[:,feat_idx])
    non_clus_no_feat = len(non_clus_feat_matrix) - non_clus_and_feat
    f_table = [[clus_and_feat, clus_no_feat], [non_clus_and_feat, non_clus_no_feat]]
    o_r, p_value = fisher_exact(f_table, alternative='greater')
    return p_value, f_table

def plot_dbscan_clusters(feature_matrix, labels):
    '''
    Plots the DBSCAN clustering results.
    '''

    # # # Map each cluster to a color in the tableau.
    # colored_clusters = [-1]
    # tableau_idx = 0
    # for label in labels:
    #     if label not in colored_clusters:
    #         colored_clusters += [label]
    #         tableau_idx += 1
    # color_map_dct = {}
    # c = Counter(labels)
    # tableau_idx = 0
    # for clus_idx, count in c.most_common(21):
    #     if clus_idx == -1:
    #         continue
    #     color_map_dct[clus_idx] = tableau20[tableau_idx]
    #     tableau_idx += 1

    # Convert the labels to colors. Only top ten clusters get colors. Everything else white.
    color_lst = [(0, 0, 0)] * len(labels)
    for i, label in enumerate(labels):
        if label != -1:
            # color_lst[i] = color_map_dct[label]
            color_lst[i] = tableau20[label]
        # # if label == -1:
        # else:
        #     color_lst[i] = (0., 0., 0.)

    plot_embeddings(feature_matrix, color_lst)
    return tableau20

def feature_analysis(feature_matrix, fname):
    # TODO: eps was 0.3.
    labels = DBSCAN(eps=1, min_samples=10, n_jobs=-1).fit_predict(feature_matrix)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)

    # Read the base feature matrix for feature analysis.
    base_feature_matrix, feature_lst, patient_lst = read_feature_matrix('unnorm')

    if args.excl_feat == None:
        enrichment_set = get_symptom_drug_set()
    elif args.excl_feat == 'biospecimen':
        enrichment_set = get_biospecimen_set()

    print enrichment_set

    clus_feat_enrich_dct = {}
    for i in range(n_clusters):
        # Initialize each cluster's set of enriched features.
        clus_feat_enrich_dct[i] = []
        # Get the cluster indices in the feature matrix.
        clus_idx_lst = [j for j, label in enumerate(labels) if label == i]
        non_clus_idx_lst = [j for j, label in enumerate(labels) if label != i]
        # Get the current cluster for the feature matrix.
        clus_feat_matrix = base_feature_matrix[clus_idx_lst]
        non_clus_feat_matrix = base_feature_matrix[non_clus_idx_lst]

        # Perform the enrichment analysis for each feature.
        for feat_idx, feature in enumerate(feature_lst):
            if feature not in enrichment_set:
                continue
            print feature
            p_value, f_table = compute_cluster_enrichment(clus_feat_matrix, non_clus_feat_matrix, feat_idx)
            # TODO: Currently writing out all things.
            # if p_value < 0.01:
            clus_feat_enrich_dct[i] += [(feature, f_table, p_value)]
    # return clus_feat_enrich_dct
    if args.excl_feat == None:
        with open('./results/biomarker_enrichments/%s.json' % fname, 'w') as fp:
            json.dump(clus_feat_enrich_dct, fp)
        out = open('./results/color_biomarker_mappings/%s.txt' % fname, 'w')
    elif args.excl_feat == 'biospecimen':
        with open('./results/biomarker_enrichments/%s_biospecimen.json' % fname, 'w') as fp:
            json.dump(clus_feat_enrich_dct, fp)
        out = open('./results/color_biomarker_mappings/%s_biospecimen.txt' % fname, 'w')
    # TODO: exiting here for now.
    exit()
    tableau20 = plot_dbscan_clusters(feature_matrix, labels)
    # TODO: figure out top symptom for each cluster.
    # for clus_idx, (r, g, b) in enumerate(tableau20):
    for clus_idx in clus_feat_enrich_dct:
        # color = tuple([rgb * 255 for rgb in color_map_dct[clus_idx]])
        color = tuple([rgb * 255 for rgb in tableau20[clus_idx]])
        biomarker_p_lst = clus_feat_enrich_dct[clus_idx]
        # Sort the list of tuples by the p-values.
        curr_p_lst = sorted(biomarker_p_lst, key=lambda x:x[2])[:10]
        biomarkers = [tup[0] for tup in curr_p_lst]
        out.write('%s\t%d\t%s\n' % (color, list(labels).count(clus_idx), '\t'.join(biomarkers)))
    out.close()

def label_visualization(feature_matrix, label_lst, suffix, test_name):
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

    # Perform feature analysis on the visualization.
    print 'currently no feature analysis...'
    feature_analysis(feature_matrix, fname)

    # accuracy_lst = knn_predict(args, fname, feature_matrix, label_lst)
    # accuracy = accuracy_lst.mean()
    # print accuracy
    plot_updrs(feature_matrix, label_lst)

    # print 'plotting...'

    # # Always plot baseline plots.
    # if args.num_dim == None: # Title for baseline is just the accuracy score.
    #     title = 'baseline accuracy (%g)' % accuracy
    #     plot_embeddings(feature_matrix, label_lst, fname, title)
    # # Plot a ProSNet plot if its accuracy is better than its corresponding baseline.
    # else:
    #     baseline_score_lst = get_baseline_accuracy(args, test_name)
    #     # TODO: plotting all plots.
    #     # if accuracy > baseline_score:
    #     title = '%g vs. baseline (%g)' % (accuracy, np.mean(baseline_score_lst))
    #     plot_embeddings(feature_matrix, label_lst, fname, title)
    #     print ttest_rel(accuracy_lst, baseline_score_lst)

# def get_baseline_accuracy(args, test_name):
#     '''
#     Get the classification accuracy of the baseline method.
#     '''
#     fname = '%s_%s_%s_%s_%s' % (args.norm_type, args.n_pca_comp,
#         args.tsne_init, test_name, args.l_rate)
#     score_lst = []
#     f = open('./results/baseline_accuracy/%s.txt' % fname, 'r')
#     for line in f:
#         score_lst += [float(line.strip())]
#     f.close()
#     return score_lst

def plot_embeddings(feature_matrix, label_lst):
    # x_points = [point[0] for point in feature_matrix]
    # y_points = [point[1] for point in feature_matrix]
    plt.figure(figsize=(10, 7.5))

    marker_lst = ('.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', '8',
        's', 'p', 'h', 'H', '+', 'x', 'd', '|', '_')
    # Match a color to each cluster.
    color_lst = tableau20[:len(set(label_lst)) - 1]
    assert len(color_lst) <= len(marker_lst)
    # Set the enrichments here.
    if args.num_dim == None:
        sympt_lst = ('Menopausal depression\nHypothyroidism',
            'REM sleep disorder\nEnlarged prostate', 'Parkinsonism\nHypertension',
            'Corrective lens user\nHypertension', 'Pregnancy\nPolycystic ovaries',
            'Hormone replacement therapy\nCaesarean section', 'Parkinsonism\nParkinson\'s disease',
            'Hypertension\nHypercholesterolaemia', 'Refraction disorder\nParkinsonism',
            'Apathy\nCoronary artery disease')
    else:
        sympt_lst = ('Extremity pain\nImpaired glucose tolerance',
            'Animal allergy\nSeasonal allergy', 'Parkinsonism\nHead injury',
            'REM sleep disorder\nParkinson\'s disease', 'Parkinsonism\nHypertension',
            'Parkinsonism\nBradykinesia', 'Corrective lens user\nHypertension',
            'Myopia\nHypertension', 'Melanocytic nevus\nAsthma',
            'Hypotension\nCholelithiasis', 'Parkinsonism\nParkinson\'s disease',
            'Parkinsonism\nBack injury', 'Polycystic ovaries\nHPV positive',
            'Atrial septal defect\nAcid reflux', 'Hypothyroidism\nMitral valve prolapse',
            'Uterine prolapse\nPneumonitis', 'Hypertension\nLeft bundle branch block',
            'Ear canal injury\nCystitis escherichia')

    for clr_idx in range(len(color_lst)):
        color, marker, label = color_lst[clr_idx], marker_lst[clr_idx], sympt_lst[clr_idx]
        # Get the points corresponding to this color.
        x_pts, y_pts = zip(*[point for i, point in enumerate(feature_matrix)
            if label_lst[i] == color])

        plt.scatter(x=x_pts, y=y_pts, c=color, s=50, edgecolors='none', marker=marker,
            label=label)

    # Plot noise points.
    x_pts, y_pts = zip(*[point for i, point in enumerate(feature_matrix)
        if label_lst[i] == (0, 0, 0)])
    plt.scatter(x=x_pts, y=y_pts, c=(0, 0, 0), s=20, edgecolors='none', marker='D',
        label='Noise points')

    # # other_x_points = [point[0] for i, point in enumerate(feature_matrix) if label_lst[i] != red]
    # # other_y_points = [point[1] for i, point in enumerate(feature_matrix) if label_lst[i] != red]
    # # Plot resulting feature matrix.
    # # edge_color = (199/255., 199/255., 199/255.)
    # color_lst = [label for label in label_lst if label != red]
    # color_to_marker = {}
    # for color in color_lst:
    #     if color not in color_to_marker:
    #         color_to_marker[color] = marker_lst.pop(0)
    # # plt.legend(loc='lower right')
    # print color_to_marker
    # # marker_lst = [color_to_marker[label] for label in label_lst if label != red]

    # for c, _x, _y in zip(color_lst, other_x_points, other_y_points):
    #     plt.scatter(x=_x, y=_y, c=c,
    #         s=20, edgecolors='none', marker=color_to_marker[c])

    if args.num_dim == None:
        plt.xlim(-16,15.5)
        plt.ylim(-10.5, 13)
    else:
        plt.xlim(-17, 12)
        plt.ylim(-9.5,13.5)

    # Legend outside of plot.
    plt.legend(loc='upper center', scatterpoints=1, bbox_to_anchor=(0.5, -0.05), ncol=3,
        frameon=False, fontsize=20)

    plt.axis('off')
    plt.tight_layout()
    plt.show()

    if args.num_dim == None:
        fname = 'baseline_dbscan_labels_with_plot'
    else:
        fname = 'prosnet_dbscan_labels_with_plot'
    pylab.savefig('./results/emr_points/%s.pdf' % fname, bbox_inches='tight')
    plt.close()

def plot_updrs(feature_matrix, label_lst):
    '''
    Plot the 2D EMR points, colored by UPDRS sums.
    '''
    plt.figure(figsize=(10, 7.5))

    # traffic_light is red, green, yellow. Plot in reverse order.
    marker_lst = ['D', 'o', 's']
    color_lst = traffic_light[::-1]
    severity_lst = ['Minor', 'Moderate', 'Severe']

    for clr_idx in range(len(color_lst)):
        color, marker, label = color_lst[clr_idx], marker_lst[clr_idx], severity_lst[clr_idx]
        # Get the points corresponding to this color.
        x_pts, y_pts = zip(*[point for i, point in enumerate(feature_matrix)
            if label_lst[i] == color])
        # Red gets alpha=1
        if clr_idx == 2:
            alpha = 1
        else:
            alpha = 0.9
        plt.scatter(x=x_pts, y=y_pts, c=color, s=30, edgecolors='none',
            marker=marker, alpha=alpha, label=label)

    if args.num_dim == None:
        plt.xlim(-16,15.5)
        plt.ylim(-10.5, 13)
    else:
        plt.xlim(-17, 12)
        plt.ylim(-9.5,13.5)

    plt.legend(loc='upper left', scatterpoints=1, fontsize=20, framealpha=0.5)

    plt.axis('off')
    plt.tight_layout()
    plt.show()

    if args.num_dim == None:
        fname = 'updrs_baseline'
    else:
        fname = 'updrs_prosnet'
    pylab.savefig('./results/emr_points/%s.pdf' % fname, bbox_inches='tight')
    plt.close()

def main():
    generate_directories()
    parse_args()

    # Process the filename for each run.
    suffix = args.norm_type # Baseline only contains normalization method.
    assert suffix in ['l1', 'l2', 'max']
    if args.num_dim != None:
        assert args.num_dim.isdigit()
        suffix += '_%s_%s_%s' % (args.num_dim, args.sim_thresh, args.where_norm)

    # feature_matrix, label_lst, feature_list = read_feature_matrix(suffix)
    feature_matrix, feature_lst, patient_lst = read_feature_matrix(suffix)

    label_dct, score_names = get_updrs_dct()
    # TODO: assertion statement. Currently using only PD patients.
    assert set(label_dct.keys()).intersection(patient_lst) == set(patient_lst)

    # Get the label list from the label dictionary.
    # 1. Sum of UPDRS scores.
    label_lst = []
    for patno in patient_lst:
        score = sum(label_dct[patno].values())
        label_lst += [score]
    for i, e in enumerate(label_lst):
        # Convert UPDRS scores to colors based on their range.
        if e < 35:
            label_lst[i] = traffic_light[2]
        elif e < 70:
            label_lst[i] = traffic_light[1]
        else:
            label_lst[i] = traffic_light[0]
    label_visualization(feature_matrix, label_lst, suffix, 'sum')

if __name__ == '__main__':
    main()