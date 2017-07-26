### Author: Edward Huang

import argparse
import numpy as np
import os
from process_loni_parkinsons import *
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize
# import sys

# This script creates the feature matrices to get ready for experiments.
# Writes an unnormalized feature matrix, in addition to a normalized one,
# depending on if we desire imputation using ProSNet or not.

def build_feature_matrix(feature_dct_list, master_feature_lst, patient_lst):
    '''
    Takes the feature list and a dictionary, and build a feature matrix.
    '''
    feature_matrix = []

    for patno in patient_lst:
        # Initialize the row for the patient.
        row = [0 for i in master_feature_lst]
        # Get the values from each of the feature dictionaries.
        for feature_dct in feature_dct_list:
            # A patient number might not be in feature_dct with ProSNet.
            if patno not in feature_dct:
                continue
            for (feature, feature_freq) in feature_dct[patno]:
                row[master_feature_lst.index(feature)] += feature_freq
        feature_matrix += [row]
    feature_matrix = np.array(feature_matrix)

    # Remove the bad columns.
    zero_column_indices = np.where(~feature_matrix.any(axis=0))[0]
    good_indices = [i for i in range(len(master_feature_lst)
        ) if i not in zero_column_indices]
    master_feature_lst = [e for i, e in enumerate(master_feature_lst
        ) if i in good_indices]

    return feature_matrix[:,good_indices], master_feature_lst

def impute_missing_data(feature_matrix, master_feature_lst, num_dim, sim_thresh):
    '''
    Given the feature matrix and the column labels (master_feature_lst), impute
    the missing feature data by getting the Prosnet vectors.
    '''
    def read_prosnet_output(master_feature_lst, num_dim):
        '''
        Reads the output low-dimensional vectors created by prosnet.
        '''
        vector_dct = {}
        # 450 is the last iteration number.
        f = open('./data/prosnet_data/embed_%s_1000.txt' % num_dim, 'r')
        f.readline()
        for line in f:
            line = line.split()
            feature, vector = line[0], map(float, line[1:])
            while '_' in feature:
                feature = feature.replace('_', ' ')
            assert len(vector) == int(num_dim) and feature not in vector_dct
            vector_dct[feature] = vector
        f.close()
        # Reorganize the matrix according to the order of master_feature_lst.
        vector_matrix = []
        for feature in master_feature_lst:
            vector_matrix += [vector_dct[feature]]
        return np.array(vector_matrix)

    vector_matrix = read_prosnet_output(master_feature_lst, num_dim)
    similarity_matrix = np.abs(cosine_similarity(vector_matrix))

    similarity_matrix[similarity_matrix < sim_thresh] = 0
    np.fill_diagonal(similarity_matrix, 1)

    # Multiply the feature matrix and the similarity matrix.
    enriched_feature_matrix = np.dot(feature_matrix, similarity_matrix)

    return enriched_feature_matrix

def write_feature_matrix(test_feat_mat, pros_feat_mat, test_feat_lst,
    pros_feat_lst, patient_lst, suffix):
# def write_feature_matrix(feature_matrix, master_feature_lst, patient_lst,
#     out_fname=''):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be patno's, and the second/third column should be the
    death/time event.
    '''
    assert len(test_feat_mat) == len(pros_feat_mat)
    out_fname = './data/feature_matrices/feature_matrix_%s.tsv' % suffix
    out = open(out_fname, 'w')
    out.write('patno\t%s\t%s\n' % ('\t'.join(test_feat_lst),
        '\t'.join(pros_feat_lst)))
    for i, row in enumerate(test_feat_mat):
        # Write out the tests not in ProSNet.
        patno = patient_lst[i]
        # TODO: %s and %f, depending on what the labels are.
        out.write('%s\t%s\t' % (patno, '\t'.join(map(str, row))))
        # Write out the ProSNet features.
        prosnet_features = pros_feat_mat[i]
        out.write('%s\n' % ('\t'.join(map(str, prosnet_features))))
    out.close()

def print_sparsity_info(feature_matrix):
    '''
    Given the feature matrix, print the average number of zeros for each
    patient.
    '''
    num_zeros = 0.0
    for row in feature_matrix:
        for ele in row:
            if ele == 0:
                num_zeros += 1
    print 'average number of zeros:', num_zeros / float(feature_matrix.shape[0])
    print 'feature matrix shape:', feature_matrix.shape

def get_non_prosnet_feat_tuples():
    '''
    Read each of the relevant spreadsheets that are not in ProSNet.
    '''
    return (read_test_score('benton'), read_epworth_scale(),
        read_family_history(), read_hvlt(), read_test_score('lns'),
        read_test_score('schwab'), read_test_score('montreal'),
        read_test_score('semantic'), read_test_score('symbol'))

def read_snp_fisher_file(fname, snp_set):
    '''
    Given a Fisher's test file on finding significant SNPs, get the list of
    SNPs.
    '''
    f = open(fname, 'r')
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        snp_set.add(line[0])
    f.close

def get_mutation_dct():
    # First, get the SNPs deemed to be significantly enriched in PD patients.
    snp_set = set([])
    for fname in ('snp_fisher_test_wes_ppmi', 'snp_fisher_test_wes_hard_ignore'):
        read_snp_fisher_file('./data/ppmi/snp_files/%s.tsv' % fname, snp_set)
        # break # TODO
    # Next, for each SNP, get the patients that have the SNP.
    with open('./data/ppmi/snp_files/patno_snp_dct_wes.json', 'r') as fp:
        patno_snp_dct = json.load(fp)
    fp.close()

    # Keep only the SNPs that were deemed significant.
    feature_set = set([])
    for patno in patno_snp_dct:
        patno_snp_set = snp_set.intersection(patno_snp_dct[patno])
        # patno_snp_dct[patno] = snp_set.intersection(patno_snp_dct[patno])
        feature_set = feature_set.union(patno_snp_set)
        patno_snp_dct[patno] = [(snp, 1) for snp in patno_snp_set]

    return patno_snp_dct, feature_set

def get_prosnet_feat_tuples():
    '''
    Read each of the relevant spreadsheets that are in ProSNet.
    '''
    code_dct = read_code_file()
    # Make sure SNP mutations comes before biospecimen, since they share feature 'APP'.
    # return (read_snp_mutations(), read_test_analysis('biospecimen'),
    return (read_test_analysis('biospecimen'), read_test_analysis('concom_medications'),
        read_clinical_diagnosis(code_dct), read_cognitive_categorizations(),
        read_medical_conditions(), read_binary_tests('neuro'),
        read_binary_tests('pd_features'), read_binary_tests('rem_disorder'),
        read_demographics(), read_pd_surgery(), read_binary_tests('medication'),
        get_mutation_dct())

def create_dct_lst(feat_tuples):
    '''
    Given a list of (dictionary, feature list) tuples, make a master list of the
    dictionaries, and compile a list of the unique features over all tuples.
    '''
    feature_dct_list, master_feature_lst = [], []
    for feature_dct, feature_list in feat_tuples:
        feature_dct_list += [feature_dct]
        # Update the master feature list.
        for feature in feature_list:
            # feature = '_'.join(feature.split())
            if feature not in master_feature_lst:
                master_feature_lst += [feature]
            else:
                # This means that some features are repeated across different
                # spreadsheets. TODO.
                print feature
    return feature_dct_list, master_feature_lst

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--norm_type', help='Normalization method.',
        required=True, choices=['max', 'l1', 'l2'])
    parser.add_argument('-d', '--num_dim', help='Number of ProSNet dimensions',
        type=int)
    parser.add_argument('-s', '--sim_thresh', type=float,
        help='Threshold for cosine similarity between ProSNet vectors')
    parser.add_argument('-w', '--where_norm', choices=['before', 'after', 'both'],
        help='Where to normalize with respect to ProSNet imputation.')
    parser.add_argument('-l', '--label_type', choices=['updrs', 'status', 'tau'],
        help='Use the patients that have this label.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    # TODO: leaving labels out of the feature matrix.
    label_type = args.label_type
    if label_type == 'updrs':
        label_dct = get_updrs_dct()[0]
    # elif label_type == 'status': # PD, SWEDD, Healthy Control.
    #     non_pd_patients = []
    #     label_dct = read_patient_status()
    #     # TODO: Take out the non-PD patients.
    #     for patno in label_dct:
    #         if label_dct[patno] != 'PD':
    #             non_pd_patients += [patno]
    #     for patno in non_pd_patients:
    #         del label_dct[patno]
    elif label_type == 'tau':
        label_dct = read_total_tau()
    # Test tuples are features that are not used in the ProSNet network.
    test_feat_tuples = get_non_prosnet_feat_tuples()
    prosnet_feat_tuples = get_prosnet_feat_tuples()

    test_feat_dct_lst, test_feat_lst = create_dct_lst(test_feat_tuples)
    pros_feat_dct_lst, pros_feat_lst = create_dct_lst(prosnet_feat_tuples)

    # Create the numpy array, and remove bad columns.
    # TODO: Using intersection of PD patients and those with UPDRS scores.
    # pd_dct = read_patient_status()
    # pd_keys = [key for key in pd_dct if pd_dct[key] == 'PD']

    # patient_lst = list(set(label_dct.keys()).intersection(pd_keys))
    patient_lst = label_dct.keys()
    test_feat_mat, test_feat_lst = build_feature_matrix(test_feat_dct_lst,
        test_feat_lst, patient_lst)
    pros_feat_mat, pros_feat_lst = build_feature_matrix(pros_feat_dct_lst,
        pros_feat_lst, patient_lst)

    matrix_folder = './data/feature_matrices'
    if not os.path.exists(matrix_folder):
        os.makedirs(matrix_folder)
    # Write out to file a unnormalized file.
    write_feature_matrix(test_feat_mat, pros_feat_mat, test_feat_lst,
        pros_feat_lst, patient_lst, 'unnorm')

    # Normalize feature matrix. TODO: tune type of normalization. Before or
    # after imputation?
    test_feat_mat = normalize(test_feat_mat, norm=args.norm_type, axis=0)
    if args.where_norm in ['before', 'both']:
        pros_feat_mat = normalize(pros_feat_mat, norm=args.norm_type, axis=0)

    # Perform either mean imputation or embedding imputation.
    if args.num_dim != None:
        pros_feat_mat = impute_missing_data(pros_feat_mat, pros_feat_lst,
            args.num_dim, float(args.sim_thresh))
        suffix = '%s_%s_%s_%s' % (args.norm_type, args.num_dim, args.sim_thresh,
            args.where_norm)
    else:
        suffix = args.norm_type

    if args.where_norm in ['after', 'both', None]:
        pros_feat_mat = normalize(pros_feat_mat, norm=args.norm_type, axis=0)

    write_feature_matrix(test_feat_mat, pros_feat_mat, test_feat_lst,
        pros_feat_lst, patient_lst, suffix)

    # print_sparsity_info(test_feat_mat)
    # print_sparsity_info(pros_feat_mat)

if __name__ == '__main__':
    main()