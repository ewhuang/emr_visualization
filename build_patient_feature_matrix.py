### Author: Edward Huang

import numpy as np
import os
from process_loni_parkinsons import *
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize
import sys

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

def impute_missing_data(feature_matrix, master_feature_lst):
    '''
    Given the feature matrix and the column labels (master_feature_lst), impute
    the missing feature data by getting the Prosnet vectors.
    '''
    def read_prosnet_output(master_feature_lst):
        '''
        Reads the output low-dimensional vectors created by prosnet.
        '''
        vector_dct = {}
        # 450 is the last iteration number.
        f = open('./data/prosnet_data/embed_%s_450.txt' % num_dim, 'r')
        f.readline()
        for line in f:
            line = line.split()
            feature, vector = line[0], map(float, line[1:])
            # TODO: convert nodes with underscores back to spaces.
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

    vector_matrix = read_prosnet_output(master_feature_lst)
    similarity_matrix = np.abs(cosine_similarity(vector_matrix))

    similarity_matrix[similarity_matrix < sim_thresh] = 0
    np.fill_diagonal(similarity_matrix, 1)

    # Multiply the feature matrix and the similarity matrix.
    enriched_feature_matrix = np.dot(feature_matrix, similarity_matrix)

    return enriched_feature_matrix

def write_feature_matrix(updrs_dct, test_feat_mat, pros_feat_mat, test_feat_lst,
        pros_feat_lst, patient_lst, suffix):
# def write_feature_matrix(feature_matrix, master_feature_lst, patient_lst,
#     out_fname=''):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be patno's, and the second/third column should be the
    death/time event.
    '''
    assert len(test_feat_mat) == len(pros_feat_mat)
    out_fname = '%s/feature_matrix_%s.txt' % (matrix_folder, suffix)
    # if out_fname == '':
    #     if isImputation:
    #         out_fname = '%s/feature_matrix_%s_%s.txt' % (matrix_folder,
    #             norm_type, num_dim)
    #     else:
    #         out_fname = '%s/feature_matrix_%s.txt' % (matrix_folder, norm_type)
    out = open(out_fname, 'w')
    out.write('patno\tupdrs\t%s\t%s\n' % ('\t'.join(test_feat_lst),
        '\t'.join(pros_feat_lst)))
    for i, row in enumerate(test_feat_mat):
        # Write out the tests not in ProSNet.
        patno = patient_lst[i]
        out.write('%s\t%f\t%s\t' % (patno, updrs_dct[patno], '\t'.join(map(str,
            row))))
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

def get_non_prosnet_feat_tuples(updrs_dct):
    '''
    Read each of the relevant spreadsheets that are not in ProSNet.
    '''
    return (read_test_score('benton', updrs_dct), read_epworth_scale(
        updrs_dct), read_family_history(updrs_dct), read_hvlt(updrs_dct),
        read_test_score('lns', updrs_dct), read_test_score('schwab', updrs_dct),
        read_test_score('montreal', updrs_dct), read_test_score('semantic',
            updrs_dct), read_test_score('symbol', updrs_dct))

def get_prosnet_feat_tuples(updrs_dct):
    '''
    Read each of the relevant spreadsheets that are in ProSNet.
    '''
    code_dct = read_code_file()
    return (read_test_analysis('biospecimen', updrs_dct), read_test_analysis(
        'concom_medications', updrs_dct), read_clinical_diagnosis(code_dct,
        updrs_dct), read_cognitive_categorizations(updrs_dct),
        read_medical_conditions(updrs_dct), read_binary_tests('neuro',
            updrs_dct), read_binary_tests('pd_features', updrs_dct),
        read_binary_tests('rem_disorder', updrs_dct), read_demographics(
            updrs_dct), read_pd_surgery(updrs_dct), read_binary_tests(
            'medication', updrs_dct), read_mutation_file(updrs_dct))

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

def main():
    if len(sys.argv) not in [2, 4]:
        print 'Usage:python %s norm_type num_dim<optional> sim_thresh<optional>' % (
            sys.argv[0])
        exit()
    # global isImputation, matrix_folder, norm_type
    global matrix_folder, norm_type
    norm_type = sys.argv[1]
    assert norm_type in ['l1', 'l2', 'max']
    isImputation = False
    if len(sys.argv) == 4:
        isImputation = True
        global num_dim, sim_thresh
        num_dim, sim_thresh = sys.argv[2], float(sys.argv[3])
        assert num_dim.isdigit()

    updrs_dct = get_updrs_dct()
    # Test tuples are features that are not used in the ProSNet network.
    test_feat_tuples = get_non_prosnet_feat_tuples(updrs_dct)
    prosnet_feat_tuples = get_prosnet_feat_tuples(updrs_dct)

    test_feat_dct_lst, test_feat_lst = create_dct_lst(test_feat_tuples)
    pros_feat_dct_lst, pros_feat_lst = create_dct_lst(prosnet_feat_tuples)

    # Create the numpy array, and remove bad columns.
    patient_lst = updrs_dct.keys()
    test_feat_mat, test_feat_lst = build_feature_matrix(test_feat_dct_lst,
        test_feat_lst, patient_lst)
    pros_feat_mat, pros_feat_lst = build_feature_matrix(pros_feat_dct_lst,
        pros_feat_lst, patient_lst)

    matrix_folder = './data/feature_matrices'
    if not os.path.exists(matrix_folder):
        os.makedirs(matrix_folder)
    # Write out to file a unnormalized file.
    write_feature_matrix(updrs_dct, test_feat_mat, pros_feat_mat, test_feat_lst,
        pros_feat_lst, patient_lst, 'unnorm')

    # Normalize feature matrix. TODO: tune type of normalization. Before or
    # after imputation?
    test_feat_mat = normalize(test_feat_mat, norm=norm_type, axis=0)
    pros_feat_mat = normalize(pros_feat_mat, norm=norm_type, axis=0)

    # Perform either mean imputation or embedding imputation.
    if isImputation:
        pros_feat_mat = impute_missing_data(pros_feat_mat, pros_feat_lst)
        
    pros_feat_mat = normalize(pros_feat_mat, norm=norm_type, axis=0)

    # Write out normalized matrix. Imputated matrix for ProSNet file.
    if isImputation:
        suffix = '%s_%s_%s' % (norm_type, num_dim, sim_thresh)
    else:
        suffix = norm_type
    write_feature_matrix(updrs_dct, test_feat_mat, pros_feat_mat, test_feat_lst,
        pros_feat_lst, patient_lst, suffix)

    print_sparsity_info(test_feat_mat)
    print_sparsity_info(pros_feat_mat)

if __name__ == '__main__':
    main()