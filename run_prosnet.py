### Author: Edward Huang

import argparse
import json
import os
from process_loni_parkinsons import *
import string
import subprocess
import sys

### This script builds the network to be used as input for ProSNet.

# Records the set of nodes already written to file. Also unique number of edges.
global_node_set, global_edge_set, num_edge_types = set([]), set([]), 0

def get_ppi_edge_set():
    '''
    Get the protein-protein edge set from a PPI network.
    '''
    ppi_edge_set = set([])
    f = open('./data/InBio-Map_Symbol.sif', 'r')
    for line in f:
        node_a, node_b = line.split()
        ppi_edge_set.add((node_a, node_b))
    f.close()
    return ppi_edge_set

def get_stitch_edge_set(fname):
    '''
    Gets the drug-protein interactions obtained from the STITCH database.
    Return each tuple as (protein, drug).
    '''
    edge_set = set([])
    f = open('./data/api_calls/%s.txt' % fname, 'r')
    for line in f:
        node, drug = line.strip().split('\t')
        edge_set.add((node, drug))
    f.close()
    return edge_set

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

def write_files(node_out, edge_out, edge_set, node_type_a, node_type_b):
    '''
    Write the edges out to file. edge_label is just a letter to differentiate
    amongst the different types of edges. We write nodes out as integers, so
    map_out contains a word on each line corresponding to these integers.
    '''
    global global_node_set, global_edge_set, num_edge_types
    # The order of node type a, b must match the edge order in edge_set.
    node_type_tup = (node_type_a, node_type_b)

    for edge in edge_set:
        # Check if edge has already been written.
        if edge in global_edge_set or edge[::-1] in global_edge_set:
            continue
        global_edge_set.add(edge)

        # Write out the edge.
        for i, node in enumerate(edge):
            # Split node and put it back together with underscore.
            node = '_'.join(node.split())
            # Write out the node if it hasn't appeared yet.
            if node not in global_node_set:
                global_node_set.add(node)
                node_out.write('%s\t%s\n' % (node, node_type_tup[i]))
            # Write out the edge.
            edge_out.write('%s\t' % node)
        # Edge weights are all = 1. Map the edge type to a letter.
        edge_label = string.ascii_lowercase[num_edge_types]
        edge_out.write('1\t%s\n' % edge_label)
        # Write the edge backwards, to make it undirected.
        edge_out.write('%s\t%s\t1\t%s\n' % (edge[1], edge[0], edge_label))
        assert edge[1] != edge[0]
    num_edge_types += 1

def get_coocc_edge_set(patient_dct_a, patient_dct_b):
    '''
    Get the symptom-herb relations from the co-occurrence counts of the patient
    records.
    '''
    coocc_edge_set = set([])
    # Get the intersecting set of patients in both dictionaries.
    patient_set = set(patient_dct_a.keys()).intersection(patient_dct_b.keys())
    for inhospital_id in patient_set:
        # for (node_a, node_a_freq) in patient_dct_a[inhospital_id]:
        for node_a in patient_dct_a[inhospital_id]:
            # for (node_b, node_b_freq) in patient_dct_b[inhospital_id]:
            for node_b in patient_dct_b[inhospital_id]:
                # We skip self edges, but keep other same-type edges. TODO.
                if node_a != node_b:
                    coocc_edge_set.add((node_a, node_b))
    return coocc_edge_set

def run_prosnet():
    os.chdir('./prosnet/model')
    network_folder = '../../data/prosnet_data'
    
    # Sort out the file names.
    if args.excl_feat == None:
        node_fname = 'prosnet_node_list'
        edge_fname = 'prosnet_edge_list'
    else:
        node_fname = 'prosnet_node_list_no_%s' % args.excl_feat
        edge_fname = 'prosnet_edge_list_no_%s' % args.excl_feat

    command = ('./embed -node "%s/%s.txt" -link "%s/%s.txt" '
        '-binary 0 -size %s -negative 5 -samples 1 '
        '-iters 1001 -threads 24 -model 2 -depth 10 -restart 0.8 '
        '-edge_type_num %d -rwr_ppi 1 -rwr_seq 1 -train_mode 2' % (
            network_folder, node_fname, network_folder, edge_fname,
            args.num_dim, num_edge_types))
    print command
    subprocess.call(command, shell=True)

def get_attributes(patient_dct_lst):
    '''
    Given a patient dictionary mapping patients to (feature, feature_freq)
    tuples, instead return a dictionary mapping patients to lists of the
    features.
    '''
    new_dct = {}
    for patient_dct, feature_set in patient_dct_lst:
        # Go through each patient dictionary.
        for patno in patient_dct:
            tuple_lst = patient_dct[patno]
            if patno not in new_dct:
                new_dct[patno] = []
            new_dct[patno] += [pair[0] for pair in tuple_lst]
    return new_dct

def get_mutation_dct():
    # First, get the SNPs deemed to be significantly enriched in PD patients.
    snp_set = set([])
    for fname in ('snp_fisher_test_wes_ppmi', 'snp_fisher_test_wes_hard_ignore'):
        read_snp_fisher_file('./data/ppmi/snp_files/%s.tsv' % fname, snp_set)
    # Next, for each SNP, get the patients that have the SNP.
    with open('./data/ppmi/snp_files/patno_snp_dct_wes.json', 'r') as fp:
        patno_snp_dct = json.load(fp)
    fp.close()

    # Keep only the SNPs that were deemed significant.
    for patno in patno_snp_dct:
        patno_snp_dct[patno] = snp_set.intersection(patno_snp_dct[patno])

    return patno_snp_dct

def get_spreadsheet_results():
    '''
    Returns a dictionary mapping some key, denoting the spreadsheet, to the 
    dictionary of that spreadsheet's results. Dictionary of dictionaries.
    '''
    f_tuples = []
    # Medical tests.
    if args.excl_feat == 'biospecimen':
        test_dct = get_attributes([read_binary_tests('neuro'),
            read_binary_tests('pd_features'), read_cognitive_categorizations(),
            read_pd_surgery()])
    else:
        test_dct = get_attributes([read_test_analysis('biospecimen'),
            read_binary_tests('neuro'), read_binary_tests('pd_features'),
            read_cognitive_categorizations(), read_pd_surgery()])
    f_tuples += [('t', test_dct)]
    # Symptoms.
    if args.excl_feat != 'symptoms':
        code_dct = read_code_file()
        symp_dct = get_attributes([read_clinical_diagnosis(code_dct),
            read_medical_conditions(), read_binary_tests('rem_disorder')])
        f_tuples += [('s', symp_dct)]
    # Demographics.
    demo_dct = get_attributes([read_demographics()])
    f_tuples += [('m', demo_dct)]

    # Drugs.
    drug_dct = get_attributes([read_test_analysis('concom_medications'),
        read_binary_tests('medication')])
    f_tuples += [('d', drug_dct)]

    # Gene mutations.
    mutation_dct = get_mutation_dct()
    f_tuples += [('g', mutation_dct)]

    return f_tuples

def get_symptom_snp_edge_set():
    '''
    Reads the symptom-SNP dictionary as computed by Fisher's exact test.
    '''
    snp_edge_set = set([])
    f = open('./data/ppmi/snp_files/snp_symptom_enrichments.txt', 'r')
    for line in f:
        line = line.split()
        snp, symptom = line[0], line[1]
        snp_edge_set.add((snp, symptom))
    f.close()
    return snp_edge_set

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_dim', type=int, required=True,
        help='Number of ProSNet dimensions.')
    parser.add_argument('-e', '--excl_feat', choices=['biospecimen', 'symptoms'],
        help='Feature types to exclude.')
    args = parser.parse_args()

def main():
    parse_args()

    input_folder = './data/prosnet_data'
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    if args.excl_feat == None:
        node_out = open('%s/prosnet_node_list.txt' % input_folder, 'w')
        edge_out = open('%s/prosnet_edge_list.txt' % input_folder, 'w')
    else:
        node_out = open('%s/prosnet_node_list_no_%s.txt' % (input_folder,
            args.excl_feat), 'w')
        edge_out = open('%s/prosnet_edge_list_no_%s.txt' % (input_folder,
            args.excl_feat), 'w')

    ppi_edge_set = get_ppi_edge_set()
    write_files(node_out, edge_out, ppi_edge_set, 'p', 'p')

    protein_drug_edge_set = get_stitch_edge_set('stitch_protein_drug')
    write_files(node_out, edge_out, protein_drug_edge_set, 'p', 'd')

    f_tuples = get_spreadsheet_results()
    # Loop through every pair of node types.
    for i in range(len(f_tuples)):
        node_type_a, patient_dct_a = f_tuples[i]
        for j in range(i, len(f_tuples)):
            node_type_b, patient_dct_b = f_tuples[j]
            # Get the co-occurrence edge set.
            edge_set = get_coocc_edge_set(patient_dct_a, patient_dct_b)

            if node_type_a == 'd' and node_type_b == 'd':
                edge_set = edge_set.union(get_stitch_edge_set('stitch_drug_drug'))
            elif node_type_a == 's' and node_type_b == 'g':
                edge_set = edge_set.union(get_symptom_snp_edge_set())

            # Write the edges out to file.
            write_files(node_out, edge_out, edge_set, node_type_a, node_type_b)

    edge_out.close()
    node_out.close()
    exit()

    # Run prosnet. Outputs the low-dimensional vectors into files.
    run_prosnet()

if __name__ == '__main__':
    main()