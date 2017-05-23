### Author: Edward Huang

import os
import process_loni_parkinsons
import string
import subprocess
import sys

### This script builds the network to be used as input for ProSNet.

# Records the set of nodes already written to file. Also unique number of edges.
global_node_set, global_edge_set, num_edge_types = set([]), set([]), 0

def get_entrez_to_hgnc_dct():
    '''
    Gets mappings from HGNC ID's to Entrez ID's.
    '''
    entrez_to_hgnc_dct = {}
    f = open('./data/hgnc_to_entrez.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip().split('\t')
        if len(line) != 2:
            continue
        hgnc_id, entrez_id = line
        assert entrez_id not in entrez_to_hgnc_dct
        entrez_to_hgnc_dct[entrez_id] = hgnc_id
    f.close()    
    return entrez_to_hgnc_dct

def get_ppi_edge_set():
    '''
    Get the protein-protein edge set from a PPI network.
    '''
    entrez_to_hgnc_dct = get_entrez_to_hgnc_dct()

    ppi_edge_set = set([])
    # Gene ID's in this PPI network are Entrez ID's.
    f = open('./data/HumanNet.v1.benchmark.txt', 'r')
    for line in f:
        line = line.split()
        assert len(line) == 2
        node_a, node_b = line
        # Skip if no HGNC analogues.
        if node_a not in entrez_to_hgnc_dct or node_b not in entrez_to_hgnc_dct:
            continue
        # Translate the Entrez ID to HGNC protein.
        ppi_edge_set.add((entrez_to_hgnc_dct[node_a],
            entrez_to_hgnc_dct[node_b]))
    f.close()
    return ppi_edge_set

def get_protein_drug_edge_set():
    '''
    Gets the drug-protein interactions obtained from the STITCH database.
    Return each tuple as (protein, drug).
    '''
    protein_drug_edge_set = set([])
    api_dir = './data/api_calls'
    fname_lst = os.listdir(api_dir)
    for fname in fname_lst:
        if 'interactions' not in fname:
            continue
        f = open('%s/%s' % (api_dir, fname), 'r')
        for line in f:
            node_a, node_b = line.strip().split('\t')
            if node_a.isupper():
                assert node_b.islower()
                protein_drug_edge_set.add((node_a, node_b))
            else:
                assert node_b.isupper()
                protein_drug_edge_set.add((node_b, node_a))
        f.close()
        break
    return protein_drug_edge_set

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
        # Edge weights are all = 1. Map the edge type to a letter. TODO.
        edge_label = string.ascii_lowercase[num_edge_types]
        # edge_label = num_edge_types + 1
        edge_out.write('1\t%s\n' % edge_label)
    num_edge_types += 1

def run_prosnet():
    # os.chdir('../simons_mouse/Sheng/prosnet/model')
    os.chdir('./prosnet/model')
    # network_folder = '../../../../cancer_survival_tcm/data/prosnet_data'
    network_folder = '../../data/prosnet_data'
    command = ('./embed -node "%s/prosnet_node_list.txt" -link "%s/prosnet_'
        'edge_list.txt" -binary 0 -size %s -negative 5 -samples 1 '
        '-iters 500 -threads 12 -model 2 -depth 10 -restart 0.8 '
        '-edge_type_num %d -rwr_ppi 1 -rwr_seq 1 -train_mode 2' % (
            network_folder, network_folder, num_dim, num_edge_types))
    print command
    subprocess.call(command, shell=True)

def get_spreadsheet_results():
    '''
    Returns a dictionary mapping some key, denoting the spreadsheet, to the 
    dictionary of that spreadsheet's results. Dictionary of dictionaries.
    '''
    line_orientation_dct = process_loni_parkinsons.read_test_score('benton')

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s num_dim' % sys.argv[0]
        exit()
    global num_dim
    num_dim = sys.argv[1]
    assert num_dim.isdigit()

    input_folder = './data/prosnet_data'
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)
    node_out = open('%s/prosnet_node_list.txt' % input_folder, 'w')
    edge_out = open('%s/prosnet_edge_list.txt' % input_folder, 'w')

    ppi_edge_set = get_ppi_edge_set()
    write_files(node_out, edge_out, ppi_edge_set, 'p', 'p')
    protein_drug_edge_set = get_protein_drug_edge_set()
    write_files(node_out, edge_out, protein_drug_edge_set, 'p', 'd')

    edge_out.close()
    node_out.close()

    get_spreadsheet_results()
    exit()

    # Run prosnet. Outputs the low-dimensional vectors into files.
    run_prosnet()

if __name__ == '__main__':
    main()