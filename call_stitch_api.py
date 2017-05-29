### Author: Edward Huang

from csv import reader
import json
import os
import urllib

# This file writes the file of WHO medicine names to upload to the STITCH 
# database to find the relevant chemical interaction network.

def read_medication_file():
    '''
    Reads the concomitant medication CSV file, and returns the list of WHO
    drug IDs.
    '''
    drug_set = set([])
    f = open('./data/parkinsons_loni/Concomitant_Medications.csv', 'r')
    f.readline()
    for line in reader(f):
        drug_id = line[23]
        if drug_id == '':
            continue
        # Add drug to the set.
        drug_set.add(drug_id)
    f.close()
    return list(drug_set)

def call_stitch_api_identifiers(i, drug_list):
    '''
    Call the STITCH API to get the identifiers corresponding to the WHO drug
    IDs in the LONI Parkinson's dataset.
    '''
    # Build the API link.
    api_link = 'http://stitch.embl.de/api/tsv/resolveList?identifiers='
    api_link += '%0D'.join(drug_list)
    api_link += '&species=9606'

    f = urllib.urlopen(api_link)
    myfile = f.read()
    out = open('%s/stitch_api_request_%d.txt' % (results_folder, i), 'w')
    out.write(myfile)
    out.close()

def map_who_to_stitch(who_to_stitch_dct, drug_list, i):
    '''
    Given a drug list and the returned API call number, map the WHO IDs to their
    corresponding STITCH identifiers.
    '''
    annotation_lst = []
    f = open('%s/stitch_api_request_%d.txt' % (results_folder, i), 'r')
    f.readline() # Skip the header line.
    for line in f:
        if 'small molecule' not in line:
            continue
        line = line.strip().split('\t')
        query_index, annotation = int(line[0]), line[5]
        # Query Index starts at -1 for the STITCH database.
        query_index += 1
        # Map the WHO drug ID to the STITCH identifier.
        who_drug_id = drug_list[query_index]
        assert who_drug_id not in who_to_stitch_dct
        who_to_stitch_dct[who_drug_id] = annotation
        annotation_lst += [annotation]
    f.close()
    return annotation_lst

def call_stitch_api_interactions(stitch_id_lst, i, protein_drug_set,
    drug_drug_set):
    '''
    Given a list of STITCH identifiers, get the interaction partners of the
    query items.
    '''
    api_link = 'http://stitch.embl.de/api/psi-mi-tab/interactionsList?identifiers='
    api_link += '%0D'.join(stitch_id_lst)
    api_link += '&species=9606&required_score=900&limit=100'

    f = urllib.urlopen(api_link)
    # out = open('%s/stitch_api_interactions_%d.txt' % (results_folder, i), 'w')
    # This if case is in situations in which there were no successful mappings.
    if stitch_id_lst != []:
        for line in f:
            line = line.strip().split('\t')
            node_a_type, node_b_type, node_a, node_b = line[:4]
            # Skip an interaction if neither node is in the query list.
            # We will fetch chemical-chemical interactions later.
            if node_a not in stitch_id_lst and node_b not in stitch_id_lst:
                continue
            # Drug-drug, protein-protein, or drug-protein.
            if 'ENSP' in node_a_type and 'ENSP' in node_b_type:
                # Skip protein-protein edges.
                continue
            elif 'ENSP' in node_a_type and 'CID' in node_b_type:
                protein_drug_set.add((node_a, node_b))
            elif 'CID' in node_a_type and 'ENSP' in node_b_type:
                protein_drug_set.add((node_b, node_a))
            else:
                assert 'CID' in node_a_type and 'CID' in node_b_type
                drug_drug_set.add((node_a, node_b))
    #         out.write('%s\t%s\n' % (node_a, node_b))
    # out.close()

def write_edges(fname, edge_set):
    out = open('%s/%s.txt' % (results_folder, fname), 'w')
    for node_a, node_b in edge_set:
        out.write('%s\t%s\n' % (node_a, node_b))
    out.close()

def main():
    # Generate directory.
    global results_folder
    results_folder = './data/api_calls'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    drug_list = read_medication_file()
    # write_stitch_upload_file(drug_set)

    who_to_stitch_dct = {}
    protein_drug_set, drug_drug_set = set([]), set([])
    # Process the drug list, 100 elements at a time.
    block_size = 10
    for i, sub_drug_lst in enumerate(zip(*[iter(drug_list)]*block_size)):
        call_stitch_api_identifiers(i, sub_drug_lst)
        
        # Map the WHO drug IDs to their corresponding STITCH names.
        annotation_lst = map_who_to_stitch(who_to_stitch_dct, sub_drug_lst, i)
        assert len(annotation_lst) <= block_size
        call_stitch_api_interactions(annotation_lst, i, protein_drug_set,
            drug_drug_set)
    with open('./data/who_to_stitch_dct.json', 'w') as fp:
        json.dump(who_to_stitch_dct, fp)
    # Write out interactions to file.
    write_edges('stitch_protein_drug', protein_drug_set)
    write_edges('stitch_drug_drug', drug_drug_set)

if __name__ == '__main__':
    main()