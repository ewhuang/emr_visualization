### Author: Edward Huang

import argparse
import json
import math
from multiprocessing import Pool
import os
from process_loni_parkinsons import read_patient_status
from scipy.stats import fisher_exact

### This script determines which of the SNPs are highly enriched in patients
### with Parkinson's disease with Fisher's exact test.

status_dct = read_patient_status()
num_pd_patients = status_dct.values().count('PD')
num_hc_patients = status_dct.values().count('Control')

def generate_directories():
    global results_folder
    results_folder = './data/ppmi/snp_files'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

def generate_file_suffix():
    global suffix
    suffix = '%s_%s' % (args.snp_type, args.healthy_control)
    if args.healthy_control == 'hard':
        suffix += '_%s' % args.bkg_strat

def read_SNP_background_freq(fname):
    '''
    Reads the background frequencies in all patients for SNPs.
    Key: SNP -> str
    Value: (frequency, function, gene, confidence) -> (float, str, str, str)
    '''
    snp_bkg_dct = {}
    assert fname in ['SNP', 'indel', 'wes']
    f = open('./data/ppmi/filter_output_merge_%s.tsv' % fname, 'r')
    f.readline() # Skip the header line.
    for line in f:
        line = line.strip().split('\t')
        assert len(line) == 12
        freq, func, gene, snp, cadd_raw, cadd_phred = line[6:]

        if snp == '.': # Skip empty SNPs.
            continue

        # For repeated SNPs, take the max frequency.
        freq = float(freq)
        if snp not in snp_bkg_dct or freq > snp_bkg_dct[snp][0]:
            snp_bkg_dct[snp] = (freq, func, gene, cadd_phred)
    f.close()
    return snp_bkg_dct

def read_mutation_file(snp_pd_count_dct, snp_ppmi_dct, fname):
    '''
    Reads the mutation file, given a filename. Updates the given dictionaries
    in place. snp_pd_count_dct contains the number of PD patients given a SNP
    key. snp_ppmi_dct contains tuples of (exonic, gene) for each SNP key.
    '''
    f = open(fname, 'r')
    if 'wes' in fname:
        f.readline() # Skip header line for WES.
    for line in f:
        if 'wes' in fname:
            patno, func, gene, exonic_func, snp = line.strip().split('\t')[1:]
        else:
            patno, gene, exonic_func, func, snp = line.strip().split('\t')
        # Skip blank SNPs. Only use exonic, nonsynonymous SNPs.
        if snp == '.' or func != 'exonic' or 'nonsynonymous' not in exonic_func:
            continue
        snp_ppmi_dct[snp] = (func, exonic_func, gene)

        # Only keep track of PD/HC patients with the SNP.
        if patno not in status_dct:
            continue
        patient_status = status_dct[patno]
        if patient_status not in ['PD', 'Control']:
            continue
        # Update the count dictionary.
        if snp not in snp_pd_count_dct:
            # First element is the set of PD patients. Second element is HC.
            snp_pd_count_dct[snp] = [set([]), set([])]
        if patient_status == 'PD':
            snp_pd_count_dct[snp][0].add(patno)
        elif patient_status == 'Control':
            snp_pd_count_dct[snp][1].add(patno)
    f.close()

def write_snp_and_count_dct():
    '''
    Writes out the number of PD and HC patients for each SNP. Also writes out
    the function and exonic information for each SNP. Both dictionaries are
    dumpd to JSON files.
    '''
    snp_pd_count_dct, snp_ppmi_dct = {}, {}
    if args.snp_type == 'wgs':
        read_mutation_file(snp_pd_count_dct, snp_ppmi_dct, './data/ppmi/PPMI_indels.txt')
        read_mutation_file(snp_pd_count_dct, snp_ppmi_dct, './data/ppmi/PPMI_mutation.txt')
    elif args.snp_type == 'wes':
        wes_folder = './data/annovar_annotate_output_wes_patient_info'
        for fname in os.listdir(wes_folder):
            read_mutation_file(snp_pd_count_dct, snp_ppmi_dct,
                '%s/%s' % (wes_folder, fname))

    # Convert sets to lengths for count dictionary.
    for snp in snp_pd_count_dct:
        snp_pd_count_dct[snp] = map(len, snp_pd_count_dct[snp])
    # Dump the resulting dictionaries to file.
    with open('%s/snp_pd_count_dct_%s.json' % (results_folder, args.snp_type), 'w') as fp:
        json.dump(snp_pd_count_dct, fp)
    with open('%s/snp_ppmi_dct_%s.json' % (results_folder, args.snp_type), 'w') as fp:
        json.dump(snp_ppmi_dct, fp)

def compute_fisher_value(snp):
    # This first if statement computes the HC patients.
    # For hard-coded healthy control mode, use the background frequency.
    if args.healthy_control == 'hard':
        if snp in snp_bkg_dct:
            freq = snp_bkg_dct[snp][0]
        elif args.bkg_strat == 'use':
            # In this mode, assume missing frequencies are 0.
            freq = 0
        elif args.bkg_strat == 'ignore':
            # In this mode, ignore SNPs with missing frequencies.
            return
        # Skip SNPs with high background frequency.
        if freq > 0.05:
            return
        # Hardcoding 65000 as the number of patients without PD.
        hc_snp = int(math.ceil(65000 * freq))
        # Number of healthy control patients without the SNP.
        hc_no_snp = 65000 - hc_snp
    elif args.healthy_control == 'ppmi':
        # Only use PPMI for SNPs that don't have background frequencies.
        if snp in snp_bkg_dct:
            return
        hc_snp = 0
        if snp in snp_pd_count_dct:
            hc_snp = snp_pd_count_dct[snp][1]
        hc_no_snp = num_hc_patients - hc_snp
        assert hc_no_snp >= 0

    # For all modes, getting PD patients with/without the SNP.
    pd_snp = 0
    if snp in snp_pd_count_dct:
        pd_snp = snp_pd_count_dct[snp][0]
    pd_no_snp = num_pd_patients - pd_snp
    assert pd_no_snp >= 0

    o_r, p_val = fisher_exact([[pd_snp, pd_no_snp], [hc_snp, hc_no_snp]],
        alternative='greater')
    # Getting the filename.
    out = open('%s/snp_fisher_test_%s_unsorted.tsv' % (results_folder, suffix), 'a')
    out.write('%s\t%d\t%d\t%d\t%d\t%g\n' % (snp, pd_snp, pd_no_snp, hc_snp,
        hc_no_snp, p_val))
    out.close()

def write_raw_fisher_snps():
    '''
    Writes out the SNPs according to mutation patterns of PD/Healthy Control
    patients in the PD dataset. Makes a table for PD/HC patients with/without
    the SNP, and then computes Fisher's exact test to determine if the SNP
    is enriched in PD patients.
    '''
    global snp_pd_count_dct
    with open('%s/snp_pd_count_dct_%s.json' % (results_folder, args.snp_type),
        'r') as fp:
        snp_pd_count_dct = json.load(fp)
    fp.close()
    with open('%s/snp_ppmi_dct_%s.json' % (results_folder, args.snp_type),
        'r') as fp:
        snp_ppmi_dct = json.load(fp)
    fp.close()

    print 'Writing out to file.'
    pool = Pool(processes=20, maxtasksperchild=1000)
    with open('%s/snp_fisher_test_%s_unsorted.tsv' % (results_folder, suffix), 'w') as out:
        pool.map(compute_fisher_value, snp_pd_count_dct.keys())
    return snp_ppmi_dct

def sort_raw_fisher_file(snp_ppmi_dct):
    '''
    Get the raw, unsorted file written by the pool task. Sort it, then write
    it back out with the SNPs information.
    '''
    # Read the unsorted file.
    snp_tuple_lst = []
    f = open('%s/snp_fisher_test_%s_unsorted.tsv' % (results_folder, suffix), 'r')
    for line in f:
        line = line.split()
        assert len(line) == 6
        # Convert the p-value to a float.
        line[-1] = float(line[-1])
        snp_tuple_lst += [line]
    f.close()
    # Sort the SNP tuples by p-value, which is the last element in the line.
    snp_tuple_lst = sorted(snp_tuple_lst, key=lambda x:x[-1])
    num_snps = len(snp_tuple_lst) # For Bonferroni.

    # Set p-value threshold, depending on the mode.
    if args.healthy_control == 'ppmi':
        p_val_thresh = 0.01
    else:
        p_val_thresh = 0.05

    out = open('%s/snp_fisher_test_%s.tsv' % (results_folder, suffix), 'w')
    out.write('SNP\tPD_SNP\tPD_no_SNP\tHC_SNP\tHC_NO_SNP\tp-value\tfreq\tfunc\texonic_func\tgene\tCADD_phred\n')
    for tup in snp_tuple_lst:
        # p_val = snp_fisher_dct[key]
        p_val = tup[-1]
        # Bonferroni test correction for hard-coded healthy patients.
        if args.healthy_control == 'hard':
            p_val *= num_snps
        if p_val < p_val_thresh:
            # Write the Fisher's table.
            tup[-1] = p_val
            tup = map(str, tup)
            out.write('%s' % '\t'.join(tup))
            # Write out the extra information regarding the SNP.
            snp = tup[0]
            if snp in snp_bkg_dct:
                freq, func, gene, cadd_phred = snp_bkg_dct[snp]
                out.write('\t%s\t%s\t%s\t%s\t%s' % (freq, func,
                    snp_ppmi_dct[snp][1], gene, cadd_phred))
            else:
                func, exonic_func, gene = snp_ppmi_dct[snp]
                out.write('\t\t%s\t%s\t%s\t' % (func, exonic_func, gene))
            out.write('\n')
    out.close()

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--snp_type', choices=['wgs', 'wes'], required=True,
        help='source of SNP data.')
    parser.add_argument('-c', '--healthy_control', choices=['hard', 'ppmi'],
        required=True, help='Hard-coded 65,000 or Control patients.')
    parser.add_argument('-b', '--bkg_strat', choices=['ignore', 'use'],
        help='use SNPs missing from background frequency files.')
    args = parser.parse_args()
    if args.healthy_control == 'ppmi':
        assert args.bkg_strat == None

def main():
    # Set the global variables.
    generate_directories()
    parse_args()
    generate_file_suffix()

    # Clear the unsorted file, since we need to append to it later.
    unsorted_fname = '%s/snp_fisher_test_%s_unsorted.tsv' % (results_folder, suffix)
    if os.path.exists(unsorted_fname):
        os.remove(unsorted_fname)

    # Set up the background frequency dictionary for SNPs.
    global snp_bkg_dct
    if args.snp_type == 'wgs':
        snp_bkg_dct = read_SNP_background_freq('indel')
        snp_bkg_dct.update(read_SNP_background_freq('SNP'))
    elif args.snp_type == 'wes':
        snp_bkg_dct = read_SNP_background_freq('wes')

    # Compute the SNP enrichments.
    if not os.path.exists('%s/snp_pd_count_dct_%s.json' % (results_folder, args.snp_type)):
        write_snp_and_count_dct()

    snp_ppmi_dct = write_raw_fisher_snps()
    sort_raw_fisher_file(snp_ppmi_dct)

if __name__ == '__main__':
    main()
