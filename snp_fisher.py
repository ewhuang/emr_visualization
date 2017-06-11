### Author: Edward Huang

import argparse
import math
import operator
from process_loni_parkinsons import read_patient_status
from scipy.stats import fisher_exact

### This script determines which of the SNPs are highly enriched in patients
### with Parkinson's disease with Fisher's exact test.

def read_SNP_background_freq(snp_bkg_dct, fname):
    '''
    Reads the background frequencies in all patients for SNPs.
    '''
    assert fname in ['SNP', 'indel', 'wes']
    f = open('./data/ppmi/filter_output_merge_%s.tsv' % fname, 'r')
    f.readline() # Skip the header line.
    for line in f:
        line = line.strip().split('\t')
        assert len(line) == 12
        freq, func, gene, snp = line[6:10]
        cadd = line[11]
        # snp, freq, func, gene, cadd = line[9], float(line[6]), line[7], line[8], line[11]
        if snp == '.': # Skip empty SNPs.
            continue
        snp_bkg_dct[snp] = (float(freq), func, gene, cadd) # TODO: some SNPs appear more than once.
    f.close()

def read_ppmi_mutations(patient_snp_dct, snp_ppmi_dct, fname):
    '''
    Reads the PPMI mutation file, and records the genes with mutations for each
    patient. No header line.
    '''
    assert fname in ['indels', 'mutation']
    f = open('./data/ppmi/PPMI_%s.txt' % fname, 'r') # new data from peng-gpu.
    for line in f:
        patno, mutated_gene, att_1, exonic, snp = line.strip().split('\t')
        if snp == '.':
            continue
        # Add the mutation tag line to the gene.
        if patno not in patient_snp_dct:
            patient_snp_dct[patno] = set([])
        patient_snp_dct[patno].add(snp)
        if ';' in exonic: # Weird lines that have stuff like exonic;exonic
            exonic = exonic.split(';')[0]
        # Check that repetitions are already in the dictionary.
        if snp in snp_ppmi_dct:
            assert snp_ppmi_dct[snp] == (exonic, mutated_gene)
        snp_ppmi_dct[snp] = (exonic, mutated_gene)
    f.close()

def write_fisher_sns(snp_bkg_dct, snp_type):
    '''
    Writes out the SNPs according to mutation patterns of PD/Healthy Control
    patients in the PD dataset. Makes a table for PD/HC patients with/without
    the SNP, and then computes Fisher's exact test to determine if the SNP
    is enriched in PD patients.
    '''
    status_dct = read_patient_status()

    patient_snp_dct, snp_ppmi_dct = {}, {}
    read_ppmi_mutations(patient_snp_dct, snp_ppmi_dct, 'indels')
    read_ppmi_mutations(patient_snp_dct, snp_ppmi_dct, 'mutation')

    out = open('./data/ppmi/snp_fisher_test_%s.tsv' % snp_type, 'w')
    out.write('SNP\tPD_SNP\tPD_no_SNP\tHC_SNP\tHC_NO_SNP\tp-value\tfreq\tfunc\tgene\tCADD_phred\n')
    for snp in snp_ppmi_dct:
        # Get the Fisher's table for each SNP.
        if snp not in snp_bkg_dct:
            continue
        freq = snp_bkg_dct[snp][0]
        if freq > 0.05:
            continue
        hc_snp = int(math.ceil(65000 * freq))
        hc_no_snp = 65000 - hc_snp
        # Get the number of PD/HC patients with/without the SNP.
        pd_snp, pd_no_snp = 0.0, 0.0
        for patno in status_dct:
            if status_dct[patno] != 'PD':
                continue
            # A patient does not have the SNP if the patient has no SNPs.
            # PD and without SNP.
            if patno not in patient_snp_dct or snp not in patient_snp_dct[patno]:
                pd_no_snp += 1
            elif snp in patient_snp_dct[patno]:
                pd_snp += 1
            # HC and without SNP.
            # elif patient_status == 'Control' and no_snp:
            #     hc_no_snp += 1
            # # PD and with SNP.
            # elif patient_status == 'PD':
            #     assert snp in patient_snp_dct[patno]
            #     pd_snp += 1
            # elif patient_status == 'Control':
            #     assert snp in patient_snp_dct[patno]
            #     hc_snp += 1
            # else:
            #     # If it's not in the table, then that means the patient is SWEDD.
            #     assert patient_status == 'SWEDD'
        o_r, p_val = fisher_exact([[pd_snp, hc_snp], [pd_no_snp, hc_no_snp]],
            alternative='greater')
        if p_val < 0.1:
            out.write('%s\t%d\t%d\t%d\t%d\t%g' % (snp, pd_snp, pd_no_snp,
                hc_snp, hc_no_snp, p_val))
            if snp in snp_bkg_dct:
                out.write('\t%s\t%s\t%s\t%s' % snp_bkg_dct[snp])
                if snp not in ['rs751730130', 'rs11272867']:
                    assert snp_bkg_dct[snp][1:3] == snp_ppmi_dct[snp]
            else:
                out.write('\t\t%s\t%s\t' % snp_ppmi_dct[snp])
            out.write('\n')
    out.close()

def write_naive_snps(snp_bkg_dct, snp_type):
    '''
    Naively writes out SNPs that occur rarely in general population, deeming
    them as "interesting".
    '''
    # Sort SNPs by their frequency.
    sorted_dct = {}
    for snp in snp_bkg_dct:
        freq = snp_bkg_dct[snp][0]
        if freq < 0.05:
            sorted_dct[snp] = freq
    sorted_dct = sorted(sorted_dct.items(), key=operator.itemgetter(1))
    # Write out the sorted SNPs, along with their full information.
    out = open('./data/ppmi/snp_healthy_freq_%s.tsv' % snp_type, 'w')
    out.write('SNP\tfreq\tfunc\tgene\tCADD_phred\n')
    for snp, freq in sorted_dct:
        out.write('%s\t%s\n' % (snp, '\t'.join(map(str, snp_bkg_dct[snp]))))
    out.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--snp_type', help='type of SNP. [\'wgs\', \'wes\']')
    return parser.parse_args()

def main():
    args = parse_args()
    snp_type = args.snp_type

    snp_bkg_dct = {}
    if snp_type == 'wgs':
        read_SNP_background_freq(snp_bkg_dct, 'indel')
        read_SNP_background_freq(snp_bkg_dct, 'SNP')
    elif snp_type == 'wes':
        read_SNP_background_freq(snp_bkg_dct, 'wes')

    write_fisher_sns(snp_bkg_dct, snp_type)
    write_naive_snps(snp_bkg_dct, snp_type)

if __name__ == '__main__':
    main()