### Author: Edward Huang

from csv import reader
import os
import process_loni_parkinsons
from scipy.stats import fisher_exact

### This script checks whether a symptom is enriched with a SNP via Fisher's
### exact test.

data_folder = './data/parkinsons_loni'

def read_snp_fisher_file(fname):
    '''
    Given a Fisher's test file on finding significant SNPs, get the list of
    SNPs.
    '''
    snp_lst = []
    f = open(fname, 'r')
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        snp_lst += [line[0]]
    f.close
    return snp_lst

def read_mutation_file(fname, snp_patient_dct, snp_lst):
    '''
    Reads the mutation file, given a filename. Updates the given dictionary in
    place. snp_patient_dct contains the patients affected by the SNP key.
    '''
    f = open(fname, 'r')
    f.readline() # Skip header line.
    for line in f:
        patno, func, gene, exonic_func, snp = line.strip().split('\t')[1:]
        if snp not in snp_lst or 'nonsynonymous' not in exonic_func or func != 'exonic':
            continue
        # Update the patient dictionary.
        if snp not in snp_patient_dct:
            snp_patient_dct[snp] = set([])
        snp_patient_dct[snp].add(patno)
    f.close()

def read_clinical_diagnosis():
    '''
    Returns the major clinical diagnosis and additional notes for each patient.
    '''
    code_dct = process_loni_parkinsons.read_code_file()

    symptom_patient_dct = {}
    feat_name_lst = ['PATNO', 'EVENT_ID', 'PAG_NAME', 'PSLVL', 'PRIMDIAG']
    test_name_lst = ['DCRTREM', 'DCRIGID', 'DCBRADY', 'DFPGDIST']
    f = open('%s/Clinical_Diagnosis_and_Management.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    test_idx_lst = [header.index(test) for test in test_name_lst]
    for line in it:
        patno, event_id, pag_name, ps_lvl, prim_diag = (line[feat_idx] for
            feat_idx in feat_idx_lst)
        # Skip patients without UPDRS scores or non-baseline visits.
        if event_id != 'BL' or ps_lvl != '1':
            continue
        # Update the patient's binary column features.
        test_val_lst = [line[test_idx] for test_idx in test_idx_lst]
        for test_name_idx, test_val in enumerate(test_val_lst):
            if test_val == '1':
                test_name = test_name_lst[test_name_idx]
                if test_name not in symptom_patient_dct:
                    symptom_patient_dct[test_name] = set([])
                symptom_patient_dct[test_name].add(patno)

        # Skip "Other neurological disorder(s) (specify)".
        if prim_diag == '97':
            continue
        # Decode the primary diagnosis.
        prim_diag = code_dct[pag_name]['PRIMDIAG'][prim_diag]
        if prim_diag not in symptom_patient_dct:
            symptom_patient_dct[prim_diag] = set([])
        symptom_patient_dct[prim_diag].add(patno)
    f.close()
    return symptom_patient_dct

def read_medical_conditions():
    '''
    Maps PATNOs to current medical conditions.
    '''
    medical_condition_dct = {}
    f = open('%s/Current_Medical_Conditions_Log.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            term_idx = line.index('PT_NAME')
            resolved_idx = line.index('RESOLVD')
            continue
        patno, resolved = line[patno_idx], line[resolved_idx]
        term = line[term_idx]
        if term == '':
            continue
        if term not in medical_condition_dct:
            medical_condition_dct[term] = set([])
        medical_condition_dct[term].add(patno)
    f.close()
    return medical_condition_dct

def read_binary_tests():
    '''
    Maps PATNOs to their exam results, whether there are any abnormalities.
    '''
    binary_test_dct = {}

    exam_name_lst = ('DRMVIVID', 'DRMAGRAC', 'DRMNOCTB', 'SLPLMBMV',
        'SLPINJUR', 'DRMVERBL', 'DRMFIGHT', 'DRMUMV', 'DRMOBJFL',
        'MVAWAKEN', 'DRMREMEM', 'SLPDSTRB', 'STROKE', 'HETRA', 'PARKISM',
        'RLS', 'NARCLPSY', 'DEPRS', 'EPILEPSY', 'BRNINFM')
    fname = 'REM_Sleep_Disorder_Questionnaire'

    f = open('%s/%s.csv' % (data_folder, fname), 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    patno_idx, event_idx = header.index('PATNO'), header.index('EVENT_ID')
    exam_idx_lst = [header.index(exam) for exam in exam_name_lst]
    for line in it:
        patno, event_id = line[patno_idx], line[event_idx]
        if event_id != 'BL':
            continue

        # Again, a duplicate entry.
        # Get abnormal values ('1') for each test.
        exam_val_lst = [line[exam_idx] for exam_idx in exam_idx_lst]
        for exam_name_idx, exam_val in enumerate(exam_val_lst):
            if exam_val == '1':
                exam_name = exam_name_lst[exam_name_idx]
                if exam_name not in binary_test_dct:
                    binary_test_dct[exam_name] = set([])
                # binary_test_dct[patno] += [(exam_name, 1)]
                binary_test_dct[exam_name].add(patno)
    f.close()
    return binary_test_dct

def compute_snp_symptom_fisher(snp_patient_dct, symptom_patient_dct):
    '''
    Between every SNP-symptom pair, compute a Fisher's test.
    '''
    # Get the patients from each dataset to build the patient universe.
    status_dct = process_loni_parkinsons.read_patient_status()
    status_patient_set = set(status_dct.keys())

    status_patient_set = status_patient_set.union(process_loni_parkinsons.get_updrs_dct()[0].keys())
    status_patient_set = status_patient_set.union(set.union(*symptom_patient_dct.values()))
    status_patient_set = status_patient_set.union(set.union(*snp_patient_dct.values()))

    # Compute the fisher dictionary.
    fisher_tups = []
    for snp in snp_patient_dct:
        snp_patients = snp_patient_dct[snp]
        for symptom in symptom_patient_dct:
            symptom_patients = symptom_patient_dct[symptom]
            # Get the Fisher table.
            snp_and_symptom = len(snp_patients.intersection(symptom_patients))
            snp_not_symptom = len(snp_patients.difference(symptom_patients))
            symptom_not_snp = len(symptom_patients.difference(snp_patients))
            neither = len(status_patient_set) - len(snp_patients.union(symptom_patients))
            f_table = [[snp_and_symptom, snp_not_symptom], [symptom_not_snp, neither]]
            o_r, p_value = fisher_exact(f_table, alternative='greater')
            if p_value < 0.01:
                # print snp, symptom, p_value, f_table
                fisher_tups += [(snp, symptom, f_table, p_value)]
    return sorted(fisher_tups, key=lambda x: x[-1])

def main():
    # First, get the SNPs deemed to be significantly enriched in PD patients.
    snp_lst = []
    for fname in ('snp_fisher_test_wes_ppmi', 'snp_fisher_test_wes_hard_ignore'):
        snp_lst += read_snp_fisher_file('./data/ppmi/snp_files/%s.tsv' % fname)
        # break # TODO
    # Next, for each SNP, get the patients that have the SNP.
    snp_patient_dct = {}
    wes_folder = './data/annovar_annotate_output_wes_patient_info'
    for fname in os.listdir(wes_folder):
        read_mutation_file('%s/%s' % (wes_folder, fname), snp_patient_dct, snp_lst)
        # break # TODO

    # Next, for each symptom, get the patients that have the symptom.
    symptom_patient_dct = read_clinical_diagnosis()
    symptom_patient_dct.update(read_medical_conditions())
    symptom_patient_dct.update(read_binary_tests())

    # Compute Fisher's test between symptom and SNP dictionaries.
    # TODO: get status dct. combine with SNP patients. union is patient universe.
    fisher_tups = compute_snp_symptom_fisher(snp_patient_dct, symptom_patient_dct)    

    # Write tuples out to file.
    out = open('./data/ppmi/snp_files/snp_symptom_enrichments.txt', 'w')
    for snp, symptom, f_table, p_value in fisher_tups:
        out.write('%s\t%s\t%d\t%d\t%d\t%d\t%f\n' % (snp, symptom, f_table[0][0],
            f_table[0][1], f_table[1][0], f_table[1][1], p_value))
    out.close()

if __name__ == '__main__':
    main()