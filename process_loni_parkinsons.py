### Author: Edward Huang

from csv import reader

### This script calls on LONI USC Parkinson's disease data.

data_folder = './data/parkinsons_loni'

def read_updrs_file(fname):
    '''
    Reads the scores for all attributes in the MDS-UPDRS file.
    '''
    f = open('%s/%s.csv' % (data_folder, fname), 'r')
    for i, line in enumerate(reader(f)):
        if i == 0:
            patno_idx = line.index('PATNO')
            infodt_idx = line.index('INFODT')
            orig_entry_idx = line.index('ORIG_ENTRY')
            continue
        # Scores are between INFODT and ORIG_ENTRY.
        # TODO: This is different for part III.
        line = line[infodt_idx + 1:orig_entry_idx]
        while '' in line:
            line.remove('')
        score_list = map(int, line)
        print score_list
    f.close()

def read_code_file():
    '''
    Maps the code of each file to the decoded name.
    Each key maps to another dictionary, depending on PAG_NAME.
    '''
    code_dct = {}
    f = open('%s/Code_List.csv' % data_folder, 'r')
    f.readline()
    for line in reader(f):
        pag_name = line[0]
        # Make a nested dictionary for the page name.
        if pag_name not in code_dct:
            code_dct[pag_name] = {}
        itm_name = line[1]
        # Make another nested dictionary for the code name.
        if itm_name not in code_dct[pag_name]:
            code_dct[pag_name][itm_name] = {}
        # Decode the code.
        code_dct[pag_name][itm_name][line[3]] = line[4]
    f.close()
    # Run some tests.
    assert code_dct['AE']['AESEVER']['3'] == 'severe'
    assert code_dct['FBIMAG']['FBATVS']['U04'] == 'Unexpected 4'
    assert code_dct['ind_upsit']['SCENT_37_CORRECT']['0'] == 'Incorrect'
    return code_dct

def read_adverse_events():
    '''
    Returns a dictionary mapping patient IDs to their adverse events.
    Key: PATNO -> str
    Value: list of (adverse event, severity) tuples -> list((str, float))
    '''
    adverse_event_dct = {}
    f = open('%s/Adverse_Event_Log.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            term_idx = line.index('AETERM')
            severity_idx = line.index('AESEVER')
            continue
        patno = line[patno_idx]
        term = line[term_idx]
        severity = float(line[severity_idx])
        # Add patient to the dictionary.
        if patno not in adverse_event_dct:
            adverse_event_dct[patno] = []
        adverse_event_dct[patno] += [(term, severity)]
    f.close()
    # Some tests.
    assert adverse_event_dct['3405'] == [('POST LP HEADACHE', 1.0), (
        'POST LP ABDOMINAL CRAMPS', 1.0)]
    assert adverse_event_dct['52215'] == [('ASYMPTOMATIC HYPERTENSION', 1
        ), ('SORENESS AT LP SITE', 1)]
    return adverse_event_dct

# def read_biospecimen_analysis():
def read_test_analysis(test_type):
    '''
    Returns a dictionary mapping PATNOs to test results.
    '''
    col_dct = {'biospecimen':('TESTNAME', 'TESTVALUE', 'UNITS'), 'hematology':(
        'LTSTNAME', 'LSIRES', 'LSIUNIT')}
    fname_dct = {'biospecimen':'Biospecimen_Analysis_Results', 'hematology':
        'Blood_Chemistry___Hematology'}

    test_result_dct = {}
    # Make sure that the units are consistent across all tests.
    test_unit_dct = {}
    f = open('%s/%s.csv' % (data_folder, fname_dct[test_type]), 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            name_idx, value_idx, unit_idx = [line.index(col
                ) for col in col_dct[test_type]]
            continue
        patno = line[patno_idx]
        test_name = line[name_idx]
        try:
            test_value = float(line[value_idx])
        except Exception:
            continue
        units = line[unit_idx]
        # Update the units for current test.
        if test_name not in test_unit_dct:
            test_unit_dct[test_name] = units
        # Skip tests that are not consistent in units.
        if units != test_unit_dct[test_name]:
            continue
        # Update the test results dictionary.
        if patno not in test_result_dct:
            test_result_dct[patno] = []
        test_result_dct[patno] += [(test_name, test_value)]
    f.close()
    return test_result_dct

def read_clinical_diagnosis(code_dct):
    '''
    Returns the major clinical diagnosis and additional notes for each patient.
    '''
    clinical_diagnosis_dct = {}
    # Relevant column names.
    col_list = ['DCNOMTR', 'DCRTREM', 'DCRIGID', 'DCBRADY', 'DFPGDIST']

    f = open('%s/Clinical_Diagnosis_and_Management.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            pag_name_idx = line.index('PAG_NAME')
            binary_columns = [line.index(col) for col in col_list]
            prim_diag_idx = line.index('PRIMDIAG')
            continue
        patno = line[patno_idx]
        # Add the patient to the dictionary.
        if patno not in clinical_diagnosis_dct:
            clinical_diagnosis_dct[patno] = set([])
        # Update the patient's binary column features.
        for col_name_idx, col_idx in enumerate(binary_columns):
            if line[col_idx] == '1':
                clinical_diagnosis_dct[patno].add((col_list[col_name_idx], 1))
        # Update the patient's primary diagnosis.
        prim_diag = line[prim_diag_idx]
        # Skip "Other neurological disorder(s) (specify)".
        if prim_diag == '97':
            continue
        # Decode the primary diagnosis.
        pag_name = line[pag_name_idx]
        prim_diag = code_dct[pag_name]['PRIMDIAG'][prim_diag]
        clinical_diagnosis_dct[patno].add((prim_diag, 1))
    f.close()
    return clinical_diagnosis_dct

def run_tests(biospecimen_dct, hematology_dct, clinical_diagnosis_dct):
    # Test the biospecimen dictionary.
    assert ('FBXO7-001', 159) in biospecimen_dct['3000']
    assert ('Serum IGF-1', 148.2) in biospecimen_dct['3386']
    assert ('HSPA8 (rep 1)', 21.50342178) in biospecimen_dct['3821']
    assert ('UBE2K (rep 2)', 24.21042061) in biospecimen_dct['4058']

    # Test hematology dictionary.
    assert ('RBC', 4.3) in hematology_dct['3404']
    assert ('Serum Bicarbonate', 26.7) in hematology_dct['3601']
    assert ('Creatinine (Rate Blanked)', 80) in hematology_dct['4011']

    # Test the clinical diagnosis and managemenet dictionary.
    assert ('Motor neuron disease with parkinsonism', 1
        ) in clinical_diagnosis_dct['3425']
    assert ('No PD nor other neurological disorder', 1
        ) in clinical_diagnosis_dct['4053']
    assert ('DCNOMTR', 1) in clinical_diagnosis_dct['3956']
    print 'ye' # TODO

def main():
    # read_updrs_file('MDS_UPDRS_Part_I')
    # read_updrs_file('MDS_UPDRS_Part_I__Patient_Questionnaire')
    # read_updrs_file('MDS_UPDRS_Part_II__Patient_Questionnaire')
    # # part III has a slightly different format.
    # read_updrs_file('MDS_UPDRS_Part_III__Post_Dose_')
    # read_updrs_file('MDS_UPDRS_Part_IV')

    adverse_event_dct = read_adverse_events()
    # biospecimen_dct = read_biospecimen_analysis()
    biospecimen_dct = read_test_analysis('biospecimen')
    hematology_dct = read_test_analysis('hematology')

    code_dct = read_code_file()
    clinical_diagnosis_dct = read_clinical_diagnosis(code_dct)

    run_tests(biospecimen_dct, hematology_dct, clinical_diagnosis_dct)

if __name__ == '__main__':
    main()