### Author: Edward Huang

from csv import reader

### This script calls on LONI USC Parkinson's disease data.

data_folder = './data/parkinsons_loni'

# TODO.
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
    '''
    adverse_event_dct = {}
    f = open('%s/Adverse_Event_Log.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            term_idx = line.index('AETERM')
            severity_idx = line.index('AESEVER')
            related_idx = line.index('AERELAT')
            continue
        isRelated = line[related_idx]
        # TODO: Skip adverse events that are unrelated to the study.
        if isRelated not in ['4', '5']:
            continue
        patno = line[patno_idx]
        term = line[term_idx]
        severity = float(line[severity_idx])
        # Create nested list mapping adverse events to their severities.
        if patno not in adverse_event_dct:
            adverse_event_dct[patno] = {}
        # Create list of the adverse event's severities.
        if term not in adverse_event_dct[patno]:
            adverse_event_dct[patno][term] = []
        adverse_event_dct[patno][term] += [severity]
    f.close()
    return adverse_event_dct

def read_test_analysis(test_type):
    '''
    Returns a dictionary mapping PATNOs to test results.
    '''
    col_dct = {'biospecimen':('TESTNAME', 'TESTVALUE', 'UNITS'), 'hematology':(
        'LTSTNAME', 'LSIRES', 'LSIUNIT'), 'medication':('CMTRT', 'CMDOSE',
        'CMDOSU')}
    fname_dct = {'biospecimen':'Biospecimen_Analysis_Results', 'hematology':
        'Blood_Chemistry___Hematology', 'medication':'Concomitant_Medications'}

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
            test_result_dct[patno] = {}
        if test_name not in test_result_dct[patno]:
            test_result_dct[patno][test_name] = []
        test_result_dct[patno][test_name] += [test_value]
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
            clinical_diagnosis_dct[patno] = {}
        # Update the patient's binary column features.
        for col_name_idx, col_idx in enumerate(binary_columns):
            if line[col_idx] == '1':
                clinical_diagnosis_dct[patno][col_list[col_name_idx]] = 1
        # Update the patient's primary diagnosis.
        prim_diag = line[prim_diag_idx]
        # Skip "Other neurological disorder(s) (specify)".
        if prim_diag == '97':
            continue
        # Decode the primary diagnosis.
        pag_name = line[pag_name_idx]
        prim_diag = code_dct[pag_name]['PRIMDIAG'][prim_diag]
        clinical_diagnosis_dct[patno][prim_diag] = 1
    f.close()
    return clinical_diagnosis_dct

def read_cognitive_assessments():
    '''
    Returns a dictionary mapping PATNOs to cognitive assessment times.
    '''
    col_list = ['HVLTRTTM', 'HVLTDRTM', 'LNORNTTM', 'SFTTM', 'LNSTM', 'SDMTM']

    cognitive_assessment_dct = {}
    f = open('%s/Cognitive_Assessments.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            time_columns = [line.index(col) for col in col_list]
            continue
        patno = line[patno_idx]
        # Add patient to dictionary.
        if patno not in cognitive_assessment_dct:
            cognitive_assessment_dct[patno] = {}
        # Update patient's cognitive assessments.
        for col_name_idx, col_idx in enumerate(time_columns):
            if line[col_idx] == '':
                continue
            test_time = line[col_idx]
            test_time = sum(x * int(t) for x, t in zip([60, 1, 0.01],
                test_time.split(":")))
            test_name = col_list[col_name_idx]
            if test_name not in cognitive_assessment_dct[patno]:
                cognitive_assessment_dct[patno][test_name] = []
            cognitive_assessment_dct[patno][test_name] += [(test_time)]
    f.close()
    return cognitive_assessment_dct

def read_cognitive_categorizations():
    '''
    Returns a dictionary mapping PATNOs to cognitive categorizations.
    '''
    cognitive_categorization_dct = {}
    f = open('%s/Cognitive_Categorization.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            COGDECLN_idx = line.index('COGDECLN')
            FNCDTCOG_idx = line.index('FNCDTCOG')
            COGSTATE_idx = line.index('COGSTATE')
            COGDXCL_idx = line.index('COGDXCL')
            continue
        patno = line[patno_idx]
        confidence = line[COGDXCL_idx] # Level of confidence of diagnosis.
        # Skip diagnoses with low confidence.
        if confidence != '1':
            continue
        # Update categorization dictionary.
        if patno not in cognitive_categorization_dct:
            cognitive_categorization_dct[patno] = {}
        if line[COGDECLN_idx] == '1':
            cognitive_categorization_dct[patno]['COGDECLN'] = 1
        if line[FNCDTCOG_idx] == '1':
            cognitive_categorization_dct[patno]['FNCDTCOG'] = 1
        cog_state = line[COGSTATE_idx]
        if 'COGSTATE' not in cognitive_categorization_dct[patno]:
            cognitive_categorization_dct[patno]['COGSTATE'] = []
        cognitive_categorization_dct[patno]['COGSTATE'] += [int(cog_state)]
    f.close()
    return cognitive_categorization_dct

def main():
    # read_updrs_file('MDS_UPDRS_Part_I')
    # read_updrs_file('MDS_UPDRS_Part_I__Patient_Questionnaire')
    # read_updrs_file('MDS_UPDRS_Part_II__Patient_Questionnaire')
    # # part III has a slightly different format.
    # read_updrs_file('MDS_UPDRS_Part_III__Post_Dose_')
    # read_updrs_file('MDS_UPDRS_Part_IV')

    adverse_event_dct = read_adverse_events()
    biospecimen_dct = read_test_analysis('biospecimen')
    hematology_dct = read_test_analysis('hematology')

    code_dct = read_code_file()
    clinical_diagnosis_dct = read_clinical_diagnosis(code_dct)

    cognitive_assessment_dct = read_cognitive_assessments()

    cognitive_categorization_dct = read_cognitive_categorizations()

    medication_dct = read_test_analysis('medication')

    ### TESTS___________________________________________________________________

    # Test the adverse event dictionary.
    assert adverse_event_dct['3226'] == {'BACK SORENESS':[1.0, 2.0]}
    assert adverse_event_dct['52215'] == {'SORENESS AT LP SITE':[1.0]}
    assert adverse_event_dct['60036'] == {'HEADACHE':[2.0, 1.0, 1.0, 3.0, 1.0],
        'NAUSEA':[1.0], 'SWEATING':[1.0], 'SHIVERS':[1.0]}

    # Test the biospecimen dictionary.
    assert biospecimen_dct['3003']['Abeta 42'] == [271.3, 277.5, 310.8]
    assert biospecimen_dct['3410']['HDL'] == [55, 46, 32]
    assert biospecimen_dct['3523']['ZNF746'] == [143, 268]

    # Test hematology dictionary.
    assert hematology_dct['3252']['Calcium (EDTA)'] == [2.57, 2.42, 2.45, 2.42,
        2.5, 2.57, 2.5]
    assert hematology_dct['3276']['Serum Sodium'] == [139, 139, 141, 138, 139]
    assert hematology_dct['10874']['Total Protein'] == [63, 73, 67]

    # Test the clinical diagnosis and managemenet dictionary.
    assert clinical_diagnosis_dct['3425']['Motor neuron disease with parkinsonism'] == 1
    assert clinical_diagnosis_dct['4053']['No PD nor other neurological disorder'] == 1
    assert clinical_diagnosis_dct['3956']['DCNOMTR'] == 1

    # Test cognitive assessments dictionary.
    assert cognitive_assessment_dct['3154']['HVLTRTTM'] == [577, 919, 593]
    assert cognitive_assessment_dct['3951']['HVLTDRTM'] == [822, 557]
    assert cognitive_assessment_dct['41749']['LNSTM'] == [685]
    
    # Test cognitive categorizations dictionary.
    assert cognitive_categorization_dct['3000']['COGSTATE'] == [1]
    assert 'FNCDTCOG' not in cognitive_categorization_dct['3057']
    assert cognitive_categorization_dct['60073']['COGDECLN'] == 1

    # Test concomitant medications dictionary.
    assert '"ESTROGEN GEL"' not in medication_dct['50157']
    assert medication_dct['51551']['ACYCLOVIR'] == [500, 400]
    assert medication_dct['52373']['ABILIFY'] == [2, 5]

    print 'Finished tests!'

if __name__ == '__main__':
    main()