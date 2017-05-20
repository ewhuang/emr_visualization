### Author: Edward Huang

from csv import reader

### This script calls on LONI USC Parkinson's disease data.

data_folder = './data/parkinsons_loni'
updrs_dct = {} # Dictionary mapping patients to the summed UPDRS scores.

# TODO.
def read_updrs_file(fname):
    '''
    Reads the scores for all attributes in the MDS-UPDRS file.
    '''
    global updrs_dct
    current_patient_set = set([])
    f = open('%s/%s.csv' % (data_folder, fname), 'r')
    for i, line in enumerate(reader(f)):
        if i == 0:
            patno_idx = line.index('PATNO')
            event_idx = line.index('EVENT_ID')
            end_idx = line.index('ORIG_ENTRY')
            if 'Part_III' in fname:
                start_idx = line.index('EXAMTM')
                end_idx = line.index('ANNUAL_TIME_BTW_DOSE_NUPDRS')
            elif 'Part_IV' in fname:
                start_idx = line.index('INFODT')
            else:
                start_idx = line.index('NUPSOURC')
            continue
        event_id, patno = line[event_idx], line[patno_idx]
        # Only use the baseline visits. Also skip duplicate entries.
        if event_id != 'BL' or patno in current_patient_set:
            continue
        # Get the score list between the start and end indices.
        score_list = line[start_idx + 1:end_idx]
        score_list = [int(score) for score in score_list if score != '' ]
        # Update the score dictionary.
        if patno not in updrs_dct:
            updrs_dct[patno] = 0
        updrs_dct[patno] += sum(score_list)
        # Update the list of patients seen in this spreadsheet.
        current_patient_set.add(patno)
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

# def read_adverse_events():
#     '''
#     Returns a dictionary mapping patient IDs to their adverse events.
#     '''
#     adverse_event_dct = {}
#     f = open('%s/Adverse_Event_Log.csv' % data_folder, 'r')
#     for i, line in enumerate(reader(f)):
#         # Process header line.
#         if i == 0:
#             patno_idx = line.index('PATNO')
#             term_idx = line.index('AETERM')
#             severity_idx = line.index('AESEVER')
#             related_idx = line.index('AERELAT')
#             continue
#         isRelated = line[related_idx]
#         # TODO: Skip adverse events that are unrelated to the study.
#         if isRelated not in ['4', '5']:
#             continue
#         patno = line[patno_idx]
#         term = line[term_idx]
#         severity = float(line[severity_idx])
#         # Create nested list mapping adverse events to their severities.
#         if patno not in adverse_event_dct:
#             adverse_event_dct[patno] = {}
#         # Create list of the adverse event's severities.
#         if term not in adverse_event_dct[patno]:
#             adverse_event_dct[patno][term] = []
#         adverse_event_dct[patno][term] += [severity]
#     f.close()
#     return adverse_event_dct

def read_line_orientation():
    '''
    Returns a dictionary mapping PATNOs to line orientation test results.
    Key: PATNO -> str
    Value: Derived-MOANS (Age and Education) score -> float
    '''
    line_orientation_dct = {}

    feat_name_lst = ['PATNO', 'EVENT_ID', 'DVS_JLO_MSSAE']
    f = open('%s/Benton_Judgment_of_Line_Orientation.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, score = [line[feat_idx] for feat_idx in feat_idx_lst]
        # Skip patients without UPDRS scores or non-baseline visits.
        if patno not in updrs_dct or event_id != 'BL' or score == '':
            continue
        # patno, score = line[patno_idx], line[score_idx]
        assert patno not in line_orientation_dct
        line_orientation_dct[patno] = [('DVS_JLO_MSSAE', float(score))]
    f.close()

    return line_orientation_dct

def read_biospecimen_analysis():
    '''
    Returns a dictionary mapping PATNOs to test results.
    Key: PATNO -> str
    Value: Another dictionary mapping test names to test values.
    '''
    biospecimen_dct = {}

    feat_name_lst = ['PATNO', 'CLINICAL_EVENT', 'TESTNAME', 'TESTVALUE', 'UNITS']
    # Make sure that the units are consistent across all tests.
    test_unit_dct = {}
    f = open('%s/Biospecimen_Analysis_Results.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, test_name, test_value, units = (line[feat_idx] for
            feat_idx in feat_idx_lst)
        # Skip patients without UPDRS scores or non-baseline visits.
        if patno not in updrs_dct or event_id != 'BL':
            continue
        # Attempt converting test value to float.
        try:
            test_value = float(test_value)
        except ValueError:
            continue
        # Update the units for current test. Use the first units we encounter.
        if test_name not in test_unit_dct:
            test_unit_dct[test_name] = units
        # Skip tests that are not consistent in units.
        if units != test_unit_dct[test_name]:
            continue
        # Update the test results dictionary with the patient.
        if patno not in biospecimen_dct:
            biospecimen_dct[patno] = {}
        biospecimen_dct[patno][test_name] = test_value
    f.close()
    # Convert each patient's dictionary to a set of tuples.
    for patno in biospecimen_dct:
        biospecimen_dct[patno] = biospecimen_dct[patno].items()
    return biospecimen_dct

def read_clinical_diagnosis(code_dct):
    '''
    Returns the major clinical diagnosis and additional notes for each patient.
    '''
    clinical_diag_dct = {}

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
        if patno not in updrs_dct or event_id != 'BL' or ps_lvl != '1':
            continue
        # Update the patient in the dictionary.
        assert patno not in clinical_diag_dct
        clinical_diag_dct[patno] = []
        # Update the patient's binary column features.
        test_val_lst = [line[test_idx] for test_idx in test_idx_lst]
        for test_name_idx, test_val in enumerate(test_val_lst):
            if test_val == '1':
                clinical_diag_dct[patno] += [(test_name_lst[test_name_idx], 1)]

        # Skip "Other neurological disorder(s) (specify)".
        if prim_diag == '97':
            continue
        # Decode the primary diagnosis.
        prim_diag = code_dct[pag_name]['PRIMDIAG'][prim_diag]
        clinical_diag_dct[patno] += [(prim_diag, 1)]
    f.close()
    return clinical_diag_dct

def read_cognitive_assessments():
    '''
    Returns a dictionary mapping PATNOs to cognitive assessment times.
    '''
    cognitive_assessment_dct = {}

    feat_name_lst = ('HVLTRTTM', 'HVLTDRTM', 'LNORNTTM', 'SFTTM', 'LNSTM',
        'SDMTM')
    f = open('%s/Cognitive_Assessments.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    patno_idx, event_idx = header.index('PATNO'), header.index('EVENT_ID')
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id = line[patno_idx], line[event_idx]
        if patno not in updrs_dct or event_id != 'BL':
            continue
        assert patno not in cognitive_assessment_dct
        cognitive_assessment_dct[patno] = []

        feat_val_lst = [line[feat_idx] for feat_idx in feat_idx_lst]
        # Update patient's cognitive assessments.
        for feat_name_idx, feat_val in enumerate(feat_val_lst):
            if feat_val == '':
                continue
            test_name = feat_name_lst[feat_name_idx]
            test_time = sum(x * int(t) for x, t in zip([60, 1, 0.01],
                feat_val.split(":")))
            cognitive_assessment_dct[patno] += [(test_name, test_time)]
    f.close()
    return cognitive_assessment_dct

def read_cognitive_categorizations():
    '''
    Returns a dictionary mapping PATNOs to cognitive categorizations.
    '''
    cognitive_categorization_dct = {}

    feat_name_lst = ('PATNO', 'EVENT_ID', 'COGDECLN', 'FNCDTCOG', 'COGSTATE',
        'COGDXCL')
    f = open('%s/Cognitive_Categorization.csv' % data_folder, 'r')
    it = reader(f)
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, decline, fncdtcog, cog_state, cogdxcl = (line[feat_idx]
            for feat_idx in feat_idx_lst)

        # Skip diagnoses with low confidence.
        if patno not in updrs_dct or event_id != 'BL' or cogdxcl != '1':
            continue
        assert patno not in cognitive_categorization_dct
        cognitive_categorization_dct[patno] = []
        if decline == '1':
            cognitive_categorization_dct[patno] += [('COGDECLN', 1)]
        if fncdtcog == '1':
            cognitive_categorization_dct[patno] += [('FNCDTCOG', 1)]
        cognitive_categorization_dct[patno] += [('COGSTATE', int(cog_state))]
    f.close()
    return cognitive_categorization_dct

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
            continue
        patno = line[patno_idx]
        if patno not in medical_condition_dct:
            medical_condition_dct[patno] = {}
        medical_condition_dct[patno][line[term_idx]] = [1]
    f.close()
    return medical_condition_dct

def read_family_history():
    '''
    Maps PATNOs to family history of PD.
    '''
    family_members = ('BIOMOM', 'BIODAD', 'FULSIB', 'HAFSIB', 'MAGPAR',
        'PAGPAR', 'MATAU', 'PATAU', 'KIDSNUM')
    family_history_dct = {}
    f = open('%s/Family_History__PD_.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            family_idx_lst = [line.index(rel) for rel in family_members]
            continue
        patno = line[patno_idx]
        if patno not in family_history_dct:
            family_history_dct[patno] = {}
        for str_idx, relative_idx in enumerate(family_idx_lst):
            # Skip relatives that have 0 in either numerator or denominator.
            num_total_rel = line[relative_idx]
            num_rel_pd = line[relative_idx + 1]
            if num_total_rel in ['', '0'] or num_rel_pd in ['', '0']:
                continue
            num_total_rel = float(num_total_rel)
            num_rel_pd = float(num_rel_pd)

            # Update the fraction of relatives with PD.
            relative = family_members[str_idx]
            if relative not in family_history_dct[patno]:
                family_history_dct[patno][relative] = []
            family_history_dct[patno][relative] = [num_rel_pd / num_total_rel]
    f.close()
    return family_history_dct

def read_rbd():
    '''
    Returns a dictionary mapping PATNOs to the drugs they are using.
    '''
    drug_lst = ('ONCLNZP', 'ONBENZ', 'ONMLATON', 'ONSSRI', 'ONNORSRI',
        'ONTRIADP', 'ONBTABLK')
    rbd_dct = {}
    f = open('%s/Features_of_REM_Behavior_Disorder.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process the header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            drug_idx_lst = [line.index(drug) for drug in drug_lst]
            continue
        patno = line[patno_idx]
        for drug_idx, col_idx in enumerate(drug_idx_lst):
            onDrug = line[col_idx]
            # Skip drug if patient is not taking it.
            if onDrug == '0':
                continue
            if patno not in rbd_dct:
                rbd_dct[patno] = {}
            drug = drug_lst[drug_idx]
            rbd_dct[patno][drug] = [1]
    f.close()
    return rbd_dct

def main():
    # Sum up the scores to compute a label for each patient.
    read_updrs_file('MDS_UPDRS_Part_I')
    read_updrs_file('MDS_UPDRS_Part_I__Patient_Questionnaire')
    read_updrs_file('MDS_UPDRS_Part_II__Patient_Questionnaire')
    read_updrs_file('MDS_UPDRS_Part_III__Post_Dose_')
    read_updrs_file('MDS_UPDRS_Part_IV')

    # TODO: Currently not doing adverse events, since they don't specify BL.
    # adverse_event_dct = read_adverse_events()
    line_orientation_dct = read_line_orientation()
    biospecimen_dct = read_biospecimen_analysis()
    # TODO: No hematology, because no baseline visits in these tests.
    # hematology_dct = read_test_analysis('hematology')

    code_dct = read_code_file()
    clinical_diag_dct = read_clinical_diagnosis(code_dct)

    cognitive_assessment_dct = read_cognitive_assessments()

    cognitive_categorization_dct = read_cognitive_categorizations()

    # TODO: no medication, because no baseline visits.
    # medication_dct = read_test_analysis('medication')

    medical_condition_dct = read_medical_conditions()

    family_history_dct = read_family_history()

    rbd_dct = read_rbd()

    ### TESTS___________________________________________________________________
    # Test the UPDRS dictionary.

    # # Test the adverse event dictionary.
    # assert adverse_event_dct['3226'] == {'BACK SORENESS':[1.0, 2.0]}
    # assert adverse_event_dct['52215'] == {'SORENESS AT LP SITE':[1.0]}
    # assert adverse_event_dct['60036'] == {'HEADACHE':[2.0, 1.0, 1.0, 3.0, 1.0],
    #     'NAUSEA':[1.0], 'SWEATING':[1.0], 'SHIVERS':[1.0]}

    # Test the line orientation test dictionary.
    assert line_orientation_dct['3400'] == [('DVS_JLO_MSSAE', 11.7)]
    assert line_orientation_dct['3552'] == [('DVS_JLO_MSSAE', 13.72)]
    assert line_orientation_dct['41412'] == [('DVS_JLO_MSSAE', 11.52)]

    # Test the biospecimen dictionary.
    assert ('Abeta 42', 310.8) in biospecimen_dct['3003']
    assert ('CSF Alpha-synuclein', 1520.36) in biospecimen_dct['3018']
    assert ('CSF Alpha-synuclein', 1036.57) in biospecimen_dct['4077']

    # # Test hematology dictionary.
    # assert hematology_dct['3252']['Calcium (EDTA)'] == [2.57, 2.42, 2.45, 2.42,
    #     2.5, 2.57, 2.5]
    # assert hematology_dct['3276']['Serum Sodium'] == [139, 139, 141, 138, 139]
    # assert hematology_dct['10874']['Total Protein'] == [63, 73, 67]

    # Test the clinical diagnosis and managemenet dictionary.
    assert clinical_diag_dct['3465'] == [('DCRTREM', 1), ('DCRIGID', 1), ('DCBRADY', 1), ('Idiopathic PD', 1)]
    assert '3082' not in clinical_diag_dct
    assert '3326' not in clinical_diag_dct
    assert clinical_diag_dct['3836'] == [('DCRTREM', 1), ('DCRIGID', 1), ('DCBRADY', 1), ('Idiopathic PD', 1)]

    # Test cognitive assessments dictionary.
    assert '3154' not in cognitive_assessment_dct
    assert '3951' not in cognitive_assessment_dct
    assert ('LNSTM', 685) in cognitive_assessment_dct['41749']
    
    # Test cognitive categorizations dictionary.
    assert cognitive_categorization_dct['3026'] == [('COGSTATE', 1)]
    assert '40709' not in cognitive_categorization_dct
    assert cognitive_categorization_dct['40755'] == [('COGDECLN', 1), (
        'FNCDTCOG', 1), ('COGSTATE', 2)]

    # # Test concomitant medications dictionary.
    # assert '"ESTROGEN GEL"' not in medication_dct['50157']
    # assert medication_dct['51551']['ACYCLOVIR'] == [500, 400]
    # assert medication_dct['52373']['ABILIFY'] == [2, 5]

    # Test medical conditions dictionary.
    assert medical_condition_dct['3001']['Urinary incontinence'] == [1]
    assert medical_condition_dct['3008']['Drug hypersensitivity'] == [1]

    # Test family history dictionary.
    assert family_history_dct['3101']['MAGPAR'] == [0.5]
    assert family_history_dct['3653']['KIDSNUM'] == [1.0 / 6]
    assert family_history_dct['52620']['HAFSIB'] == [3.0 / 9]

    # Test REM behavior disorder dictionary.
    assert rbd_dct['60033'] == {'ONMLATON':[1], 'ONSSRI':[1], 'ONNORSRI':[1]}
    assert rbd_dct['60006'] == {'ONCLNZP':[1]}

    print 'Finished tests!'

if __name__ == '__main__':
    main()