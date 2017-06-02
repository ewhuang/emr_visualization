### Author: Edward Huang

from csv import reader
import json

### This script calls on LONI USC Parkinson's disease data.

data_folder = './data/parkinsons_loni'

def read_updrs_file(fname, updrs_dct):
    '''
    Reads the scores for all attributes in the MDS-UPDRS file.
    '''
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

def get_updrs_dct():
    '''
    Returns a dictionary mapping patients to the summed UPDRS scores.
    '''
    updrs_dct = {}
    # Sum up the scores to compute a label for each patient.
    for fname in ('MDS_UPDRS_Part_I', 'MDS_UPDRS_Part_I__Patient_Questionnaire',
        'MDS_UPDRS_Part_II__Patient_Questionnaire',
        'MDS_UPDRS_Part_III__Post_Dose_', 'MDS_UPDRS_Part_IV'):
        read_updrs_file(fname, updrs_dct)
    return updrs_dct

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

def read_test_score(test_type):
    '''
    Returns a dictionary mapping PATNOs to a variety of test results.
    Includes Benton, letter-number sequencing, modified Schwab-England, and
    MoCA scores.
    Key: PATNO -> str
    Value: Derived-MOANS (Age and Education) score -> float
    Not to be used with ProSNet.
    '''
    assert test_type in ('benton', 'lns', 'schwab', 'montreal', 'semantic',
        'symbol')
    if test_type == 'benton':
        fname = 'Benton_Judgment_of_Line_Orientation'
        test_name = 'DVS_JLO_MSSAE'
    elif test_type == 'lns':
        fname = 'Letter_-_Number_Sequencing__PD_'
        test_name = 'DVS_LNS'
    elif test_type == 'schwab':
        fname = 'Modified_Schwab_+_England_ADL'
        test_name = 'MSEADLG'
    elif test_type == 'montreal':
        fname = 'Montreal_Cognitive_Assessment__MoCA_'
        test_name = 'MCATOT'
    elif test_type == 'semantic':
        fname = 'Semantic_Fluency'
        test_name = 'DVS_SFTANIM'
    elif test_type == 'symbol':
        fname = 'Symbol_Digit_Modalities'
        test_name = 'DVSD_SDM'

    test_dct = {}

    feat_name_lst = ['PATNO', 'EVENT_ID', test_name]
    f = open('%s/%s.csv' % (data_folder, fname), 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, score = [line[feat_idx] for feat_idx in feat_idx_lst]
        # Skip patients without UPDRS scores or non-baseline visits.
        if event_id != 'BL' or score == '':
            continue
        # patno, score = line[patno_idx], line[score_idx]
        assert patno not in test_dct or patno == '54186'
        test_dct[patno] = [(test_name, float(score))]
    f.close()
    return test_dct, set([test_name])

def read_test_analysis(test_type):
    '''
    Returns a dictionary mapping PATNOs to test results.
    Key: PATNO -> str
    Value: Another dictionary mapping test names to test values.
    '''
    assert test_type in ['biospecimen', 'concom_medications']
    test_dct, feature_set = {}, set([])

    if test_type == 'biospecimen':
        feat_name_lst = ('PATNO', 'CLINICAL_EVENT', 'TESTNAME', 'TESTVALUE',
            'UNITS')
        fname = 'Biospecimen_Analysis_Results'
    elif test_type == 'concom_medications':
        feat_name_lst = ('PATNO', 'EVENT_ID', 'WHODRUG', 'CMDOSE', 'CMDOSU')
        fname = 'Concomitant_Medications'
        # Load the dictionary mapping WHO terms to STITCH identifiers.
        with open('./data/who_to_stitch_dct.json', 'r') as fp:
            who_to_stitch_dct = json.load(fp)
    # Make sure that the units are consistent across all tests.
    test_unit_dct = {}
    f = open('%s/%s.csv' % (data_folder, fname), 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, test_name, test_value, units = (line[feat_idx] for
            feat_idx in feat_idx_lst)
        if test_type == 'concom_medications':
            if test_name not in who_to_stitch_dct:
                continue
            test_name = who_to_stitch_dct[test_name]
        # Only biospecimen data contains BL visits.
        if test_type == 'biospecimen' and event_id != 'BL':
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
        if units != test_unit_dct[test_name] or test_value == 0:
            continue
        # Update the test results dictionary with the patient.
        if patno not in test_dct:
            test_dct[patno] = {}
        # Skip tests that already occur for the patient.
        if test_name in test_dct[patno]:
            continue
        test_dct[patno][test_name] = test_value
        feature_set.add(test_name)
    f.close()
    # Convert each patient's dictionary to a set of tuples.
    for patno in test_dct:
        test_dct[patno] = test_dct[patno].items()
    return test_dct, feature_set

def read_clinical_diagnosis(code_dct):
    '''
    Returns the major clinical diagnosis and additional notes for each patient.
    '''
    clinical_diag_dct, feature_set = {}, set([])

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
        # Update the patient in the dictionary.
        assert patno not in clinical_diag_dct
        clinical_diag_dct[patno] = []
        # Update the patient's binary column features.
        test_val_lst = [line[test_idx] for test_idx in test_idx_lst]
        for test_name_idx, test_val in enumerate(test_val_lst):
            if test_val == '1':
                test_name = test_name_lst[test_name_idx]
                clinical_diag_dct[patno] += [(test_name, 1)]
                feature_set.add(test_name)

        # Skip "Other neurological disorder(s) (specify)".
        if prim_diag == '97':
            continue
        # Decode the primary diagnosis.
        prim_diag = code_dct[pag_name]['PRIMDIAG'][prim_diag]
        clinical_diag_dct[patno] += [(prim_diag, 1)]
        feature_set.add(prim_diag)
    f.close()
    return clinical_diag_dct, feature_set

# def read_cognitive_assessments():
#     '''
#     Returns a dictionary mapping PATNOs to cognitive assessment times.
#     '''
#     cognitive_assessment_dct = {}

#     feat_name_lst = ('HVLTRTTM', 'HVLTDRTM', 'LNORNTTM', 'SFTTM', 'LNSTM',
#         'SDMTM')
#     f = open('%s/Cognitive_Assessments.csv' % data_folder, 'r')
#     it = reader(f)
#     # Process the header line.
#     header = it.next()
#     patno_idx, event_idx = header.index('PATNO'), header.index('EVENT_ID')
#     feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
#     for line in it:
#         patno, event_id = line[patno_idx], line[event_idx]
#         if patno not in updrs_dct or event_id != 'BL':
#             continue
#         assert patno not in cognitive_assessment_dct
#         cognitive_assessment_dct[patno] = []

#         feat_val_lst = [line[feat_idx] for feat_idx in feat_idx_lst]
#         # Update patient's cognitive assessments.
#         for feat_name_idx, feat_val in enumerate(feat_val_lst):
#             if feat_val == '':
#                 continue
#             test_name = feat_name_lst[feat_name_idx]
#             test_time = sum(x * int(t) for x, t in zip([60, 1, 0.01],
#                 feat_val.split(":")))
#             cognitive_assessment_dct[patno] += [(test_name, test_time)]
#     f.close()
#     return cognitive_assessment_dct

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
        if event_id != 'BL' or cogdxcl != '1':
            continue
        assert patno not in cognitive_categorization_dct
        cognitive_categorization_dct[patno] = []
        if decline == '1':
            cognitive_categorization_dct[patno] += [('COGDECLN', 1)]
        if fncdtcog == '1':
            cognitive_categorization_dct[patno] += [('FNCDTCOG', 1)]
        cognitive_categorization_dct[patno] += [('COGSTATE', int(cog_state))]
    f.close()
    return cognitive_categorization_dct, ['COGDECLN', 'FNCDTCOG', 'COGSTATE']

def read_medical_conditions():
    '''
    Maps PATNOs to current medical conditions.
    '''
    medical_condition_dct, feature_set = {}, set([])
    f = open('%s/Current_Medical_Conditions_Log.csv' % data_folder, 'r')
    for i, line in enumerate(reader(f)):
        # Process header line.
        if i == 0:
            patno_idx = line.index('PATNO')
            term_idx = line.index('PT_NAME')
            resolved_idx = line.index('RESOLVD')
            continue
        patno, resolved = line[patno_idx], line[resolved_idx]
        if patno not in medical_condition_dct:
            medical_condition_dct[patno] = []
        term = line[term_idx]
        if term == '':
            continue
        medical_condition_dct[patno] += [(term, 1)]
        feature_set.add(term)
    f.close()
    return medical_condition_dct, feature_set

def read_epworth_scale():
    '''
    Maps PATNOs to epworth sleepiness scales. These are not to be used as
    features in the ProSNet network.
    '''
    epworth_dct, feature_set = {}, set([])

    feat_name_lst = ['ESS' + str(i) for i in range(1, 9)]
    f = open('%s/Epworth_Sleepiness_Scale.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    patno_idx, event_idx = header.index('PATNO'), header.index('EVENT_ID')
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id = line[patno_idx], line[event_idx]
        # Skip non-baseline visits.
        if event_id != 'BL':
            continue
        feat_val_lst = [line[feat_idx] for feat_idx in feat_idx_lst]
        # Update the dictionary with the current patient.
        assert patno not in epworth_dct
        epworth_dct[patno] = []
        for feat_name_idx, feat_val in enumerate(feat_val_lst):
            if feat_val == '':
                feat_val = '0'
            feat_name = feat_name_lst[feat_name_idx]
            epworth_dct[patno] += [(feat_name, float(feat_val))]
            feature_set.add(feat_name)
    f.close()

    return epworth_dct, feature_set

def read_family_history():
    '''
    Maps PATNOs to family history of PD. These are not to be used as features
    in the ProSNet network.
    '''
    family_members = ('BIOMOM', 'BIODAD', 'FULSIB', 'HAFSIB', 'MAGPAR',
        'PAGPAR', 'MATAU', 'PATAU', 'KIDSNUM')

    family_history_dct, feature_set = {}, set([])
    f = open('%s/Family_History__PD_.csv' % data_folder, 'r')
    it = reader(f)
    header = it.next()
    patno_idx = header.index('PATNO')
    family_idx_lst = [header.index(rel) for rel in family_members]
    for i, line in enumerate(reader(f)):
        patno = line[patno_idx]

        # PATNO 54186 has duplicate entries of BL.
        assert patno not in family_history_dct or patno == '54186'
        family_history_dct[patno] = []
        for name_idx, value_idx in enumerate(family_idx_lst):
            # Skip relatives that have 0 in either numerator or denominator.
            num_total_rel = line[value_idx]
            num_rel_pd = line[value_idx + 1]
            if num_total_rel in ['', '0'] or num_rel_pd in ['', '0']:
                continue
            num_total_rel = float(num_total_rel)
            num_rel_pd = float(num_rel_pd)

            # Update the fraction of relatives with PD.
            relative = family_members[name_idx]
            family_history_dct[patno] += [(relative, num_rel_pd /
                num_total_rel)]
            feature_set.add(relative)
    f.close()
    return family_history_dct, feature_set

# def read_rbd():
#     '''
#     Returns a dictionary mapping PATNOs to the drugs they are using.
#     '''
#     drug_lst = ('ONCLNZP', 'ONBENZ', 'ONMLATON', 'ONSSRI', 'ONNORSRI',
#         'ONTRIADP', 'ONBTABLK')
#     rbd_dct = {}
#     f = open('%s/Features_of_REM_Behavior_Disorder.csv' % data_folder, 'r')
#     for i, line in enumerate(reader(f)):
#         # Process the header line.
#         if i == 0:
#             patno_idx = line.index('PATNO')
#             drug_idx_lst = [line.index(drug) for drug in drug_lst]
#             continue
#         patno = line[patno_idx]
#         for drug_idx, col_idx in enumerate(drug_idx_lst):
#             onDrug = line[col_idx]
#             # Skip drug if patient is not taking it.
#             if onDrug == '0':
#                 continue
#             if patno not in rbd_dct:
#                 rbd_dct[patno] = {}
#             drug = drug_lst[drug_idx]
#             rbd_dct[patno][drug] = [1]
#     f.close()
#     return rbd_dct

def read_binary_tests(exam_type):
    '''
    Maps PATNOs to their exam results, whether there are any abnormalities.
    '''
    assert exam_type in ['neuro', 'pd_features', 'rem_disorder', 'medication']
    binary_test_dct, feature_set = {}, set([])

    if exam_type == 'neuro':
        exam_name_lst = ('MSRARSP', 'MSLARSP', 'MSRLRSP', 'MSLLRSP', 'COFNRRSP',
            'COFNLRSP', 'COHSRRSP', 'COHSLRSP', 'SENRARSP', 'SENLARSP',
            'SENRLRSP', 'SENLLRSP', 'RFLRARSP', 'RFLLARSP', 'RFLRLRSP',
            'RFLLLRSP', 'PLRRRSP', 'PLRLRSP')
        fname = 'General_Neurological_Exam'
    elif exam_type == 'pd_features':
        exam_name_lst = ('DXTREMOR', 'DXRIGID', 'DXBRADY', 'DXPOSINS')
        fname = 'PD_Features'
    elif exam_type == 'rem_disorder':
        exam_name_lst = ('DRMVIVID', 'DRMAGRAC', 'DRMNOCTB', 'SLPLMBMV',
            'SLPINJUR', 'DRMVERBL', 'DRMFIGHT', 'DRMUMV', 'DRMOBJFL',
            'MVAWAKEN', 'DRMREMEM', 'SLPDSTRB', 'STROKE', 'HETRA', 'PARKISM',
            'RLS', 'NARCLPSY', 'DEPRS', 'EPILEPSY', 'BRNINFM')
        fname = 'REM_Sleep_Disorder_Questionnaire'
    elif exam_type == 'medication':
        exam_name_lst = ['ONLDOPA', 'ONDOPAG']
        fname = 'Use_of_PD_Medication'

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
        assert patno not in binary_test_dct or patno == '54186'
        binary_test_dct[patno] = []
        # Get abnormal values ('1') for each test.
        exam_val_lst = [line[exam_idx] for exam_idx in exam_idx_lst]
        for exam_name_idx, exam_val in enumerate(exam_val_lst):
            if exam_val == '1':
                exam_name = exam_name_lst[exam_name_idx]
                binary_test_dct[patno] += [(exam_name, 1)]
                feature_set.add(exam_name)
    f.close()
    return binary_test_dct, feature_set

def read_hvlt():
    '''
    Maps PATNOs to Hopkins Test dictionary. Not to be used with ProSNet.
    '''
    hvlt_dct, feature_set = {}, set([])

    score_name_lst = ('DVT_TOTAL_RECALL', 'DVT_DELAYED_RECALL', 'DVT_RETENTION',
        'DVT_RECOG_DISC_INDEX')
    f = open('%s/Hopkins_Verbal_Learning_Test.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    score_idx_lst = [header.index(score_name) for score_name in score_name_lst]
    patno_idx, event_idx = header.index('PATNO'), header.index('EVENT_ID')
    comm_idx, age_idx = header.index('comm'), header.index('AGE_ASSESS_HVLT')
    for line in it:
        patno, event_id = line[patno_idx], line[event_idx]
        if event_id != 'BL':
            continue
        # Skip patients that have comments or do not have age recorded.
        comment, age = line[comm_idx], line[age_idx]
        if comment != '' or age == '':
            continue

        assert patno not in hvlt_dct
        hvlt_dct[patno] = []
        # Update patient with scores.
        score_lst = map(float, [line[score_idx] for score_idx in score_idx_lst])
        for score_name_idx, score in enumerate(score_lst):
            score_name = score_name_lst[score_name_idx]
            hvlt_dct[patno] += [(score_name, score)]
            feature_set.add(score_name)
    f.close()

    return hvlt_dct, feature_set

def read_demographics():
    demographics_dct, feature_set = {}, set([])

    gender_dct = {'0':'Female of child bearing potential',
        '1':'Female of non-child bearing potential', '2':'Male'}
    race_lst = ('HISPLAT', 'RAINDALS', 'RAASIAN', 'RABLACK', 'RAHAWOPI',
        'RAWHITE')
    f = open('%s/Screening___Demographics.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    race_idx_lst = [header.index(race) for race in race_lst]
    patno_idx, gender_idx = header.index('PATNO'), header.index('GENDER')
    for line in it:
        patno, gender = line[patno_idx], line[gender_idx]
        # Skip patients with no gender information.
        if gender == '':
            continue
        gender = gender_dct[gender]
        # Update the patient.
        if patno in demographics_dct:
            continue
        # Initialize patient with gender.
        demographics_dct[patno] = [(gender, 1)]
        feature_set.add(gender)
        # Add race information.
        race_val_lst = [line[race_idx] for race_idx in race_idx_lst]
        for race_idx, race_val in enumerate(race_val_lst):
            if race_val == '1':
                race = race_lst[race_idx]
                demographics_dct[patno] += [(race, 1)]
                feature_set.add(race)
    f.close()
    return demographics_dct, feature_set

def read_pd_surgery():
    pd_surgery_dct = {}

    feat_name_lst = ('PATNO', 'EVENT_ID', 'PDSURG', 'PDSURGTP', 'PDSLGPI',
        'PDSLSTN')
    pdsurg_type_dct = {'1':'DBS (Deep Brain Stimulation)',
        '2':'Levodopa intestinal gel infusion'}
    feature_set = set(pdsurg_type_dct.values())
    f = open('%s/Surgery_for_Parkinson_Disease.csv' % data_folder, 'r')
    it = reader(f)
    # Process the header line.
    header = it.next()
    feat_idx_lst = [header.index(feat) for feat in feat_name_lst]
    for line in it:
        patno, event_id, pdsurg, pdsurgtp, pdslgpi, pdslstn = (line[feat_idx]
            for feat_idx in feat_idx_lst)
        if event_id != 'BL' or pdsurg != '1':
            continue
        # Skip surgery of type "other".
        if pdsurgtp == '3':
            continue
        assert patno not in pd_surgery_dct

        pd_surgery_dct[patno] = [(pdsurg_type_dct[pdsurgtp], 1)]

        if pdslgpi == '1':
            pd_surgery_dct[patno] += [('PDSLGPI', 1)]
        if pdslstn == '1':
            pd_surgery_dct[patno] += [('PDSLSTN', 1)]
    f.close()
    feature_set.add('PDSLGPI')
    feature_set.add('PDSLSTN')

    return pd_surgery_dct, feature_set

def read_mutation_file():
    '''
    Reads the PPMI mutation file, and records the genes with mutations for each
    patient. No header line.
    '''
    mutation_dct, feature_set = {}, set([])
    # f = open('./data/PPMI_mutation.txt', 'r')
    f = open('./data/PPMI_indels.txt', 'r')
    for line in f:
        patno, mutated_gene, att_1, att_2 = line.strip().split('\t')
        # Add the mutation tag line to the gene.
        mutated_gene += ' MUT'
        if patno not in mutation_dct:
            mutation_dct[patno] = set([])
        mutation_dct[patno].add((mutated_gene, 1))
        feature_set.add(mutated_gene)
    f.close()
    return mutation_dct, feature_set

def read_patient_status():
    '''
    Reads the patient status for each patient number.
    '''
    status_dct = {}
    f = open('%s/Biospecimen_Analysis_Results.csv' % data_folder, 'r')
    f.readline()
    for line in reader(f):
        patno, diagnosis = line[0], line[2]
        if diagnosis not in ['PD', 'SWEDD', 'Control'] or patno in status_dct:
            continue
        status_dct[patno] = diagnosis
    f.close()
    return status_dct

def main():
    updrs_dct = get_updrs_dct()
    status_dct = read_patient_status()

    biospecimen_dct = read_test_analysis('biospecimen')[0]
    concom_medication_dct = read_test_analysis('concom_medications')[0]

    code_dct = read_code_file()
    clinical_diag_dct = read_clinical_diagnosis(code_dct)[0]

    cognitive_categorization_dct = read_cognitive_categorizations()[0]
    medical_condition_dct = read_medical_conditions()[0]

    neuro_exam_dct = read_binary_tests('neuro')[0]
    pd_feat_dct = read_binary_tests('pd_features')[0]
    rem_disorder_dct = read_binary_tests('rem_disorder')[0]

    demographics_dct = read_demographics()[0]
    pd_surgery_dct = read_pd_surgery()[0]
    pd_medication_dct = read_binary_tests('medication')[0]
    mutation_dct = read_mutation_file()[0]

    # This block not to be used in ProSNet network.
    line_orientation_dct = read_test_score('benton')[0]
    epworth_dct = read_epworth_scale()[0]
    family_history_dct = read_family_history()[0]
    hvlt_dct = read_hvlt()[0]
    lns_dct = read_test_score('lns')[0]
    schwab_dct = read_test_score('schwab')[0]
    moca_dct = read_test_score('montreal')[0]
    semantic_dct = read_test_score('semantic')[0]
    symbol_dct = read_test_score('symbol')[0]
    # End block not to be used in ProSNet network.

    print 'Running tests...'
    # Test the UPDRS dictionary.
    assert updrs_dct['3400'] == 51
    assert '3210' not in updrs_dct

    # Test the patient status dictionary.
    assert status_dct['3538'] == 'SWEDD'
    assert '53286' not in status_dct
    assert status_dct['3401'] == 'Control'
    assert status_dct['3314'] == 'PD'
    
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

    # Test concomitant medications dictionary.
    assert 'estrogen' not in concom_medication_dct['50157']
    assert ('acyclovir', 500) in concom_medication_dct['51551']
    assert ('aripiprazole', 2) in concom_medication_dct['52373']
    assert ('gabapentin', 600) in concom_medication_dct['4027']

    # # Test hematology dictionary.
    # assert hematology_dct['3252']['Calcium (EDTA)'] == [2.57, 2.42, 2.45, 2.42,
    #     2.5, 2.57, 2.5]
    # assert hematology_dct['3276']['Serum Sodium'] == [139, 139, 141, 138, 139]
    # assert hematology_dct['10874']['Total Protein'] == [63, 73, 67]

    # Test the clinical diagnosis and managemenet dictionary.
    assert clinical_diag_dct['3465'] == [('DCRTREM', 1), ('DCRIGID', 1
        ), ('DCBRADY', 1), ('Idiopathic PD', 1)]
    assert '3082' not in clinical_diag_dct
    assert '3326' not in clinical_diag_dct
    assert clinical_diag_dct['3836'] == [('DCRTREM', 1), ('DCRIGID', 1
        ), ('DCBRADY', 1), ('Idiopathic PD', 1)]

    # # Test cognitive assessments dictionary.
    # assert '3154' not in cognitive_assessment_dct
    # assert '3951' not in cognitive_assessment_dct
    # assert ('LNSTM', 685) in cognitive_assessment_dct['41749']
    
    # Test cognitive categorizations dictionary.
    assert cognitive_categorization_dct['3026'] == [('COGSTATE', 1)]
    assert '40709' not in cognitive_categorization_dct
    assert cognitive_categorization_dct['40755'] == [('COGDECLN', 1), (
        'FNCDTCOG', 1), ('COGSTATE', 2)]

    # # Test medical conditions dictionary.
    assert ('Urinary incontinence', 1) in medical_condition_dct['3001']
    assert ('Drug hypersensitivity', 1) in medical_condition_dct['3008']
    assert ('Postmenopause', 1) in medical_condition_dct['3453']

    # Test the Epworth sleepiness scale dictionary.
    assert ('ESS1', 2) in epworth_dct['3000']
    assert ('ESS2', 1) in epworth_dct['3000']
    assert ('ESS5', 2) in epworth_dct['3009']
    assert ('ESS5', 0) in epworth_dct['3636']
    assert ('ESS4', 3) in epworth_dct['60065']

    # Test family history dictionary.
    assert ('MAGPAR', 0.5) in family_history_dct['3101']
    assert ('KIDSNUM', 1.0 / 6) in family_history_dct['3653']
    assert ('HAFSIB', 3.0 / 9) in family_history_dct['52620']

    # # Test REM behavior disorder dictionary.
    # assert rbd_dct['60033'] == {'ONMLATON':[1], 'ONSSRI':[1], 'ONNORSRI':[1]}
    # assert rbd_dct['60006'] == {'ONCLNZP':[1]}

    # Test the general neurological exam dictionary.
    assert ('COFNLRSP', 1) in neuro_exam_dct['42449']
    assert ('MSRARSP', 0) not in neuro_exam_dct['50813']
    assert ('MSRARSP', 1) not in neuro_exam_dct['50813']

    # Test Hopkins verbal learning test dictionary.
    assert ('DVT_TOTAL_RECALL', 37) in hvlt_dct['3502']
    assert ('DVT_RETENTION', 55) in hvlt_dct['3010']
    assert ('DVT_DELAYED_RECALL', 32) in hvlt_dct['56680']

    # Test the letter-number sequencing dictionary.
    assert lns_dct['3400'] == [('DVS_LNS', 9)]
    assert lns_dct['40594'] == [('DVS_LNS', 5)]
    assert lns_dct['55124'] == [('DVS_LNS', 14)]

    # Test the Schwab-England dictionary.
    assert schwab_dct['40612'] == [('MSEADLG', 70)]
    assert schwab_dct['51918'] == [('MSEADLG', 100)]
    assert schwab_dct['41395'] == [('MSEADLG', 20)]

    # Test the Montreal cognitive assessment dictionary.
    assert moca_dct['40612'] == [('MCATOT', 15)]
    assert moca_dct['51918'] == [('MCATOT', 24)]
    assert moca_dct['41395'] == [('MCATOT', 12)]

    # Test the PD feature dictionary.
    assert pd_feat_dct['40551'] == [('DXTREMOR', 1), ('DXBRADY', 1)]
    assert pd_feat_dct['40938'] == [('DXTREMOR', 1), ('DXRIGID', 1), ('DXBRADY',
        1), ('DXPOSINS', 1)]
    assert pd_feat_dct['40753'] == [('DXRIGID', 1), ('DXBRADY', 1)]

    # Test the REM disorder dictionary.
    assert ('DRMVIVID', 1) in rem_disorder_dct['3401']
    assert ('DRMFIGHT', 1) in rem_disorder_dct['3400']
    assert ('HETRA', 1) in rem_disorder_dct['3451']
    assert ('SLPLMBMV', 1) in rem_disorder_dct['54186']

    # Test the demographics dictionary.
    assert demographics_dct['3400'] == [('Female of child bearing potential',
        1), ('RAWHITE', 1)]

    # Test semantic fluency dictionary.
    assert semantic_dct['3411'] == [('DVS_SFTANIM', 19)]
    assert semantic_dct['3501'] == [('DVS_SFTANIM', 7)]
    assert semantic_dct['41985'] == [('DVS_SFTANIM', 8)]

    # Test PD surgery dictionary.
    assert '3410' not in pd_surgery_dct
    assert '50418' not in pd_surgery_dct
    assert pd_surgery_dct['52599'] == [('DBS (Deep Brain Stimulation)',
        1), ('PDSLSTN', 1)]

    # Test the symbol digit modalities dictionary.
    assert symbol_dct['3403'] == [('DVSD_SDM', -1)]
    assert symbol_dct['50621'] == [('DVSD_SDM', 1.25)]
    assert symbol_dct['42164'] == [('DVSD_SDM', -0.167)]

    # Test the PD medication dictionary.
    assert pd_medication_dct['40595'] == [('ONLDOPA', 1)]
    assert pd_medication_dct['41400'] == [('ONLDOPA', 1), ('ONDOPAG', 1)]
    assert pd_medication_dct['41288'] == [('ONDOPAG', 1)]

    # Test the PPMI mutation dictionary.
    assert ('RAB29 MUT', 1) in mutation_dct['3800']
    assert ('MAPT MUT', 1) in mutation_dct['3004']

    print 'Finished tests!'

if __name__ == '__main__':
    main()