### Author: Edward Huang

from csv import reader
from datetime import datetime

### This script calls on the MIMIC III data files and extracts the relevant
### information.
### Files to be read from: ADMISSIONS, CHARTEVENTS(?), DIAGNOSES_ICD,
### PROCEDURES_ICD, 

def get_first_visit_ids():
    '''
    Opens ADMISSIONS.csv and gets the first visit HADM_ID for each patient.
    Returns the list of HIDM_IDs that correspond to a first visit.
    '''
    first_visit_dct = {}
    hadm_to_date_dct = {}
    f = open('./data/mimic/ADMISSIONS.csv', 'r')
    f.readline() # Skip header line
    for line in reader(f):
        assert len(line) == 19
        subject_id, hadm_id = line[1:3]
        current_date = datetime.strptime(line[3], '%Y-%m-%d %H:%M:%S')
        # Map visit ID to date.
        assert hadm_id not in hadm_to_date_dct
        hadm_to_date_dct[hadm_id] = current_date
        # Add subject to dictionary.
        if subject_id not in first_visit_dct:
            first_visit_dct[subject_id] = hadm_id
        else:
            stored_hadm = first_visit_dct[subject_id]
            stored_date = hadm_to_date_dct[stored_hadm]
            # Check if stored date is actually the first date.
            if current_date < stored_date:
                first_visit_dct[subject_id] = hadm_id
    f.close()
    return set(first_visit_dct.values())

def get_icd_dct(fname):
    '''
    Opens DIAGNOSES_ICD.csv or PROCEDURES_ICD.csv, and maps the HADM_IDs to the
    diagnosis/procedure information.
    Return a dictionary mapping HADM IDs to ICD9 codes.
    Key: HADM ID -> str
    Value: ICD9 code -> str
    '''
    assert fname in ['DIAGNOSES', 'PROCEDURES']
    icd_dct = {}
    f = open('./data/mimic/%s_ICD.csv' % fname, 'r')
    f.readline() # Skip header line
    for line in reader(f):
        assert len(line) == 5
        hadm_id, icd_code = line[2], line[4]
        # Only add ICD codes for first visits.
        if hadm_id in first_visit_set:
            if hadm_id not in icd_dct:
                icd_dct[hadm_id] = set([])
            icd_dct[hadm_id].add(icd_code.strip('"'))
    f.close()
    return icd_dct

def get_lab_event_dct():
    '''
    Returns dictionary mapping HADM IDs to laboratory measurements. Only keep
    measurements with numerical values.
    Key: HADM ID -> str
    Value: list of (lab test, value) -> list((str, float))
    '''
    lab_event_dct = {}
    f = open('./data/mimic/LABEVENTS.csv', 'r')
    f.readline() # Skip header line
    for line in reader(f):
        assert len(line) == 9
        hadm_id, item_id, value = line[2], line[3], line[6]
        # Skip non-numerical features and non-first-visits.
        if value == '' or hadm_id not in first_visit_set:
            continue
        if hadm_id not in lab_event_dct:
            lab_event_dct[hadm_id] = []
        lab_event_dct[hadm_id] += [(item_id, float(value))]
    f.close()
    return lab_event_dct

def get_prescription_dct():
    prescription_dct = {}
    f = open('./data/mimic/PRESCRIPTIONS.csv', 'r')
    f.readline()
    for line in reader(f):
        assert len(line) == 19
        hadm_id, drug, dosage = line[2], line[10], line[14]
        if drug == '':
            drug = line[7]
        if hadm_id not in first_visit_set or dosage == '' or '-' in dosage:
            continue
        if ',' in dosage:
            dosage = dosage.replace(',', '')
        try:
            dosage = float(dosage)
        except ValueError:
            pass
        if hadm_id not in prescription_dct:
            prescription_dct[hadm_id] = []
        prescription_dct[hadm_id] += [(drug, dosage)]
    f.close()
    return prescription_dct

def main():
    global first_visit_set
    first_visit_set = get_first_visit_ids()

    # diagnosis_dct = get_icd_dct('DIAGNOSES')
    # procedure_dct = get_icd_dct('PROCEDURES')

    # lab_event_dct = get_lab_event_dct()
    prescription_dct = get_prescription_dct()

if __name__ == '__main__':
    main()