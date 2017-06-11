### Author: Edward Huang

from process_loni_parkinsons import *


def main():
    updrs_dct = get_updrs_dct()
    code_dct = read_code_file()

    tupe_lst = [read_cognitive_categorizations(updrs_dct),
            read_clinical_diagnosis(code_dct, updrs_dct),
            read_medical_conditions(updrs_dct)]
    feat_set = set([])
    for patient_dct, feature_lst in tupe_lst:
        feat_set = feat_set.union(feature_lst)
    
    out = open('./results/symptom_file.txt', 'w')
    out.write('\n'.join(feat_set))
    out.close()

if __name__ == '__main__':
    main()