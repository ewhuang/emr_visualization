import numpy as np
import subprocess


# Baseline stuff.
for norm_type in ['l2', 'max']:
    command = 'python build_patient_feature_matrix.py -n %s -l updrs' % norm_type
    print command
    subprocess.call(command, shell=True)

    for pca in [100, 250, 500, 750, 1000]:
        for tsne in ['pca', 'random']:
            for learning_rate in [100, 200, 300, 400, 500]:
                for k in [5, 10, 15, 20]:                    
                    command = 'python visualize_emrs.py -n %s -t %s -p %d -r %d -k %d' % (
                        norm_type, tsne, pca, learning_rate, k)
                    print command
                    subprocess.call(command, shell=True)

# ProSNet stuff.
for norm_type in ['l2', 'max']:
    for sim_thresh in [0.01, 0.05, 0.1, 0.15]:
        command = 'python build_patient_feature_matrix.py -n %s -d 500 -s %g -w after' % (
            norm_type, sim_thresh, where_norm)
        print command
        subprocess.call(command, shell=True)

        for pca in [100, 250, 500, 750, 1000]:
            for tsne in ['pca', 'random']:
                for learning_rate in [100, 200, 300, 400, 500]:
                    for k in [5, 10, 15, 20]:                    
                        command = 'python visualize_emrs.py -n %s -d 500 -s %g -w after -t %s -p %d -r %d -k %d' % (
                            norm_type, sim_thresh, tsne, pca, learning_rate, k)
                        print command
                        subprocess.call(command, shell=True)