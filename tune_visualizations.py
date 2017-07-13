import numpy as np
import subprocess


# Baseline stuff.
# for norm_type in ['l2', 'max']:
for norm_type in ['max']:
    command = 'python build_patient_feature_matrix.py -n %s -l updrs' % norm_type
    print command
    subprocess.call(command, shell=True)

    # for pca in [100, 250, 500, 750, 1000]:
    for pca in [500]:
        # for tsne in ['pca', 'random']:
        for tsne in ['pca']:
            command_lst = [] # Reset the command list.
            # for learning_rate in [100, 200, 300, 400, 500]:
            for learning_rate in [200]:
                command = 'python visualize_emrs.py -n %s -t %s -p %d -r %d' % (
                    norm_type, tsne, pca, learning_rate)
                print command
                command_lst += [command]
            processes = [subprocess.Popen(cmd, shell=True) for cmd in command_lst]
            print 'waiting...'
            for p in processes: p.wait()

# ProSNet stuff.
# for norm_type in ['l2', 'max']:
for norm_type in ['max']:
    command_lst = []
    for sim_thresh in np.arange(0.0, 0.51, 0.01):
        command = 'python build_patient_feature_matrix.py -n %s -d 500 -s %g -w after -l updrs' % (
            norm_type, sim_thresh)
        print command
        # subprocess.call(command, shell=True)
        command_lst += [command]
    processes = [subprocess.Popen(cmd, shell=True) for cmd in command_lst]
    print 'waiting...'
    for p in processes: p.wait()

# for norm_type in ['l2', 'max']:
for norm_type in ['max']:
    command_lst = []
    for sim_thresh in np.arange(0.0, 0.51, 0.01):
        # for pca in [100, 250, 500, 750, 1000]:
        for pca in [500]:
            # for tsne in ['pca', 'random']:
            for tsne in ['pca']:
                # for learning_rate in [100, 200, 300, 400, 500]:
                for learning_rate in [200]:
                    command = 'python visualize_emrs.py -n %s -d 500 -s %s -w after -t %s -p %d -r %d' % (
                        norm_type, sim_thresh, tsne, pca, learning_rate)
                    print command
                    # subprocess.call(command, shell=True)
                    command_lst += [command]
    processes = [subprocess.Popen(cmd, shell=True) for cmd in command_lst]
    print 'waiting...'
    for p in processes: p.wait()