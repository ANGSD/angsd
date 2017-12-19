#!/usr/bin/env python
import sys
import numpy as np
import scipy.optimize as optim
import logging
from copy import deepcopy
# from collections import namedtuple
# from copy import deepcopy
# import multiprocessing as mp

# NAMED_ROW = namedtuple("ROW", "relpos prime strand qual d t s counts")
"""
input:9 columns:
bam_idx, relpos, prime, strand, qual, base_out, base_ref, base_sampe, counts
bam_idx: index of the bam in the bamlist
relpos pos on read from the near prime, default (0-19) and 20 is the center
prime 0:5p; 1:3p; 2:c
strand: 0:+ and 1:-
qual: userdefined. defualt 0:0-19; 1:20-29; 2:30-100
A=0; C=1; G=2; T=3

output:
bam_idx, relpos, prime, strand, qual, base_d, base_s, errorrate
#paste <(sort -g anc_likes.errors.est ) <( sort -g anc_likes.errors.est2) | awk '{print $0, $8-$NF}' |  sort -gk17,17  |less -S
"""

logging.basicConfig(filename=sys.argv[1]+'.log', level=logging.INFO)
sys.stderr.write("please check {} for error matrices "
                 "replaced by global estimates\n".format(sys.argv[1]+'.log'))

# all indecies without diagonal
INDECES = [i for i in range(16) if i not in (0, 5, 10, 15)]
MAT_INDECES = ['_'.join([str(row), str(col)]) for row in range(4) for col in range(4)]
MINVAL = 10e-10
PSEUDOCOUNT = 0


def loglik2(par, pm, sm):
    emat = np.zeros(16)
    emat[INDECES] = par
    emat = emat.reshape((4, 4))
    np.fill_diagonal(emat, 1-emat.sum(axis=1))
    P = np.matmul(emat, pm)
    P[P < MINVAL] = MINVAL
    return -np.sum(np.log(P) * sm)


def ll_wrapper(out_per, out_sam):
    par2 = np.random.uniform(MINVAL, 0.01, size=12)
    conv = optim.minimize(loglik2, par2, method='L-BFGS-B',
                          bounds=[(MINVAL, 1-MINVAL) for x in range(12)],
                          args=(out_per, out_sam))

    while not conv.success:
        par2 = np.random.uniform(MINVAL, 0.01, size=12)
        conv = optim.minimize(loglik2, par2, method='L-BFGS-B',
                              bounds=[(MINVAL, 1-MINVAL) for x in range(12)],
                              args=(out_per, out_sam))
    results = np.zeros(16)
    results[INDECES] = conv.x
    results = results.reshape((4, 4))
    np.fill_diagonal(results, 1-results.sum(axis=1))
    return results


def write_mat_to_file(fh, group, mat):
    for i in range(4):
        for j in range(4):
            fh.write("{} {} {} {}\n".format(
                " ".join(map(str, group)), i, j, mat[i, j])
            )


global_counts_q = {}
all_data_per_group = {}
mykeys = set()
global_estimate_out_per_dic = {}
global_estimate_out_sam_dic = {}
global_estimate_out_per = np.zeros((4, 4), dtype=np.float)
global_estimate_out_sam = np.zeros((4, 4), dtype=np.float)
with open(sys.argv[1], 'r') as fh:
    data = fh.read()
    for line in data.split("\n"):
        if not line:
            continue
        bam_idx, relpos, prime, strand, qual, out, per, sam, counts = [int(x) for x in line.split()]
        group = (bam_idx, relpos, prime, strand, qual)
        mykeys.add(group)
        if group not in all_data_per_group.iterkeys():
            all_data_per_group[group] = []

        all_data_per_group[group].append((out, per, sam, counts))

        if bam_idx not in global_counts_q.iterkeys():
            global_counts_q[bam_idx] = {}
        if qual not in global_counts_q[bam_idx].iterkeys():
            global_counts_q[bam_idx][qual] = {}

        try:
            global_counts_q[bam_idx][qual][(out, per, sam)] += counts
        except KeyError:
            global_counts_q[bam_idx][qual][(out, per, sam)] = counts

        if bam_idx not in global_estimate_out_per_dic.iterkeys():
            global_estimate_out_per_dic[bam_idx] = np.zeros((4, 4), dtype=np.float)
            global_estimate_out_sam_dic[bam_idx] = np.zeros((4, 4), dtype=np.float)
        global_estimate_out_per_dic[bam_idx][out, per] += counts
        global_estimate_out_sam_dic[bam_idx][out, sam] += counts
sys.stderr.write("Estimate {} 4x4 matrices\n".format(len(mykeys)))
global_error_rates_perindi = {}
global_error_rates_perindi_qual = {}
for bam_idx in global_estimate_out_per_dic.keys():
    global_estimate_out_per = global_estimate_out_per_dic[bam_idx]
    global_estimate_out_per /= global_estimate_out_per.sum(axis=1)[:, np.newaxis]  # rowsum

    global_estimate_out_sam = global_estimate_out_sam_dic[bam_idx]

    global_estimate = ll_wrapper(global_estimate_out_per, global_estimate_out_sam)
    logging.info("Global estimates; Indi: {}; base observations {}".format(bam_idx, int(np.sum(global_estimate_out_sam))))
    logging.info(" ".join([coord+":"+str(round(x, 5)) for coord, x in zip(MAT_INDECES, global_estimate.flatten())]))
    global_error_rates_perindi[bam_idx] = deepcopy(global_estimate)

for bam_idx in global_counts_q.keys():
    for q, global_counts in global_counts_q[bam_idx].iteritems():
        out_per = np.zeros((4, 4), dtype=np.float) + PSEUDOCOUNT
        out_sam = np.zeros((4, 4), dtype=np.float) + PSEUDOCOUNT

        for (out, per, sam), counts in global_counts.iteritems():
            out_per[out, per] += counts
            out_sam[out, sam] += counts
        if np.any(np.diag(out_sam) < 5e3) or np.sum(out_sam) < 5e4:
            logging.warning("Indi: {}; Binqual {} does not contain sufficient"
                            " data to get proper estimates."
                            " Using global estimates "
                            "instead: {}. Insufficient Count mat:".format(bam_idx,q, q))
            logging.warning(" ".join([coord+":"+str(int(x))
                                      for coord, x in zip(MAT_INDECES, out_sam.flatten())]))
            global_error_rates_perindi_qual[(bam_idx, q)] = global_error_rates_perindi[bam_idx]
        else:
            out_per /= out_per.sum(axis=1)[:, np.newaxis]  # rowsum
            results = ll_wrapper(out_per, out_sam)
            global_error_rates_perindi_qual[(bam_idx, q)] = results
            logging.info("Indi: {}; Binqual: {}; base observations: {}".format(bam_idx, q, int(np.sum(out_sam))))
            logging.info(" ".join([coord+":"+str(round(x, 5)) for coord, x in zip(MAT_INDECES, results.flatten())]))

est_done=0

with open(sys.argv[1]+".est", 'w') as fhout:
    for group_key, group_data in all_data_per_group.iteritems():
        if est_done % 50 == 0:
            sys.stderr.write("{} done out of {}\r".format(est_done, len(mykeys)))
        out_per = np.zeros((4, 4), dtype=np.float) + PSEUDOCOUNT
        out_sam = np.zeros((4, 4), dtype=np.float) + PSEUDOCOUNT

        for (out, per, sam, counts) in group_data:
            out_per[out, per] += counts
            out_sam[out, sam] += counts

        if np.any(np.diag(out_sam) < 5000):
            logging.warning("{} does not contain sufficient"
                            " data to get proper estimates."
                            " Using estimates of quality bin {}"
                            " instead. Insufficient Count mat:".format(
                                " ".join(map(str, group_key)), group_key[-1]))
            logging.warning(" ".join([coord+":"+str(int(x))
                                      for coord, x in zip(MAT_INDECES, out_sam.flatten())]))

            # k[-1] == qualbin
            # k[0] == bam_idx
            write_mat_to_file(fhout, group_key,
                              global_error_rates_perindi_qual[(group_key[0],
                                                               group_key[-1])])
        else:
            out_per /= out_per.sum(axis=1)[:, np.newaxis]  # rowsum
            results = ll_wrapper(out_per, out_sam)
            write_mat_to_file(fhout, group_key, results)
        est_done += 1
sys.stderr.write("Estimated {} 4x4 matrices\n".format(est_done))
