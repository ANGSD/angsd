#!/usr/bin/env python
import sys
import numpy as np
import scipy.optimize as optim
from collections import namedtuple
from copy import deepcopy
import multiprocessing as mp

NAMED_ROW = namedtuple("ROW", "relpos prime strand qual d t s counts")
MINVAL = 10e-10
PSEUDOCOUNT = 1

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
"""
class GeneralData(object):
    def __init__(self):
        self.groups = set()
        self.all_data = list()
        self.donkey_perfect = {}
        self.donkey_sample = {}
        self.best_params = {}
        self._mat_cache = {}
        self.bin_combinations = []

    def _get_bases(self):
        for outgroup in range(4):
            for per_or_sample in range(4):
                yield outgroup, per_or_sample

    def update(self, row):
        self.groups.add((row.relpos, row.prime, row.strand, row.qual))
        self.all_data.append(row)

    def load_data_to_dic(self):
        temp_dic = {}
        for combi in self.groups:
            temp_dic[combi] = {}
            for outgroup, per_or_sample in self._get_bases():
                temp_dic[combi][(outgroup, per_or_sample)] = PSEUDOCOUNT

        self.donkey_perfect = deepcopy(temp_dic)
        self.donkey_sample = deepcopy(temp_dic)

        for r in self.all_data:
            self.donkey_perfect[(r.relpos, r.prime, r.strand, r.qual)][(r.d, r.t)] += r.counts
            self.donkey_sample[(r.relpos, r.prime, r.strand, r.qual)][(r.d, r.s)] += r.counts

    def yield_general(self):
        for mygroup in self.groups:
            mat_d_t, mat_d_s = self._get_matrix(mygroup)
            if mat_d_t is None:
                continue
            yield mygroup, mat_d_t, mat_d_s

    def _get_matrix(self, mygroup):
        if mygroup in self._mat_cache.iterkeys():
            return self._mat_cache[mygroup]

        perfectmatrix = list()
        samplematrix = list()
        for outgroup, per_or_sample in self._get_bases():
            perfectmatrix.append(
                self.donkey_perfect[mygroup][(outgroup, per_or_sample)]
            )
            samplematrix.append(
                self.donkey_sample[mygroup][(outgroup, per_or_sample)]
            )

        assert len(perfectmatrix) == 16, "FUCKED"
        assert len(samplematrix) == 16, "FUCKED"

        samplematrix = np.asarray(samplematrix).reshape((4, 4))
        perfectmatrix = np.asarray(perfectmatrix, dtype=np.float).reshape((4, 4))
        ## this is with pseudocount is included
        if np.all(perfectmatrix == 1) or np.all(samplematrix == 1):
            # sys.stderr.write("Excluding: {} {} {}\n".format(relpos, prime, strand, qual))
            self._mat_cache[mygroup] = (None, None)
            return None, None
        if np.all(perfectmatrix == 0) or np.all(samplematrix == 0):
            # sys.stderr.write("Excluding: {} {} {}\n".format(relpos, prime, strand, qual))
            self._mat_cache[mygroup] = (None, None)
            return None, None
        if np.sum(perfectmatrix) == 16 or np.sum(samplematrix) == 16:
            # sys.stderr.write("Excluding: {} {} {}\n".format(relpos, prime, strand, qual))
            self._mat_cache[mygroup] = (None, None)
            return None, None
        # https://stackoverflow.com/a/19602209
        perfectmatrix /= perfectmatrix.sum(axis=1)[:, np.newaxis]  # rowsum
        self._mat_cache[mygroup] = (perfectmatrix, samplematrix)
        return perfectmatrix, samplematrix


def loglik(par, pm, sm):
    emat = par.reshape(4, 4)
    np.fill_diagonal(emat, 0)
    np.fill_diagonal(emat, 1-emat.sum(axis=1))
    P = np.matmul(emat, pm)
    if np.any(P < MINVAL):
        P[P < MINVAL] = MINVAL
    return -np.sum(np.log(P) * sm)


def worker_helper(args):
    return worker(*args)


def worker(mat_d_t, mat_d_s, my_id):
    par = np.random.uniform(MINVAL, 0.01, size=16).reshape((4, 4))
    weird_rows = []
    for idx, val in enumerate(np.diag(mat_d_s)):
        # if mat_d_s[idx,idx] < (np.sum(mat_d_s[idx,:])-mat_d_s[idx,idx]):
        if (mat_d_s[idx,idx]) < (2*(np.sum(mat_d_s[idx,:])-mat_d_s[idx,idx])):
            mat_d_s[idx, :] = MINVAL
            mat_d_t[idx, :] = MINVAL
            weird_rows.append(idx)
        
    conv = optim.minimize(loglik, par, method='L-BFGS-B',
                          bounds=[(MINVAL, 1-MINVAL) for x in range(16)],
                          args=(mat_d_t, mat_d_s))


    # print(my_id)
    # print(mat_d_t)
    # print(mat_d_s)
    # print np.round(best_conv.x.reshape(4, 4), 5)
    while not conv.success:
        par = np.random.uniform(MINVAL, 0.1, size=16).reshape((4, 4))
        conv = optim.minimize(loglik, par, method='L-BFGS-B',
                            bounds=[(MINVAL, 1-MINVAL) for x in range(16)],
                              args=(mat_d_t, mat_d_s))

    optim_par = conv.x.reshape((4, 4))

    for idx in weird_rows:
        optim_par[idx, : ] = MINVAL

    return (my_id, optim_par)


bams = {}
if len(sys.argv) == 3:
    mypool = mp.Pool(processes=int(sys.argv[2]))
else:
    mypool = mp.Pool(processes=1)

with open(sys.argv[1], 'r') as fh:
    for line in fh.read().split("\n"):
        if not line:
            continue
        bam_idx, relpos, prime, strand, qual, d, t, s, counts = [int(x) for x in line.split()]

        if(counts < 0):
            sys.stderr.write("Counts are negative. Not Good. Exiting\n")
            exit(1)

        if bam_idx not in bams.keys():
            bams[bam_idx] = GeneralData()
        bams[bam_idx].update(NAMED_ROW(relpos, prime, strand, qual, d, t, s, counts))

todos = []
for bam_idx, my_data in sorted(bams.iteritems()):
    my_data.load_data_to_dic()
    for (relpos, prime, strand, qual), mat_d_t, mat_d_s in my_data.yield_general():
        todos.append((mat_d_t, mat_d_s,
                      (bam_idx, relpos, prime, strand, qual)))

for mygroup, optim_par in mypool.imap_unordered(worker_helper, todos, chunksize=100):
    bams[mygroup[0]].best_params[mygroup] = optim_par

with open(sys.argv[1]+".est", 'w') as fhout:
    for bam_idx, my_data in sorted(bams.iteritems()):
        for (relpos, prime, strand, qual), mat_d_t, mat_d_s in my_data.yield_general():
            optim_par = my_data.best_params[(bam_idx, relpos, prime, strand, qual)]
            for base_d in range(4):
                for base_s in range(4):
                    fhout.write("{} {} {} {} {} {} {} {}\n".format(bam_idx, relpos, prime, strand, qual, base_d, base_s, optim_par[base_d, base_s]))
