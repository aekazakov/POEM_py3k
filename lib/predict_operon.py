#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
This is a copy of deep_operon.py, with only code for operon prediction in a genome
'''
import sys
import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils.CodonUsage import SynonymousCodons
import math
from collections import Counter
import os
# set backend
os.environ['KERAS_BACKEND'] = 'tensorflow'
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"

from keras import backend as K
try:
    from tensorflow import set_random_seed
except Exception:
    set_random_seed = np.random.seed
import keras
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint

# global parameters:
global_kmr = 3
global_len = 128
global_split_rate = 1./3
global_epoch = 32


def f1_score(y_true, y_pred):
    def recall(y_true, y_pred):
        """Recall metric.

        Only computes a batch-wise average of recall.

        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

    def precision(y_true, y_pred):
        """Precision metric.

        Only computes a batch-wise average of precision.

        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2 * ((precision * recall) / (precision + recall))


##########################################################################
# share kmer between 2 sequence
##########################################################################
def pearson(x, y):
    N, M = len(x), len(y)
    assert N == M
    x_m, y_m = sum(x) * 1. / N, sum(y) * 1. / M
    a, b, c = 0., 0., 0.
    for i in range(N):
        xi, yi = x[i] - x_m, y[i] - y_m
        a += xi * yi
        b += xi ** 2
        c += yi ** 2
    try:
        return a / math.sqrt(b * c)
    except Exception:
        return 0


def sharekmer(s1, s2):
    # SynonymousCodons
    n1, n2 = list(map(len, [s1, s2]))
    k1 = [s1[elem: elem + 3] for elem in range(0, n1, 3)]
    k2 = [s1[elem: elem + 3] for elem in range(0, n2, 3)]
    fq1 = Counter(k1)
    fq2 = Counter(k2)

    kmers = []
    for i in SynonymousCodons:
        j = SynonymousCodons[i]
        if len(j) < 2:
            continue
        c1 = [[fq1[elem], elem] for elem in j]
        best1 = max(c1, key=lambda x: x[0])
        c2 = [[fq2[elem], elem] for elem in j]
        best2 = max(c2, key=lambda x: x[0])

        if best1[1] == best2[1]:
            kmers.append(best1[1])

    vec1 = [fq1[elem] for elem in kmers]
    vec2 = [fq2[elem] for elem in kmers]

    return pearson(vec1, vec2)


##########################################################################
# the motif found
##########################################################################
box_up10 = ['TATAAT', [77, 76, 60, 61, 56, 82]]
box_up35 = ['TTGACA', [69, 79, 61, 56, 54, 54]]
# find the best region that may be a candidate of a motif


def find_motif(seq, motif, bg=None):
    if bg is None:
        bg = {}
    motif_len = len(motif[0])
    best = -100
    idx = -1
    for i in range(0, len(seq) - motif_len + 1):
        lmer = seq[i: i + motif_len]
        score = 0
        for a, b, c in zip(lmer, motif[0], motif[1]):
            if a == b:
                score += math.log(float(c) / bg.get(a, 1.))
            else:
                score += math.log((100. - c) / bg.get(a, 1.))

        if score >= best:
            idx = i
            best = score

    return [seq[idx: idx + motif_len], len(seq) - idx, best]


##########################################################################
# cai, from biopython
##########################################################################
index = Counter({'GCT': 1, 'CGT': 1, 'AAC': 1, 'GAC': 1, 'TGC': 1, 'CAG': 1,
                 'GAA': 1, 'GGT': 1, 'CAC': 1, 'ATC': 1, 'CTG': 1, 'AAA': 1,
                 'ATG': 1, 'TTC': 1, 'CCG': 1, 'TCT': 1, 'ACC': 1, 'TGG': 1,
                 'TAC': 1, 'GTT': 1, 'ACT': 0.965, 'TCC': 0.744, 'GGC': 0.724,
                 'GCA': 0.586, 'TGT': 0.5, 'GTA': 0.495, 'GAT': 0.434,
                 'GCG': 0.424, 'AGC': 0.41, 'CGC': 0.356, 'TTT': 0.296,
                 'CAT': 0.291, 'GAG': 0.259, 'AAG': 0.253, 'TAT': 0.239,
                 'GTG': 0.221, 'ATT': 0.185, 'CCA': 0.135, 'CAA': 0.124,
                 'GCC': 0.122, 'ACG': 0.099, 'AGT': 0.085, 'TCA': 0.077,
                 'ACA': 0.076, 'CCT': 0.07, 'GTC': 0.066, 'AAT': 0.051,
                 'CTT': 0.042, 'CTC': 0.037, 'TTA': 0.02, 'TTG': 0.02,
                 'GGG': 0.019, 'TCG': 0.017, 'CCC': 0.012, 'GGA': 0.01,
                 'CTA': 0.007, 'AGA': 0.004, 'CGA': 0.004, 'CGG': 0.004,
                 'ATA': 0.003, 'AGG': 0.002
                 })

def cai(seq):
    if seq.islower():
        seq = seq.upper()

    N = len(seq)
    cai_value, cai_length = 0, 0
    for i in range(0, N, 3):
        codon = seq[i: i + 3]
        if codon in index:
            if codon not in ['ATG', 'TGG']:
                cai_value += math.log(index[codon])
                cai_length += 1
        elif codon not in ['TGA', 'TAA', 'TAG']:
            continue
        else:
            continue

    if cai_length > 0:
        return math.exp(cai_value / cai_length)
    else:
        return 0


##########################################################################
# get the features
##########################################################################
# convert ATCG based kmer number
#code = {'A': 1, 'a': 1, 'T': 2, 't': 2, 'G': 3, 'g': 3, 'C': 4, 'c': 4}
code = [0] * 256
code5 = [0] * 256
flag = 0
for i in 'ATGC':
    code[ord(i.lower())] = code[ord(i)] = flag
    code5[ord(i.lower())] = code5[ord(i)] = flag + 1
    flag += 1


# convert string to number
def s2n(s, code=code, scale=None):
    if scale is None:
        scale = max(code) + 1
    N = 0
    output = 0
    for i in s[::-1]:
        output += code[ord(i)] * scale ** N
        N += 1

    return output


# reverse of s2n
def n2s(n, length, alpha='ATGC', scale=None):
    if scale is None:
        scale = max(code) + 1
    N = n
    s = []
    for i in range(length):
        s.append(alpha[N % scale])
        N /= scale

    return ''.join(s[::-1])


# convert the dna sequence to kmer-position matrix.
# if length of dna < given, then add NNN in the center of the sequence.
# else if length of dna > given, then trim the center of the sequence.

# the new kpm, reshape
def kpm(S, d=64, k=3, code=code, scale=None):
    if scale is None:
        scale = max(code) + 1

    N = scale ** k
    assert isinstance(d, int)
    L = len(S)
    if d < L:
        F = d // 2
        R = d - F
        seq = ''.join([S[: F], S[-R:]])
    elif d > L:
        F = L // 2
        R = L - F
        seq = ''.join([S[: F], 'N' * (d - L), S[-R:]])
    else:
        seq = S

    mat = [[0] * (d // k) for elem in range(N * k)]
    for i in range(0, d - k + 1):
        kmer = seq[i: i + k]
        if 'N' in kmer or 'n' in kmer:
            continue
        R = s2n(kmer, code=code, scale=scale)
        mat[R + i % k * N][i // k] = 1

    mat = np.asarray(mat, 'int8')
    return mat


# get features by give loc1, start and end:
# get xx
def get_xx(j, seq_dict, kmer=2, dim=128, mode='train', context=False):
    loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = j[: 10]
    if scf1 != scf2 or std1 != std2:
        if context:
            X0 = np.ones((4 ** kmer * kmer, dim // kmer * kmer))
        else:
            X0 = np.ones((4 ** kmer * kmer, dim // kmer))
        X1 = [10**4] * 11
        X2 = [127] * dim
        return [X0], X1, X2

    # get the sequence
    st1, ed1, st2, ed2 = list(map(int, [st1, ed1, st2, ed2]))
    st1 -= 1
    st2 -= 1

    if st1 > st2:
        loc1, scf1, std1, st1, ed1, loc2, scf2, std2, st2, ed2 = \
        loc2, scf2, std2, st2, ed2, loc1, scf1, std1, st1, ed1

    seq1 = seq_dict[scf1][st1: ed1]
    seq1 = std1 == '+' and seq1 or seq1.reverse_complement()
    seq2 = seq_dict[scf2][st2: ed2]
    seq2 = std1 == '+' and seq2 or seq2.reverse_complement()

    start, end = ed1, st2
    seq12 = seq_dict[scf1][start: end]

    seq12 = std1 == '+' and seq12 or seq12.reverse_complement()
    seq1, seq2, seq12 = list(map(str, [seq1.seq, seq2.seq, seq12.seq]))
    seq1, seq2, seq12 = seq1.upper(), seq2.upper(), seq12.upper()

    # 1D features such as gc, dist
    cai1, cai2, cai12 = list(map(cai, [seq1, seq2, seq12]))
    distn = (st2 - ed1) * 1. / (ed2 - st1)
    ratio = math.log((ed1 - st1) * 1. / (ed2 - st2))
    ratio = std1 == '+' and ratio or -ratio
    idx = -100
    bgs = Counter(seq12[idx:])
    up10, up35 = find_motif(seq12[idx:], box_up10, bgs), find_motif(
        seq12[idx:], box_up35, bgs)
    if seq12[idx:]:
        gc = SeqUtils.GC(seq12[idx:])
        try:
            skew = SeqUtils.GC_skew(seq12[idx:])[0]
        except Exception:
            skew = 0.
    else:
        gc = skew = 0.

    bias = sharekmer(seq1, seq2)
    if st1 == st2 == '+':
        X1 = [cai1, cai2, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]
    else:
        X1 = [cai2, cai1, bias, distn, ratio, gc, skew] + up10[1:] + up35[1:]

    # 2D features of kmer matrix
    if context:
        seqmat12 = kpm(seq12, d=dim, k=kmer, scale=4)
        seqmat1 = kpm(seq1, d=dim, k=kmer, scale=4)
        seqmat2 = kpm(seq2, d=dim, k=kmer, scale=4)
        seqmat = np.concatenate((seqmat1, seqmat12, seqmat2), 1)
    else:
        seqmat = kpm(seq12, d=dim, k=kmer, scale=4)

    if ed1 > st2:
        seqmat[:] = 0
    X0 = [seqmat]
    n12 = len(seq12)
    X2 = [s2n(seq12[elem: elem + kmer], code5)
          for elem in range(n12 - kmer + 1)]

    return X0, X1, X2


# get single line of features
def get_xx_one(j, seq_dict, kmer=2, dim=128, mode='train'):
    X0, X1, X2 = get_xx(j, seq_dict, kmer, dim, mode)
    x0, x1, x2 = list(map(np.asarray, [[X0], [X1], [X2]]))
    return x0, x1, X2

##########################################################################
# the CNN class
##########################################################################
class CNN:

    def __init__(self, nb_filter=32, nb_pool=3, nb_conv=2, nb_epoch=10,
                 batch_size=64, maxlen=128, save_path='./weights.hdf5'):
        # def __init__(self, nb_filter=64, nb_pool=3, nb_conv=2, nb_epoch=10,
        # batch_size=32, maxlen=128, save_path='./weights.hdf5'):

        self.nb_filter = nb_filter
        self.nb_pool = nb_pool
        self.nb_conv = nb_conv
        self.nb_epoch = nb_epoch
        self.batch_size = batch_size
        self.maxlen = maxlen
        self.opt = Adam(lr=5e-4, beta_1=0.995, beta_2=0.999, epsilon=1e-09)
        self.checkpointer = [ModelCheckpoint(
            filepath=save_path, verbose=1,
            save_best_only=True, mode='max', monitor='val_f1_score'
        )]
        self.metric = f1_score
        self.cross_val = 1. / 3

    def predict_2d(self, X):
        return self.model_2d.predict(X).argmax(1)

    # load an training model
    def load(self, name, mode='2d'):
        model = keras.models.load_model(
            name,
            custom_objects={'f1_score': f1_score, 'fbeta_score': f1_score},
            compile=False
        )
        if mode == '2d' or mode == 'lstm':
            self.model_2d = model
        else:
            pass


# run the whole genome prediction

# generate adjacent gene pairs from the gene list
def adjacent_genes(f):
    locus_list = []
    for i in f:
        j = i[: -1].split('\t')
        if len(j) < 7:
            j.extend([0] * 7)
        locus, scaf, strand, start, end = j[: 5]
        start, end = list(map(int, [start, end]))
        locus_list.append([locus, scaf, strand, start, end])

    locus_list.sort(key=lambda x: x[1: 5])
    return locus_list


def run_genome_predict(genome, seq_dict, model, clf, mode='2d'):
    #genome, model = sys.argv[3: 5]
    #seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    clf.load(model, mode)

    # get the locus of genes
    f = open(genome, 'r')
    locus_list = adjacent_genes(f)
    f.close()
    for a, b in zip(locus_list[: -1], locus_list[1:]):
        j = a + b
        #x0, x1, x2 = get_xx_one(j, seq_dict, 3, 128, 'test')
        x0, x1, x2 = get_xx_one(j, seq_dict, global_kmr, global_len, 'test')
        if mode == '2d':
            if a[1] == b[1] and a[2] == b[2]:
                res = clf.predict_2d(x0)[0]
            else:
                res = 0

        res = res == 1 and 'True' or 'False'
        i = '\t'.join(map(str, j))
        print(i + '\t' + str(res))


def print_usage():
    print('#' * 79)
    print('# To make a adjacent genes prediction')
    print('#' * 79)
    print('python this.py predict foo.fasta foo.adjacent.txt foo.model [mode]\n')
    print('foo.adjacent.txt is the gene location in the format:')
    print('       locus1\tscf1\tstrand1\tstart1\tend1\t'
          'locus2\tscf2\tstrand2\tstart2\tend2\n')


if __name__ == '__main__':

    if len(sys.argv[1:]) < 4:
        print_usage()
        sys.exit(0)
    
    command, fasta = sys.argv[1: 3]
    # save the genome to an dict
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

    if command.startswith('predict'):
        test, model = sys.argv[3: 5]
        mode = '2d'
        clf = CNN(nb_epoch=32, maxlen=128)
        run_genome_predict(test, seq_dict, model, clf, mode)
    else:
        print_usage()
        sys.exit(0)
