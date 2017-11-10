#!/usr/bin/env python
import os, sys
import itertools
import numpy as np
import copy, time
import time
import cPickle as pickle

import RBNS_utils

##### Implements the Streaming Kmer Assignment (SKA) algorithm

def main():
    try:
        k, reads_file, num_passes, out_dir = sys.argv[1:]
        k = int(k)
        assert k in range(3,12)
        assert os.path.exists(reads_file)
        num_passes = int(num_passes)
        assert num_passes in range(2,20)
        assert os.direxists(out_dir)
    except:
        print 'Input argument error!'
        print ''
        print 'Usage:\npython streaming_convergence.py k',
        print 'reads_file num_passes output_directory'
        print 'k must be in [3,11]'
        print 'num_passes must be in [2, 19]'
        print 'readsfile must exist and only contain the read sequences (not fastq or fasta).'
        print 'out_dir must exist'
    SKA_weights = stream_counts(k, reads_file, num_passes)
    of = open(os.path.join(out_dir, 'SKA_weights.txt'), 'w')
    of.write('kmer\tweight\n')
    for kmer, weight in zip(yield_kmers(k), SKA_weights):
        of.write('%s\t%g\n' % (kmer, weight))
    of.close()



def stream_counts(
        k,
        inFile,
        num_iterations ):
    """
    - Main function to perform the SKA algorithm
    """
    current_weights = np.ones(4 ** k)
    for iteration_i in range(num_iterations):
        if iteration_i == 0:
            current_weights = stream_continual_update(k,
              current_weights, inFile)
        else:
            current_weights = stream_without_continual_update(k,
              current_weights, inFile)
        current_weights = copy.copy(current_weights)
    assert len(current_weights) == 4 ** k
    return current_weights



def make_table(
        k,
        weight_history,
        out_file ):
    """ Makes a table of kmer weight histories """
    of = open(out_file, 'w')
    of.write('kmer\t' + '\t'.join(
      ['round_%i' % i for i in range(len(weight_history))]) + '\n')
    for kmer_i, kmer in enumerate(yield_kmers(k)):
        of.write('%s\t' % kmer)

        for col_i in range(len(weight_history)):
            assert len(weight_history[col_i]) == 4 ** k
            of.write('%g\t' % weight_history[col_i][kmer_i])
        of.write('\n')
    of.close()




def make_sorted_table(
        k,
        weight_history,
        out_file ):
    """
    - Makes a table of kmer weight histories sorted by decreasing final weight
    """
    of = open(out_file, 'w')
    of.write('kmer\t' + '\t'.join(
      ['round_%i' % i for i in range(len(weight_history))]) + '\n')
    # first get the ordering of the kmers from highest to lowest
    kmers_ordered_L = []
    for kmer_i, kmer in enumerate(yield_kmers(k)):
        final_weight = weight_history[len(weight_history)-1][kmer_i]
        kmers_ordered_L.append( (kmer_i, kmer, final_weight) )
    kmers_ordered_L.sort(key = lambda tupl: -1*tupl[2])
    # now go through the ordered kmers and write to the out_file
    for tupl in kmers_ordered_L:
        kmer_i = tupl[0]
        kmer = tupl[1]
        of.write('%s\t' % kmer)

        for col_i in range(len(weight_history)):
            assert len(weight_history[col_i]) == 4 ** k
            of.write('%g\t' % weight_history[col_i][kmer_i])
        of.write('\n')
    of.close()





def stream_continual_update_with_convergence_table(
        k,
        weights,
        inFile,
        out_file,
        how_often_to_write = 10000):
    """
    - Performs streaming kmer assignment (SKA) algorithm in which the weights
        are continually updated after each read (typically, this is used for
        just the first pass); also makes an ouptput summary table at the end
    """
    internal_history = []
    for linei,line in enumerate( RBNS_utils.aopen( inFile ) ):
        if linei % how_often_to_write == 0:
            norm_weights = copy.copy( weights )
            norm_weights = normalize_mean_1( norm_weights )
            internal_history.append( norm_weights )

        read_seq = line.strip()
        pk = get_kmers( read_seq, k )
        assigned_weights = assign_kmer_weights( pk, weights )
        for kmer, weight in zip( pk, assigned_weights ):
            kmeri = get_index_from_kmer( kmer )
            weights[kmeri] += weight

    of = open(out_file, 'w')
    of.write('kmer\t' + '\t'.join(
      ['reads_read_%i' % (i * how_often_to_write) for i in\
              range(len(internal_history))]) + '\n')
    for kmer_i, kmer in enumerate( yield_kmers(k) ):

        of.write('%s\t' % kmer)
        for col_i in range(len(internal_history)):
            assert len(internal_history[col_i]) == 4 ** k
            of.write('%g\t' % internal_history[col_i][kmer_i])
        of.write('\n')

    of.close()

    return weights



def stream_continual_update(
        k,
        weights,
        inFile ):
    """
    - Performs streaming kmer assignment (SKA) algorithm in which the weights
        are continually updated after each read (typically, this is used for
        just the first pass)
    """
    total_lines = count_lines(inFile) * 2
    start_time = time.time()
    for linei, line in enumerate( RBNS_utils.aopen( inFile ) ):
        if linei % 10000 == 0 and linei:
            elapsed_time = time.time() - start_time
            print 'Predicted time remaining for stream_continual_update:',\
              (total_lines - linei) / linei * elapsed_time / 3600,\
              'hours'
        read_seq = line.strip()
        pk = get_kmers(read_seq, k)
        assigned_weights = assign_kmer_weights(pk, weights)
        for kmer, weight in zip(pk, assigned_weights):
            kmeri = get_index_from_kmer(kmer)
            weights[kmeri] += weight
    return weights


def stream_without_continual_update(
        k,
        in_weights,
        inFile ):
    """
    - Performs streaming kmer assignment (SKA) algorithm in which the weights
        are NOT continually updated after each read, just after going through
        all of the reads (typically used from the second pass onward)
    """
    new_weights = np.ones(4 ** k)
    for linei, line in enumerate( RBNS_utils.aopen( inFile ) ):
        read_seq = line.strip()
        pk = get_kmers(read_seq, k)
        additional_weights = assign_kmer_weights(pk, in_weights)
        assert sum(additional_weights) - 1.0 < 0.001
        for kmer, weight in zip(pk, additional_weights):
            kmeri = get_index_from_kmer(kmer)
            new_weights[kmeri] += weight
    return new_weights


#################################### < UTILS > ################################

def count_lines( inFile ):
    t = 0
    for l in open(inFile).xreadlines():
        t += 1
    return t

def assign_kmer_weights( pk, weights ):
    kmer_weight = np.array(map(
      lambda s: float(weights[get_index_from_kmer(s)]), pk))
    kmer_weight /= float(sum(kmer_weight))
    return kmer_weight


def get_kmers_no_homopolymers( seq, k ):
    """ Returns a set of all non-homopolymeric kmers within seq"""
    pk = []
    for i in range(0, len(seq) - k):
        kmer = seq[i:i + k]
        if len(set(kmer)) > 1: # no homo polymers
            pk.append(kmer)
    return pk


def get_kmers( seq, k ):
    """ Returns a set of all of the kmers within seq """
    pk = set()
    for i in range(0, len(seq) - k):
        kmer = seq[i:i + k]
        pk.add(kmer)
    return pk


def get_index_from_kmer( kmer ):
    """ Returns the base10 version of the base 4 DNA representation """
    index = 0
    base2face = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(kmer):
        if not base in 'ACGT':
            return -1
        power = len(kmer) - 1 - i
        index += base2face[base] * (4 ** power)
    return index


def yield_kmers( k ):
    """ An iterater to all kmers of length k in alphabetical order """
    bases = 'ACGT'
    for kmer in itertools.product(bases, repeat=k):
        yield ''.join(kmer)

def normalize_mean_1( ws ):
    """ Normalizes a list such that the average entry value is 1 """
    total = float(sum(ws))
    l = float(len(ws))
    ans = [w * l / total for w in ws]
    return ans


def normalize_sum_1( ws ):
    """ Normalizes a list to sum to 1 """
    total = float(sum(ws))
    ans = [w / total for w in ws]
    return ans

################################### </ UTILS > ################################


if __name__ == '__main__':
    main()




