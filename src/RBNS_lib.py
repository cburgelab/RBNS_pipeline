#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
import subprocess
import os, time
import cPickle
import RBNS_utils
import numpy as np
import math

import RBNS_cluster_utils

class RBNS_Lib:
    def __init__( self,
            experiment_settings,
            lib_settings ):
        """
        An RBNS library's settings & functions to return information
            about that library
        """
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.load_counts()
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir

    def get_conc(self):
        return self.lib_settings.conc

    def return_conc_annot_like_20_nM(self):
        return self.lib_settings.return_conc_annot_like_320_nM()

    def is_input(self):
        return self.lib_settings.is_input()

    def is_0nM(self):
        return self.lib_settings.is_0nM()

    def is_denom_for_enrichments(self):
        """
        - Returns True if this library is used as the denominator for the
            main enrichments (either the input or 0 nM library)
        """
        if (self.get_property( "input_library_for_logos" ) == "input"):
            if self.is_input():
                return True
            else:
                return False
        elif (self.get_property("input_library_for_logos") == "0") or\
                (self.get_property("input_library_for_logos") == "0 nM") or\
                (self.get_property("input_library_for_logos") == "0 nm"):
            if self.is_0nM():
                return True
            else:
                return False
        else:
            print "ERROR: NO INPUT LIBRARY FOR LOGOS"

    def get_barcode(self):
        return self.lib_settings.barcode

    def get_conc_for_fastq(self):
        return self.lib_settings.conc_for_fastq

    def get_conc_for_title(self):
        conc_for_fastq = self.get_conc_for_fastq()
        if (conc_for_fastq.find( 'input' ) == -1):
            conc_for_fastq += ' nM'
        return conc_for_fastq

    def get_stream_libfrac(self, k, kmer_i):
        return self.type2k2counts['stream'][k].get_libfrac(kmer_i)

    def get_stream_libfrac_kmer(self, kmer):
        k = len(kmer)
        return self.type2k2counts['stream'][k].kmer_value(kmer)

    def load_counts(self):
        self.type2k2counts = {}
        for count_type in ['naive', 'naive_max_once', 'stream']:
            if not self.experiment_settings.isCountRequested(count_type):
                continue
            self.type2k2counts[count_type] = {}
            for k in self.experiment_settings.get_ks(count_type):
                self.type2k2counts[count_type][k] =\
                  RBNS_profile(self.lib_settings, k, count_type)

    def get_naive_counts(self, k):
        return self.type2k2counts['naive'][k].get_profile()

    def get_naive_libfracs(self, k):
        return self.type2k2counts['naive'][k].get_libfracs()

    def get_stream_libfracs(self,k):
        return self.type2k2counts['stream'][k].get_libfracs()

    def get_naive_count_by_index(self, k, kmer_i):
        return self.type2k2counts['naive'][k].kmeri_value(kmer_i)

    def get_naive_count_kmer(self, kmer):
        k = len(kmer)
        return self.type2k2counts['naive'][k].kmer_value(kmer)

    def get_enrichments(self, k):
        return self.type2k2counts['naive'][k].get_enrichments()

    def get_enrichments_max_once(self, k):
        return self.type2k2counts['naive_max_once'][k].get_enrichments()

    def get_enrichment(self, k, kmer_i):
        return self.type2k2counts['naive'][k].get_enrichment(kmer_i)

    def get_enrichment_kmer(self, kmer):
        k = len(kmer)
        kmer_i = RBNS_utils.get_index_from_kmer(kmer)
        return self.get_enrichment(k, kmer_i)

    def calculate_enrichment(self, k, input_lib):
        enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'enrichment_Ds',
          '%s_%s_to_input.%imer.enrichments.pkl'  %
          (self.experiment_settings.get_property('protein_name'),
          self.lib_settings.get_conc_string(), k))
        if RBNS_utils.file_exists(enrich_pkl):
            self.type2k2counts['naive'][k].load_enrichments( enrich_pkl )
        else:
            input_profile = input_lib.type2k2counts['naive'][k]
            self.type2k2counts['naive'][k].calculate_enrichments( input_profile )
            self.type2k2counts['naive'][k].save_enrichments( enrich_pkl )


    def calculate_enrichment_max_once(self, k, input_lib):
        enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'enrichment_Ds',
          '%s_%s_to_input.%imer.enrichments.max_once.pkl'  %
          (self.experiment_settings.get_property('protein_name'),
          self.lib_settings.get_conc_string(), k))
        if RBNS_utils.file_exists(enrich_pkl):
            self.type2k2counts['naive_max_once'][k].load_enrichments( enrich_pkl )
        else:
            input_profile = input_lib.type2k2counts['naive_max_once'][k]
            self.type2k2counts['naive_max_once'][k].calculate_enrichments( input_profile )
            self.type2k2counts['naive_max_once'][k].save_enrichments( enrich_pkl )


    def return_enrich_D_pkl_to_input(self, k):
        enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'enrichment_Ds',
          '%s_%s_to_input.%imer.enrichments.pkl'  %
          (self.experiment_settings.get_property('protein_name'),
          self.lib_settings.get_conc_string(), k))
        return enrich_pkl

    def return_frequencies_D(self, k):
        for i in range( 10 ):
            freqs_D_F = os.path.join( self.experiment_settings.get_rdir(),
                    "frequency_Ds/{0}_{1}.{2}mer.frequencies.pkl".format(
                        self.experiment_settings.get_property('protein_name'),
                        self.get_conc_for_fastq(),
                        k ) )
            try:
                freqs_D = cPickle.load( open( freqs_D_F ) )
                return freqs_D
            except EOFError:
                time.sleep( 2 )

    def get_0nM_enrichments(self, k):
        return self.type2k2counts['naive'][k].get_0nM_enrichments()

    def get_0nM_enrichments_max_once(self, k):
        return self.type2k2counts['naive_max_once'][k].get_0nM_enrichments()

    def get_0nM_enrichment(self, k, kmer_i):
        return self.type2k2counts['naive'][k].get_0nM_enrichment(kmer_i)

    def get_0nM_enrichment_kmer(self, kmer):
        k = len(kmer)
        kmer_i = RBNS_utils.get_index_from_kmer(kmer)
        return self.get_0nM_enrichment(k, kmer_i)

    def calculate_enrichment_to_0nM(self, k, zero_nM_lib):
        zero_enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'enrichment_Ds',
          '%s_%s_to_0nM.%imer.enrichments.pkl'  %
          (self.experiment_settings.get_property('protein_name'),
          self.lib_settings.get_conc_string(), k))

        zero_nM_profile = zero_nM_lib.type2k2counts['naive'][k]
        self.type2k2counts['naive'][k].calculate_enrichments_to_0nM(zero_nM_profile)

        if RBNS_utils.file_exists(zero_enrich_pkl) == False:
            self.type2k2counts['naive'][k].save_0nM_enrichments(zero_enrich_pkl)

    def return_enrich_D_pkl_to_0nM(self, k):
        enrich_pkl = os.path.join(
          self.experiment_settings.get_rdir(),
          'enrichment_Ds',
          '%s_%s_to_0nM.%imer.enrichments.pkl'  %
          (self.experiment_settings.get_property('protein_name'),
          self.lib_settings.get_conc_string(), k))
        return enrich_pkl

    def return_main_enrichments_D(self, k):
        """
        - Assuming this is the not the "input_library_for_logos", returns
            a dictionary of the enrichments (to the input or 0nM library, as
            appropriate)
        """
        if (self.settings.get_property( "input_library_for_logos" ) == "input"):
            return self.return_enrich_D_pkl_to_input( k )
        elif (self.settings.get_property( "input_library_for_logos" ) == "0"):
            return self.return_enrich_D_pkl_to_0nM( k )

    def get_max_enrichment(self, k_for_max_enrichment):
        """
        returns the largest enrichment kmer of size k
        """
        return max(self.get_enrichments(k_for_max_enrichment))

    def get_max_0nM_enrichment(self, k_for_max_enrichment):
        """
        returns the largest enrichment kmer of size k
        """
        return max(self.get_0nM_enrichments(k_for_max_enrichment))

    def split_reads_exist(self):
        """
        returns true if the split reads file for this library exists
        and is non empty
        does not check if it is complete
        """
        return RBNS_utils.file_exists(self.get_split_reads())

    def get_split_reads(self):
        """
        returns the full path of this library's split reads file
        """
        split_reads_file = os.path.join(
          self.get_rdir(),
          'split_reads',
          '%s_%s.reads' % (
          self.experiment_settings.get_property('protein_name'),
          self.get_conc_for_fastq()))
        assert split_reads_file == self.lib_settings.get_split_reads()
        return self.lib_settings.get_split_reads()

    def get_split_fastqs(self):
        """
        returns the full path of this library's split reads file
        """
        split_fastq_file = os.path.join(
          self.get_rdir(),
          'split_reads',
          '%s_%s.fastq.gz' % (
          self.experiment_settings.get_property('protein_name'),
          self.get_conc_for_fastq()))
        assert split_fastq_file == self.lib_settings.get_split_fastqs()
        return self.lib_settings.get_split_fastqs()

    def get_split_whandle(self):
        """
        returns a write file handle to the split reads
        """
        return RBNS_utils.aopen(self.get_split_reads(), 'w')

    def get_naive_enrichment_dict(self, k):
        """
        returns a dictionary of kmer -> naive enrichment
        """
        return self.type2k2counts['naive'][k].get_enrichment_dict()

    def get_naive_0nM_enrichment_dict(self, k):
        """
        returns a dictionary of kmer -> naive enrichment
        """
        return self.type2k2counts['naive'][k].get_0nM_enrichment_dict()

    def get_B_kmer(self, kmer, read_len):
        return self.type2k2counts['naive'][len(kmer)].get_B_kmer(kmer, read_len)

    def calcB(self, kmer):
        """
         calculates the B value
        """
        read_len = self.experiment_settings.get_property('read_len')
        return self.type2k2counts['naive'][len(kmer)].get_B_kmer(kmer, read_len)






class RBNS_profile:
    """ The 'profile' (enrichments & counts) of an RBNS library """
    def __init__(self, lib_settings, k, count_type):
        self.count_type = count_type
        self.k = k
        counts_pkl = lib_settings.counts_file(count_type, k)
        self.profile = cPickle.load(open(counts_pkl, 'rb'))
        assert len(self.profile) == 4 ** self.k
        if (self.count_type == 'naive') or\
                (self.count_type == 'naive_max_once') or\
                (self.count_type == 'stream'):
            self.calculate_libfracs()
            assert 0.999 < sum(self.libfrac) < 1.001

    def kmer_value(self, kmer):
        kmeri = RBNS_utils.get_index_from_kmer(kmer)
        return self.kmeri_value(kmeri)

    def kmeri_value(self, kmeri):
        return self.profile[kmeri]

    def calculate_libfracs(self):
        counts = np.array(self.profile, dtype=float)
        total_counts = sum(counts)
        if not total_counts:
            raise ValueError('Naive counts are 0')
        # self.libfrac is a Numpy array of length 4^k
        self.libfrac = counts / total_counts

    def get_libfracs(self):
        assert self.count_type == 'naive' or self.count_type == 'stream' or\
                self.count_type == 'naive_max_once'
        return self.libfrac

    def get_libfrac(self, kmer_i):
        assert self.count_type == 'naive' or self.count_type == 'stream' or\
                self.count_type == 'naive_max_once'
        return self.libfrac[kmer_i]

    def get_B(self, kmer_i, read_len):
        assert kmer_i
        enrichment = self.get_enrichment(kmer_i)
        B = RBNS_utils.B_factor(enrichment, self.k, read_len)
        return B

    def get_B_kmer(self, kmer, read_len):
        assert len(kmer) == self.k
        enrichment = self.get_enrichment_kmer(kmer)
        B = RBNS_utils.B_factor(enrichment, self.k, read_len)
        return B

    def get_B_values(self, read_len):
        return [RBNS_utils.B_factor(enrich, self.k, read_len) for enrich in self.enrichments]

    def weight_dict(self):
        kmer2weight = {}
        for kmer, weight in zip(RBNS_utils.yield_kmers(self.k), self.profile):
            kmer2weight[kmer] = weight
        return kmer2weight

    def get_profile(self):
        return self.profile

    def save_enrichments(self, enrich_pkl):
        RBNS_utils.make_dir( os.path.dirname( enrich_pkl ) )
        enriches_by_kmer_D = {}
        num_kmers = len( self.enrichments )
        k = int( math.log( num_kmers, 4. ) )
        for kmer_num, kmer in enumerate( RBNS_utils.yield_kmers( k )):
            enriches_by_kmer_D[kmer] = self.enrichments[kmer_num]
        cPickle.dump(enriches_by_kmer_D, open(enrich_pkl, 'wb'))

    def load_enrichments(self, enrich_pkl):
        enriches_by_kmer_D = cPickle.load(open(enrich_pkl, 'rb'))
        k = int( math.log( len( enriches_by_kmer_D ), 4. ) )
        enriches_L = []
        for kmer in RBNS_utils.yield_kmers( k ):
            enriches_L.append( enriches_by_kmer_D[kmer] )
        self.enrichments = np.array( enriches_L )

    def calculate_enrichments(self, input_profile):
        self.enrichments = norm_libfracs(self, input_profile)

    def save_0nM_enrichments(self, enrich_pkl):
        enriches_by_kmer_D = {}
        num_kmers = len( self.enrichments_to_0nM )
        k = int( math.log( num_kmers, 4. ) )
        for kmer_num, kmer in enumerate( RBNS_utils.yield_kmers( k )):
            enriches_by_kmer_D[kmer] = self.enrichments_to_0nM[kmer_num]
        cPickle.dump( enriches_by_kmer_D, open(enrich_pkl, 'wb') )

    def load_0nM_enrichments(self, enrich_pkl):
        enriches_by_kmer_D = cPickle.load(open(enrich_pkl, 'rb'))
        k = int( math.log( len( enriches_by_kmer_D ), 4. ) )
        enriches_L = []
        for kmer in RBNS_utils.yield_kmers( k ):
            enriches_L.append( enriches_by_kmer_D[kmer] )
        self.enrichments_to_0nM = np.array( enriches_L )

    def calculate_enrichments_to_0nM(self, zero_nM_profile):
        self.enrichments_to_0nM = norm_libfracs(self, zero_nM_profile)

    def get_enrichment_kmer(self, kmer):
        kmer_i = RBNS_utils.get_index_from_kmer(kmer)
        return self.get_enrichment(kmer_i)

    def get_enrichments(self):
        return self.enrichments

    def get_enrichment(self, kmer_i):
        return self.enrichments[kmer_i]

    def get_enrichment_dict(self):
        assert len(self.enrichments) == 4 ** self.k
        return {kmer: enrich
                for kmer, enrich in
                zip(RBNS_utils.yield_kmers(self.k),
                self.enrichments)}


    def get_0nM_enrichment_kmer(self, kmer):
        kmer_i = RBNS_utils.get_index_from_kmer(kmer)
        return self.get_0nM_enrichment(kmer_i)

    def get_0nM_enrichments(self):
        return self.enrichments_to_0nM

    def get_0nM_enrichment(self, kmer_i):
        return self.enrichments_to_0nM[kmer_i]


    def get_0nM_enrichment_dict(self):
        assert len(self.enrichments_to_0nM) == 4 ** self.k
        return {kmer: enrich
                for kmer, enrich in
                zip(RBNS_utils.yield_kmers(self.k),
                self.enrichments_to_0nM )}




############################ < HELPER FUNCTIONS > #############################

def norm_libfracs( profile, input_profile ):
    assert profile.k == input_profile.k
    return profile.get_libfracs() / input_profile.get_libfracs()

########################### </ HELPER FUNCTIONS > #############################






