#!/usr/bin/env python
import os
import ConfigParser
try:
    import simplejson
except ImportError:
    import json
import itertools
import pprint

import RBNS_utils

class RBNS_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)

    def get_fastq(self):
        return os.path.abspath(self.settings['fastq'])

    def set_fastq(self, fastq_file):
        self.settings['fastq'] = fastq_file

    def get_ks(self, count_type):
        assert 'ks_to_test_%s' % count_type in self.settings
        settings_ks = self.settings['ks_to_test_%s' % count_type]
        return settings_ks

    def get_naiveks(self):
        return self.settings['ks_to_test_naive']

    def get_naiveks_max_once(self):
        return self.settings['ks_to_test_naive_max_once']

    def get_force_recount(self, count_type):
        return self.settings['force_%s_recount' % count_type]

    def get_settings_file(self):
        return self.settings_file

    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)

    def get_rdir(self):
        RBNS_utils.make_dir( self.rdir )
        return self.rdir

    def get_edir(self):
        RBNS_utils.make_dir(self.edir)
        return self.edir

    def get_input_barcode(self):
        return self.settings['input_barcode']

    def get_0nM_barcode(self):
        return self.settings['0_nm_barcode']

    def isCountRequested(self, count_type):
        if count_type == 'naive':
            return self.settings['naive_count']
        elif count_type == 'stream':
            return self.settings['stream_count']
        elif count_type == 'by_position':
            return self.settings['by_position_count']

    def get_conc_for_mostenriched_analyses(self):
        return self.settings['conc_for_mostenriched_analyses']

    def iter_lib_settings(self):
        for i in range(len(self.barcodes)):
            yield RBNS_lib_settings(self,
              self.barcodes[i],
              self.settings['concentrations'][i])

    def process_settings(self, settings_file):
        """
        - reads the settings file and converts str to float, list, etc.
        - stores result in self.settings as a dict()
        """
        #### Defaults for arguments not specified in the settings file
        defaults_D = {
                'begin_barcode_symb': '#',
                'by_position_count': False,
                'end_barcode_symb': '/',
                'fold_all_or_mostenrichedconc_only': 'most_enriched_only',
                'force_by_position_recount': False,
                'force_naive_max_once_recount': False,
                'force_naive_recount': False,
                'force_redo_logos': False,
                'force_stream_recount': False,
                'input_library_for_logos': 'input',
                'kmers_to_ignore': [],
                'ks_to_test_by_position': [],
                'ks_to_test_naive_max_once': [],
                'lower_ks_to_test_logos': [4],
                'mismatches_allowed_in_barcode': 0,
                'naive_count': True,
                'naive_max_once_count': False,
                'nt_freqs_by_position': True,
                'num_reads_for_logos': 0.5,
                #### adapter_sequences_l should be the DNA, not RNA, template adapters
                'adapter_sequences_l': ["CCTTGACACCCGAGAATTCCA",
                    "GATCGTCGGACTGTAGAACTCCCTATAGTGAGTCGT"],
                #### rna_5p_adapter / rna_5p_adapter are used in RNA folding
                'rna_5p_adapter': "GGGGAGTTCTACAGTCCGACGATC",
                'rna_3p_adapter': "TGGAATTCTCGGGTGTCAAGG",
                'stream_count': False,
                'temp': 4,
                'weblogo_path': '',
                'z_scores_for_logos': [3.] }

        int_keys = [
                'read_len',
                'mismatches_allowed_in_barcode',
                'temp']

        int_or_float_keys = ['num_reads_for_logos']

        float_keys = []

        boolean_keys = [
                'by_position_count',
                'force_by_position_recount',
                'force_stream_recount',
                'force_naive_recount',
                'force_redo_logos',
                'naive_count',
                'naive_max_once_count',
                'nt_freqs_by_position',
                'stream_count']

        list_str_keys = [
                'barcodes',
                'kmers_to_ignore']

        list_int_keys =[
                'ks_to_test_by_position',
                'ks_to_test_naive',
                'ks_to_test_naive_max_once',
                'ks_to_test_stream',
                'ks_to_test_logos',
                'lower_ks_to_test_logos']

        list_float_keys = ['concentrations',
                "z_scores_for_logos"]

        list_of_tuples_keys = []

        extant_files = ['fastq']

        config = ConfigParser.ConfigParser()
        config.read(settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True

        for key in int_keys:
            try:
                settings[key] = int(settings[key])
            except KeyError:
                settings[key] = int( defaults_D[key] )

        for key in float_keys:
            settings[key] = float(settings[key])

        for key in int_or_float_keys:
            try:
                int_or_float_key = settings[key]
                if (int( int_or_float_key ) == int_or_float_key):
                    settings[key] = int( int_or_float_key )
                else:
                    settings[key] = float( int_or_float_key )
            except KeyError:
                settings[key] = defaults_D[key]

        for key in boolean_keys:
            try:
                settings_bool = settings[key]
                if (settings_bool == 'true') or (settings_bool == 'True') or\
                        (settings_bool.strip('"') == 'true') or\
                        (settings_bool.strip('"') == 'True') or\
                        (settings_bool.strip("'") == 'true') or\
                        (settings_bool.strip("'") == 'True'):
                    settings[key] = True
                else:
                    settings[key] = False
            except KeyError:
                settings[key] = defaults_D[key]

        for key in list_float_keys:
            try:
                settings[key] = map(float, simplejson.loads(settings[key]))
            except KeyError:
                settings[key] = defaults_D[key]

        for key in list_int_keys:
            try:
                settings[key] = map(int, simplejson.loads(settings[key]))
            except KeyError:
                settings[key] = defaults_D[key]

        for key in list_str_keys:
            try:
                settings[key] = simplejson.loads(settings[key])
            except KeyError:
                settings[key] = defaults_D[key]

        for key in list_of_tuples_keys:
            try:
                settings[key] = simplejson.loads(settings[key])
            except KeyError:
                settings[key] = defaults_D[key]

        #### Check to make sure the .fastq exists
        try:
            assert RBNS_utils.file_exists(settings['fastq'])
        except AssertionError:
            curr_DIR = os.path.dirname( os.path.realpath(__file__) )
            fastq_F = os.path.join( curr_DIR, 'test_data', os.path.basename(settings['fastq']) )
            settings['fastq'] = fastq_F
            try:
                assert( RBNS_utils.file_exists(settings['fastq']) )
            except:
                print "\nERROR: FASTQ ({}) does not exist".format(
                            settings['fastq'] )
                print "\t-> Change fastq path in {}\n\n".format( self.settings_file )
                exit(1)

        #### See if an adapter sequence was passed in to look for
        ####    randomer reads of the correct length
        if ("start_of_adapter_seq" not in settings):
            settings["start_of_adapter_seq"] = None
        else:
            settings["start_of_adapter_seq"] =\
                    settings["start_of_adapter_seq"].strip('"')

        #### See if a protein_name_for_plotting was included; if not, make it
        ####    the protein_name
        if ("protein_name_for_plotting" not in settings):
            settings["protein_name_for_plotting"] =\
                    settings["protein_name"]
        else:
            settings["protein_name_for_plotting"] =\
                    settings["protein_name_for_plotting"].strip( "'" ).\
                    strip('"')

        #### Go through all of the defaults_D keys, and if they aren't yet in
        ####    settings, add them
        for key, default_value in defaults_D.iteritems():
            if key not in settings:
                settings[key] = default_value

        #### Make sure that the input_library_for_logos is properly
        ####    escaped if it is the 0 nM library (i.e., it should be entered
        ####    in the settings file as input_library_for_logos = 0 nM,
        ####    NOT input_library_for_logos = "0 nM"
        assert( settings['input_library_for_logos'] in ['input', '0 nM', '0 nm'] )

        conc_str_for_Fnames_L = []
        for conc, barcode in zip( settings['concentrations'], settings['barcodes'] ):
            if ( conc == 0 ):
                if ( settings['input_barcode'] == barcode ):
                    conc_str_for_Fnames_L.append( 'input' )
                else:
                    conc_str_for_Fnames_L.append( '0' )
            else:
                if float( int( conc ) ) == float( conc ):
                    conc_str_for_Fnames_L.append( str( int( conc ) ) )
                else:
                    conc_str_for_Fnames_L.append( "{:.1f}".format( conc ) )
        settings['conc_str_for_Fnames_L'] = conc_str_for_Fnames_L

        self.settings = settings
        #### Make sure that the results dir. can be made; if not, make
        ####    it one directory up from the pipeline directory
        self.rdir = settings['results_dir']
        made_DIR = RBNS_utils.make_dir( self.rdir )
        print "made_DIR: {}".format( made_DIR )
        if not made_DIR:
            curr_DIR = os.path.dirname( os.path.realpath(__file__) )
            parent_DIR = os.path.dirname( curr_DIR )
            rDIR = os.path.join( parent_DIR, os.path.basename( self.rdir.rstrip('/') ) )
            settings['results_dir'] = rDIR
            self.rdir = settings['results_dir']
            RBNS_utils.make_dir(self.rdir)
            print self.rdir
        self.edir = os.path.join( settings['results_dir'], "errors" )
        RBNS_utils.make_dir(self.edir)
        self.check_barcode_lens()
        self.barcodes = self.settings['barcodes']
        self.check_barcodes_are_separated()

        if not 1 == len(set(map(len, [self.barcodes,
            settings['concentrations']]))):
            print 'not all library descriptions are the same length'
            print 'barcodes: %i' % len(self.barcodes)
            print 'concentrations: %i' % len(settings['concentrations'])
            raise ValueError('bad input')

        self.protein_name = settings['protein_name']
        if settings['stream_count'] and not set(settings['ks_to_test_stream']).issubset(
          settings['ks_to_test_naive']):
            raise ValueError('All ks for stream must also be in naive')

        #### If conc_for_mostenriched_analyses is passed in
        try:
            self.settings['conc_for_mostenriched_analyses'] =\
                    float( settings['conc_for_mostenriched_analyses'].strip('"').strip("'") )
        except KeyError:
            self.settings['conc_for_mostenriched_analyses'] = None

        ##### Add the command line executable weblogo to the $PATH, so it
        ####    can be called from the command line when making logos
        if ( self.settings['weblogo_path'] != '' ) and\
                self.settings['weblogo_path'] not in os.environ["PATH"]:
            os.environ["PATH"] = self.settings['weblogo_path'] +\
                    ":" + os.environ["PATH"]

        print os.environ["PATH"]

        print "SETTINGS:\n"
        pprint.pprint( settings )
        print "\n\n"


    def check_barcode_lens( self ):
        """
        verifies that all the barcodes are the same length
        """
        barcode_lens = set(map(len, self.settings['barcodes']))
        if 1 != len(barcode_lens):
            raise ValueError('all barcodes must be the same length')
        self.barcode_len = barcode_lens.pop()
        self.settings['barcode_len'] = self.barcode_len

    def check_barcodes_are_separated( self,
            min_hamming_distance = 2):
        """
        - Makes sure the barcodes are all totally distinguishable
        """
        for b1, b2 in itertools.combinations(self.settings['barcodes'], 2):
            hamming_dist = RBNS_utils.hamming_distance(b1, b2)
            if hamming_dist < min_hamming_distance:
                raise ValueError('The barcodes supplied are not well '
                  'separated: %s-%s' % (b1, b2))





class RBNS_lib_settings:
    def __init__(self, experiment_settings, barcode, conc):
        self.experiment_settings = experiment_settings
        self.barcode = barcode
        self.conc = conc
        if (barcode == self.experiment_settings.get_input_barcode()):
            self.conc_for_fastq = "input"
        else:
            self.conc_for_fastq = '{0:.3f}'.format(conc).rstrip("0").rstrip(".")

    def is_input(self):
        return self.barcode == self.experiment_settings.get_input_barcode()

    def is_0nM(self):
        try:
            return self.barcode == self.experiment_settings.get_0nM_barcode()
        except KeyError:
            return False

    def return_annot(self):
        if (self.is_input()):
            return "Input"
        else:
            conc = self.get_conc()
            if (int(conc) == conc):
                return "{} nM".format( int(conc) )
            else:
                return "{0:.1f} nM".format( conc )

    def return_conc_annot_like_320_nM(self):
        return self.return_annot().replace(" ", "_")

    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_barcode(self):
        return self.barcode

    def get_conc(self):
        return self.conc

    def counts_file(self, count_type, k):
        counts_file = os.path.join(
          self.experiment_settings.get_rdir(),
          'counts',
          count_type,
          '{0}_{1}_{2}mers.pkl'.format(
                self.experiment_settings.get_property('protein_name'),
                self.conc_for_fastq,
                int( k ) ) )
        return counts_file

    def counts_exist(self, count_type, k):
        return RBNS_utils.file_exists( self.counts_file(count_type, k) )

    def get_split_reads(self):
        split_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'split_reads',
          '%(protein)s_%(conc_for_fastq)s.reads' %
           {'protein': self.experiment_settings.get_property('protein_name'),
          'conc_for_fastq': self.conc_for_fastq})
        return split_reads

    def get_split_reads_wrong_insert_len(self):
        split_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'split_reads',
          '%(protein)s_%(conc_for_fastq)s.wrong_insert_len.reads' %
           {'protein': self.experiment_settings.get_property('protein_name'),
          'conc_for_fastq': self.conc_for_fastq})
        return split_reads

    def get_conc_string(self):
        return self.conc_for_fastq

    def get_split_fastqs(self):
        try:
            conc_string = int(self.conc_for_fastq)
        #### if it's the 'input' conc.
        except ValueError:
            conc_string = self.conc_for_fastq
        split_reads_fastq = os.path.join(
          self.experiment_settings.get_rdir(),
          'split_reads',
          '%(protein)s_%(conc)s.fastq.gz' %
           {'protein': self.experiment_settings.get_property('protein_name'),
          'conc': conc_string})
        return split_reads_fastq

    def split_reads_exist(self):
        split_reads = self.get_split_reads()
        return RBNS_utils.file_exists(split_reads)

    def split_fastqs_exist(self):
        split_fastqs = self.get_split_fastqs()
        return RBNS_utils.file_exists(split_fastqs)


