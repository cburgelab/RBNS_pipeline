#!/usr/bin/env python
import os, pprint
import operator
import itertools
import gzip
import subprocess
import numpy as np
#import time
#from scipy import stats
from math import pow as power
from math import sqrt, log
import cPickle as pickle
import datetime




def file_exists( F, min_size = 0 ):
    """
    makes sure a given file exists
    """
    if not os.path.exists( F ):
        return False
    fstats = os.stat( F )
    if fstats.st_size < min_size:
        return False
    if not os.access( F, os.R_OK ):
        raise ValueError( 'Input File %s cannot be read' % F )
    return True




def make_dir( dirname ):
    """
    Makes the directory; doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
            return True
        except:
            print 'The directory was made by another thread extremely recently.'
            return False
    return True


def aopen( F, mode = 'r' ):
    if ( F[-3:] == '.gz' ):
        return gzip.open( F, mode + 'b')
    else:
        return open( F, mode)



def hamming_distance(str1, str2):
    """
    - Returns the Hamming distance of str1 and str2
    """
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))



def iterNlines(
        inFile,
        N,
        strip_newlines = True ):
    """
    yields N lines at a time
    """
    assert ( N >= 1 )
    f = aopen(inFile)
    if strip_newlines:
        lines = [f.readline().strip() for i in range(N)]
    else:
        lines = [f.readline() for i in range(N)]
    try:
        while True:
            yield lines
            if strip_newlines:
                lines = [f.readline().strip() for i in range( N )]
            else:
                lines = [f.readline() for i in range(N)]
            if lines[0] == '':
                break
    except IOError: pass
    f.close()



def get_barcode(
        line,
        begin_barcode_symb,
        end_barcode_symb ):
    """
    - Extracts the barcode from the first line of a fastq quartet
        e.g., if begin_barcode_symb = #
                 end_barcode_symb = /
            @D5FF8JN1:4:1101:1220:2099#ACTTGA/1

            returns: ACTTGA
    """
    return line.split(begin_barcode_symb)[-1].split(end_barcode_symb)[0]




def yield_kmers(k):
    """
    An iterater to all kmers of length k in alphabetical order
    """
    bases = 'ACGT'
    for kmer in itertools.product( bases, repeat = k ):
        yield ''.join( kmer )



def normalize_D(
        D,
        vals_sum_to = 1.):
    """
    - Returns a normalized D, in which the values all add up to vals_sum_to
    """
    D_sum = float( sum( D.values() ) )
    norm_D = {}
    for key, val in D.iteritems():
        try:
            norm_D[key] = (val / D_sum) * vals_sum_to
        except ZeroDivisionError:
            norm_D[key] = 0.
    return norm_D




def get_kmer_from_index( kmax, index ):
    """
    takes a number (base 4)
    and returns the kmer it corresponds to in alphabetical order
    eg.
    AAAA = 0*1
    CA = 4*4 + 0*1
    GC = 3*4 + 1 * 1
    """
    bases = 'ACGT'
    out = ''
    for k in range( kmax - 1, -1, -1 ):
        face, index = divmod( index, 4 ** k )
        out += bases[face]
    return out


def kmer_in_adapters(
        kmer,
        adapters_L ):
    """
    - Returns TRUE if kmer can be found in one of the adapters in adapters_L;
        otherwise, returns False

    - Generally, the adapters_L passed in are the REVERSE-COMPLEMENT of
        the RNA adapters (i.e., they match the DNA T7 adapter sequences),
        since we're interested in possible folding between the RNA adapters
        and the random region
    """
    for adapter in adapters_L:
        if ( adapter.find( kmer ) != -1 ):
            return True
    return False


def return_Zscore(
        val,
        mean,
        st_dev ):
    """ Returns the Z-score """
    try:
        Zscore = ( val - mean ) / st_dev
    except ZeroDivisionError:
        Zscore = 0.
    return Zscore


def mean_std(
        L,
        denom = 0 ):
    """
    - Takes a list L and returns a tuple of (mean, standard deviation)

    - If denom = 0, std will have a denominator of n. If denom = -1, it will
            have a denominator of n-1
    """
    try:
        L = [float(x) for x in L]
        avg = sum(L)/len(L)
        squared = [power(x-avg,2) for x in L]
        std = sqrt(sum(squared) / (len(L) + denom))
    except ZeroDivisionError:
        print "NO ENTRIES IN LIST!"
        avg = None
        std = None
    return avg, std



def kmer_to_ignore_in_kmer(
        kmer,
        kmers_to_ignore_L ):
    """
    - Returns False if none of the kmers_to_ignore_L are contained within
        the kmer; otherwise, returns True
    """
    for sub_kmer in kmers_to_ignore_L:
        if (kmer.find( sub_kmer ) != -1):
            return True
    return False


def pkl_vals_by_kmer_D_w_formatfile(
        vals_by_kmer_D,
        abs_pkl_F ):
    """
    - Pickles the dictionary with a .format file detailing the contents
    """
    with open( abs_pkl_F, "wb" ) as f_pkl:
        pickle.dump( vals_by_kmer_D, f_pkl )

    #### make a format file
    kmer_vals_T_L = [(kmer, vals_by_kmer_D[kmer]) for kmer in vals_by_kmer_D]
    kmer_vals_T_L.sort( key = lambda x: -1 * x[1] )
    format_F = abs_pkl_F + ".format"
    with open( format_F, "w" ) as f_format:
        for kmer, val in kmer_vals_T_L:
            f_format.write( "{0}\t{1:.3f}\n".format( kmer, val ) )




def get_index_from_kmer(kmer):
    """
    returns the base10 version of the base 4 DNA representation
    """
    index = 0
    base2face = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(kmer):
        if not base in 'ACGT':
            return -1
        power = len(kmer) - 1 - i
        index += base2face[base] * (4 ** power)
    return index



def B_factor(R, k, read_len):
    """
    - Implements the B-factor formula (eq. 18 of the
        Supp. Material of Lambert et al. Mol Cell 2014)
    """
    return ((4. ** k) - 1) * (read_len - k + 1 ) * (R - 1) /\
            ((4. ** k) + read_len - k - (R * (read_len - k + 1)))




def yield_kmers_not_containing_kmer_to_ignore(
        k,
        kmers_to_ignore_L ):
    """
     An iterater to all kmers of length k in alphabetical order, but returning
        only those kmers that DON'T contain any of the kmers_to_ignore_L
    """
    for kmer in yield_kmers( k ):
        found_any = False
        for sub_kmer in kmers_to_ignore_L:
            if (kmer.find( sub_kmer ) != -1):
                found_any = True
                break
        if (found_any == False):
            yield kmer




def KL_divergence(
        p_L,
        q_L ):
    """
    - Takes two lists (must be the same length and composed of int/floats)
        and returns the KL divergence: KL(p||q)
    """
    sum_p = float( sum( p_L ) )
    sum_q = float( sum( q_L ) )

    p = [x/sum_p for x in p_L]
    q = [x/sum_q for x in q_L]

    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)

    return np.sum(np.where(p != 0, p * np.log(p / q), 0))





def return_num_lines_in_F( F ):
    """
    - Returns an integer of the number of lines in the file F
    """
    if ( F[-3:] == '.gz' ):
        p = subprocess.Popen(['zgrep', '-Ec', '$', F],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE )
    else:
        p = subprocess.Popen(['wc', '-l', F],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE )
    output, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(output.strip().split()[0])


month_str_by_num_D = {
        "01": "Jan",
        "02": "Feb",
        "03": "Mar",
        "04": "Apr",
        "05": "May",
        "06": "Jun",
        "07": "Jul",
        "08": "Aug",
        "09": "Sep",
        "10": "Oct",
        "11": "Nov",
        "12": "Dec" }



def return_nice_datetime_str_for_filename():
    """ '2015_02_15__14h_00m_36s' """
    month_str = datetime.datetime.now().strftime("%m")

    return datetime.datetime.now().strftime("%Y_{0}_%d__%Hh_%Mm_%Ss".format(
        month_str_by_num_D[month_str] ) )




def split_splitreadsF_blockidx_T_L_into_lists_max_X(
        splitreadsF_blockidx_T_L,
        num_max_jobs_at_one_time ):
    """
    - Splits up the list of tuples in splitreadsF_blockidx_T_L into multiple
        lists, each of which is max length num_max_jobs_at_one_time

    - Called by fold_each_reads_by_block_F() in RBNS_main.py
    """
    splitreadsF_blockidx_T_Ls_L = []

    this_L = []
    for T in splitreadsF_blockidx_T_L:
        this_L.append( T )
        if ( len( this_L ) == num_max_jobs_at_one_time ):
            splitreadsF_blockidx_T_Ls_L.append( this_L )
            this_L = []

    if ( len( this_L ) > 0 ):
        splitreadsF_blockidx_T_Ls_L.append( this_L )

    return splitreadsF_blockidx_T_Ls_L


def copy_lines_lower_through_upper_to_another_F_without_N(
        orig_F,
        out_F,
        lower_idx,
        upper_idx ):
    """
    - Copies all lines that don't contain 'N' from orig_F to out_F
    """
    w_N_F = out_F + ".w_N"
    #### Now remove the N's
    remove_N_cmd = "grep -v 'N' {0} > {1}".format( orig_F, w_N_F )
    os.system( remove_N_cmd )

    lines_str = '"{0},{1} p"'.format( lower_idx + 1, upper_idx )
    cmd = 'sed -n {0} {1} > {2}'.format( lines_str, w_N_F, out_F )
    os.system( cmd )

    os.system( "rm {}".format( w_N_F ) )



def return_num_lines_in_F_WITHOUT_pattern( F, patt ):
    """
    - Returns an integer of the number of lines in the file F
    """
    if ( F[-3:] == '.gz' ):
        num_lines = 0
        with aopen( F ) as f:
            for line in f:
                if ( line.find( patt ) == -1 ):
                    num_lines += 1
        return num_lines
    else:
        p = subprocess.Popen(['grep', '-v', '-c', patt, F],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE )
        output, err = p.communicate()
        if p.returncode != 0:
            raise IOError(err)
        return int(output.strip().split()[0])









