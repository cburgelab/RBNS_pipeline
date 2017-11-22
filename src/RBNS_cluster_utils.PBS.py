#!/usr/bin/env python
import os, sys, pprint
import subprocess

import RBNS_utils


#### NOTE: these functions were written for a TORQUE / PBS scheduling manager
####    (e.g., submitting jobs to a queue through qsub)

####    - Compatability with SLURM workload manager may be forthcoming



def launch(
        command,
        script_options = False,
        q = 'short',
        jobname = 'ajob',
        out_file = '',
        error_DIR = '~',
        username = 'pfreese',
        ppn = '1',
        mem = '1gb' ):
    """
    - Launches the command on the cluster using the given script options
    """
    error_F = os.path.join( error_DIR, "{}.error".format( jobname ) )
    output_F = os.path.join( error_DIR, "{}.out".format( jobname ) )

    print 'will launch: '
    print command
    if not script_options:
        script_options = {'nodes': '1', 'ppn': str(ppn),
          'jobname': jobname,
          'queue': q, 'workingdir': os.getcwd()}
    script_options['error_DIR'] = "wiley:" + os.path.expanduser(error_DIR)
    script_options['command'] = command
    script_options['username'] = username
    script_options['output_F'] = output_F
    script_options['error_F'] = error_F
    outtext = """#!/bin/bash
    #PBS -m a
    #PBS -M %(username)s@mit.edu
    #PBS -N %(jobname)s
    #PBS -q %(queue)s
    #PBS -e %(error_F)s
    #PBS -o %(output_F)s
    #PBS -S /bin/bash
    echo $HOSTNAME
    echo Working directory is %(workingdir)s
    cd %(workingdir)s
    %(command)s 1> %(output_F)s 2> %(error_F)s
    echo "===== command finished =====" """ % script_options
    call = "qsub -"
    #error_f = open( os.path.join( error_DIR, "{}.error".format( jobname )), "w" )
    #output_f = open( os.path.join( error_DIR, "{}.out".format( jobname )), "w" )
    qsub = subprocess.Popen(
        call,
        shell = True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        stdin = subprocess.PIPE )
    qsub.stdin.write( outtext )

    out_F = os.path.join( error_DIR, "{}.output".format( jobname ) )
    open( out_F, 'w').write( outtext )
    output, error = qsub.communicate()

    print "output is: {}".format( output )

    if output.strip().endswith('.coyote.mit.edu'):
        job_id = int(output.strip().split('.')[0])
        print job_id, "\n"
        return job_id
    else:
        print output
        print outtext
        print 'Failure to launch'
        raise ValueError('Failure to launch')



def launch_counter(
        lib_settings,
        count_type,
        k,
        error_DIR ):
    """
    - Launches a job to perform kmer counts of count_type, calling the
        'counter' function below
    """
    split_reads = lib_settings.get_split_reads()
    out_pkl = lib_settings.counts_file( count_type, k )
    RBNS_utils.make_dir( os.path.dirname(out_pkl) )
    cluster_python_script = os.path.abspath( __file__ )
    barcode = lib_settings.get_barcode()
    out_file = os.path.join(error_DIR, 'count.%s.%s.%i.out' % (barcode, count_type, k))
    err_file = os.path.join(error_DIR, 'count.%s.%s.%i.err' % (barcode, count_type, k))
    command = ('hostname ; python %(cluster_python_script)s '
               'counter '
               '%(count_type)s '
               '%(split_reads)s '
               '%(k)i '
               '%(out_pkl)s '
               '1> %(out_file)s '
               '2> %(err_file)s ' % locals())
    conc = lib_settings.get_conc()
    jobname = '%s.%s.%i.%g' % (os.path.basename(split_reads), count_type, k, conc)
    return launch(
            command,
            jobname = jobname,
            ppn='1',
            error_DIR = error_DIR )





if __name__ == '__main__':
    fxn = sys.argv[1]
    args = ['"' + arg + '"' for arg in sys.argv[2:]]
    python_command = fxn+'('+','.join(args) + ')'
    eval( python_command )



