#!/usr/bin/env python
import os, sys, pprint
import subprocess

import RBNS_main
import RBNS_utils


#### NOTE: these functions were written for a SLURM scheduling manager,
####    with the "#SBATCH -C centos7" argument being specific to the
####    cluster on which this pipeline was developed. You may need to
####    modify the launch() function below so it is compatible with your
####    cluster

####    NOTE: compatibility with a PBS/Torque scheduler is provided in
####        the RBNS_cluster_utils.PBS.py file



def launch(
        command,
        script_options = False,
        q = 'sched_engaging_default',
        jobname = 'ajob',
        out_file = '',
        error_DIR = '~',
        username = 'kkrick',
        time_mins = 'default',
        mem = '1gb' ):
    """
    - Launches the command on the cluster using the given script options
    %(command)s 1> %(output_F)s 2> %(error_F)s
    call = "sbatch -o %(output_F)s -e %(error_F)s -p %(partition)s -J %(jobname)s -N %(nodes)s $(command)" % script_options
    """
    error_F = os.path.join( error_DIR, "{}.error".format( jobname ) )
    output_F = os.path.join( error_DIR, "{}.out".format( jobname ) )

    print 'will launch: '
    print command
    if not script_options:
        script_options = {
                'nodes': '1',
                'jobname': jobname,
                'partition': q,
                'workingdir': os.getcwd() }
    script_options['error_DIR'] = os.path.expanduser( error_DIR )
    script_options['command'] = command
    script_options['username'] = username
    script_options['output_F'] = output_F
    script_options['error_F'] = error_F

    cmd_L = [
            "#!/bin/bash",
            "#SBATCH -o %(output_F)s" % script_options,
            "#SBATCH -e %(error_F)s" % script_options,
            "#SBATCH -p %(partition)s" % script_options,
            "#SBATCH -J %(jobname)s" % script_options,
            "#SBATCH -N %(nodes)s" % script_options,
            "#SBATCH -C centos7"]
    if ( time_mins != 'default' ):
        cmd_L.append( "#SBATCH -t {}".format( time_mins ) )
    cmd_L += ["cd %(workingdir)s" % script_options,
              "%(command)s" % script_options ]

    sh_script = os.path.join( error_DIR, "{}.sh".format( jobname ) )
    print sh_script
    with open( sh_script, 'w' ) as f:
        f.write( "\n".join( cmd_L ) )

    call = "sbatch {0}".format( sh_script )
    #error_f = open( os.path.join( error_DIR, "{}.error".format( jobname )), "w" )
    #output_f = open( os.path.join( error_DIR, "{}.out".format( jobname )), "w" )
    qsub = subprocess.Popen(
        call,
        shell = True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        stdin = subprocess.PIPE )

    output, error = qsub.communicate()

    print "output is: {}".format( output )

    if ( output.strip().find( 'batch job' ) != -1 ):
        job_id = int(output.strip().split(' ')[-1])
        print job_id, "\n"
        return job_id
    else:
        print output
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
    out_file = os.path.join( error_DIR, 'count.%s.%s.%i.out' % (barcode, count_type, k) )
    err_file = os.path.join( error_DIR, 'count.%s.%s.%i.err' % (barcode, count_type, k) )
    command = ('python %(cluster_python_script)s '
               'counter '
               '%(count_type)s '
               '%(split_reads)s '
               '%(k)i '
               '%(out_pkl)s ' % locals() )
               #'1> %(out_file)s '
               #'2> %(err_file)s ' % locals())
    conc = lib_settings.get_conc()
    jobname = '%s.%s.%i.%g' % (os.path.basename(split_reads), count_type, k, conc)
    return launch(
            command,
            jobname = jobname,
            error_DIR = error_DIR )



def counter(
        count_type,
        split_reads,
        k,
        out_pkl ):
    """
    - Performs the counts of count_type on the split_reads, pickling them to
        out_pkl
    """
    k = int( k )
    if count_type == 'naive':
        RBNS_main.count_naive(split_reads, k, out_pkl)
    elif count_type == 'stream':
        RBNS_main.count_stream(split_reads, k, out_pkl)
    elif count_type == 'by_position':
        RBNS_main.count_by_position(split_reads, k, out_pkl)
    else:
        raise ValueError( 'Unknown count type: %s ' % count_type )




def test_launch():

    launch(
        "python /home/pfreese/py_lib/test.py",
        script_options = False,
        q = 'sched_mit_hill',
        error_DIR = "/net/eofe-data010/data001/burgelab/nevermind/data/nm/pfreese/ttt",
        jobname = 'tm' )


if __name__ == '__main__':
    fxn = sys.argv[1]
    args = ['"' + arg + '"' for arg in sys.argv[2:]]
    python_command = fxn+'('+','.join(args) + ')'
    eval( python_command )



