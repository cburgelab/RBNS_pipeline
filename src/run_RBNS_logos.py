#!/usr/bin/env python
import os, sys
import subprocess
import pprint

###### A helper function called within RBNS_main that directs a call to the
#####   logo pipeline


def run_multiple_logos(
        multiple_runs_txt_F,
        config_Fs_DIR,
        RBNS_logos_script_basename = "RBNS_logos.py" ):
    """
    - INPUT:
        - a .txt file containing a separate config_F basename (with
            all config_F's needing to be in the config_Fs_DIR) and its options
            on a separate line, e.g.:
            contains the lines:
                MBNL1.config_F --num-reads 2000000 --starting-k 6 --Zscore-kmers-to-keep 3.
                CELF1.config_F --num-reads 2000000 --starting-k 7 --Zscore-kmers-to-keep 3.
                RBFOX2.config_F --num-reads 2000000 --starting-k 5 --Zscore-kmers-to-keep 3.

    - Note that the argument after --num-reads can be either an int (if
        the # of reads to use) or a float (a proportion of the total # reads in
        the file)
    """
    config_Fs_DIR = config_Fs_DIR.rstrip( "/" )

    RBNS_logos_scripts_DIR = os.path.dirname(
            os.path.abspath(__file__) )
    RBNS_logos_script = os.path.join( RBNS_logos_scripts_DIR,
            RBNS_logos_script_basename )

    with open( multiple_runs_txt_F ) as f:
        for line in f:
            if (len( line.strip() ) == 0):
                continue
            config_basename = line.strip().split( " " )[0]
            print config_basename
            args_L = line.strip().split( " " )[1:]
            pprint.pprint( args_L )
            print "\n\n"
            cmd = 'python {0} {1} {2}'.format(
                RBNS_logos_script,
                os.path.join( config_Fs_DIR, config_basename ),
                " ".join( args_L ) )
            print cmd
            os.system( cmd )





if __name__ == '__main__':

    multiple_runs_txt_F = sys.argv[1]
    config_Fs_DIR = sys.argv[2]

    run_multiple_logos(
            multiple_runs_txt_F,
            config_Fs_DIR )





