#! /usr/bin/env python3
"""
Tests srclist_by_beam.py script
"""
import os
import filecmp
import subprocess

def test_srclist_by_beam():
    """Test the entire srclist_by_beam.py script by running it and comparing outputs"""
    obsids = [1173782232, 1205099768, 1224243968, 1240826896, 1253482168]
    source_numbers = [100, 1000]
    source_lists = ['srclist_pumav3_EoR0aegean_EoR1pietro+ForA_phase1+2.txt',
                    'srclist_pumav3_EoR0aegean_EoR1pietro+ForA_TGSSgalactic.txt',
                    'srclist_puma-v2_complete.txt',
                    'srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt']
    
    # Loop over all inputs then compare the files with reference files
    for obs in obsids:
        for sn in source_numbers:
            for slist in source_lists:
                # Run srclist
                subprocess.call(['python3', 'srclist_by_beam.py', '-m',
                                 'tests/test_files/{}.metafits'.format(obs), '-n',
                                 str(sn), '-s', slist])

                # Compare the outputs
                output_filename = '{}_{}_patch{}.txt'.format(slist.split('.')[0], obs, sn)
                if obs == 1240826896:
                    # No source above 10.00Jy within initial cutoff distance
                    # so make sure no file is made
                    if os.path.isfile(output_filename):
                        raise AssertionError()
                else:
                    if not filecmp.cmp(output_filename,
                                       "tests/test_files/test_{}".format(output_filename)):
                        raise AssertionError()
                    os.remove(output_filename)



if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
