#!/usr/bin/env python
"""
Integration tests for the srclist_by_beam.py script.
"""

import os
import filecmp

import pytest

from srclists.srclist_by_beam import _parse_args, main as srclist_by_beam


# To get whatever default command-line arguments are being used in
# srclist_by_beam.py, just run the _parse_args function there. We can overwrite
# fields in the returned args for different tests.
default_args = _parse_args()

source_lists = ['srclist_pumav3_EoR0aegean_EoR1pietro+ForA_phase1+2.txt',
                'srclist_pumav3_EoR0aegean_EoR1pietro+ForA_TGSSgalactic.txt',
                'srclist_puma-v2_complete.txt',
                'srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt']

path_to_test_files_dir = os.path.dirname(os.path.realpath(__file__)) + "/test_files"

# These obsids should not generate any outputs.
bad_obsids = [1240826896]


def srclist_by_beam_shell_string(obsid, num_sources, source_list):
    return "srclist_by_beam.py -m tests/test_files/{}.metafits -n {} -s {}" \
        .format(obsid, num_sources, source_list).split(" ")


# The following tests run srclist_by_beam.py and compare the outputs to
# previously-generated outputs.
def compare_output_with_old_output(obsid, num_sources):
    for slist in source_lists:
        args = default_args
        args.srclist = slist
        args.num_sources = num_sources
        args.metafits = "{}/{}.metafits".format(path_to_test_files_dir, obsid)
        new_output = '{}_{}_patch{}.txt'.format(slist.split('.')[0], obsid, num_sources)

        if obsid not in bad_obsids:
            srclist_by_beam(args)
            old_output = "{}/test_{}".format(path_to_test_files_dir, new_output)
            try:
                assert(filecmp.cmp(new_output, old_output))
            finally:
                os.remove(new_output)
        else:
            # Obsids that don't make an output will emit a message like "No
            # source above 10.00Jy within initial cutoff distance" and run
            # sys.exit()
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                srclist_by_beam(args).will_exit_somewhere_down_the_stack()
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 1
            assert not os.path.isfile(new_output)


def test_srclist_by_beam_1173782232_100_sources():
    compare_output_with_old_output(1173782232, 100)


def test_srclist_by_beam_1173782232_1000_sources():
    compare_output_with_old_output(1173782232, 1000)


def test_srclist_by_beam_1205099768_100_sources():
    compare_output_with_old_output(1205099768, 100)


def test_srclist_by_beam_1205099768_1000_sources():
    compare_output_with_old_output(1205099768, 1000)


def test_srclist_by_beam_1224243968_100_sources():
    compare_output_with_old_output(1224243968, 100)


def test_srclist_by_beam_1224243968_1000_sources():
    compare_output_with_old_output(1224243968, 1000)


def test_srclist_by_beam_1240826896_100_sources():
    compare_output_with_old_output(1240826896, 100)


def test_srclist_by_beam_1240826896_1000_sources():
    compare_output_with_old_output(1240826896, 1000)


def test_srclist_by_beam_1253482168_100_sources():
    compare_output_with_old_output(1253482168, 100)


def test_srclist_by_beam_1253482168_1000_sources():
    compare_output_with_old_output(1253482168, 1000)


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
