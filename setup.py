#! /usr/bin/env python3
"""
Setup for srclists
"""
import os
from setuptools import setup


def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name="srclists",
      url="https://github.com/JLBLine/srclists",
      #description="",
      long_description=read('README.md'),
      scripts=['srclist_by_beam.py',
               'srclist_by_beam_lib.py',
               'srclist_pumav3_EoR0aegean_EoR1pietro+ForA_phase1+2.txt',
               'srclist_pumav3_EoR0aegean_EoR1pietro+ForA_TGSSgalactic.txt',
               'srclist_puma-v2_complete.txt',
               'srclist_pumav3_EoR0aegean_EoR1pietro+ForA.txt'],
      #data_files=[],
      setup_requires=['pytest-runner'],
      tests_require=['pytest']
)

