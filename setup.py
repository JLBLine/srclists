#!/usr/bin/env python
"""Setup for srclists."""

import os
import re
from setuptools import setup


here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "srclists/__init__.py")) as f:
    contents = f.read()
    version_number = re.search(r"__version__ = \"(\S+)\"", contents).group(1)


def read_file(fname):
    """Read the contents of a file given a relative path."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name="srclists",
      version=version_number,
      description="A location for all of the current srclists",
      long_description=read_file("README.md"),
      url="https://github.com/JLBLine/srclists",
      author="Jack Line",
      author_email="jack.line@curtin.edu.au",
      license="MPL 2.0",
      packages=["srclists"],
      scripts=["srclists/srclist_by_beam.py"],
      install_requires=["numpy", "astropy", "future", "mwa_pb"],
      )
