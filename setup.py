#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
from itertools import chain
import sys

from Cython.Distutils import build_ext
from setuptools import setup, Extension
import pip
from pip.req import parse_requirements

import versioneer

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

base_reqs = parse_requirements(
    'requirements.txt', session=pip.download.PipSession())

cmdclass = versioneer.get_cmdclass()

# Add cyordereddict for Python <=3.5 for performance boost
if sys.version_info[:2] < (3, 6):
    pre36_reqs = parse_requirements(
        'requirements_pre36.txt', session=pip.download.PipSession())
else:
    pre36_reqs = []

requirements = [str(ir.req) for ir in chain(base_reqs, pre36_reqs)]

test_requirements = [
    str(ir.req)
    for ir in parse_requirements(
        'requirements_test.txt', session=pip.download.PipSession())
]

# Building of HTSlib Cython wrapper, chunk taken from cyvcf
excludes = ['irods', 'plugin']
sources = [x for x in glob.glob('htslib/*.c')
           if not any(e in x for e in excludes)] + glob.glob('htslib/cram/*.c')
# these have main()'s or require other libraries
sources = [x for x in sources
           if not x.endswith((
               'htsfile.c', 'tabix.c', 'bgzip.c',  # main()'s
               'hfile_s3.c', 'hfile_libcurl.c', # more dependencies
            ))]

extension = Extension('vcfpy.cyhtslib.cyhtslib',
                      ['vcfpy/cyhtslib/cyhtslib.pyx'] + sources,
                      libraries=['z'],
                      extra_compile_args=['-O0'],
                      include_dirs=['htslib', 'vcfpy/cyhtslib'])

# Update cmdclass dict
cmdclass['build_ext'] = build_ext

setup(
    name='vcfpy',
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    description=(
        'Python 3 VCF library with good support for both reading and writing'),
    long_description=readme + '\n\n' + history,
    author="Manuel Holtgrewe",
    author_email='manuel.holtgrewe@bihealth.de',
    url='https://github.com/bihealth/vcfpy',
    packages=[
        'vcfpy',
    ],
    package_dir={'vcfpy': 'vcfpy'},
    ext_modules=[extension],
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='vcfpy',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
