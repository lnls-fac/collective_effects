#!/usr/bin/env python3

from setuptools import setup

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='cppcolleff',
    version=__version__,
    author='lnls-fac',
    description='cppcolleff python package',
    url='https://github.com/lnls-fac/collective_effects',
    download_url='https://github.com/lnls-fac/collective_effects',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['cppcolleff'],
    package_data={'cppcolleff': ['_cppcolleff.so', 'VERSION']},
    zip_safe=False
)
