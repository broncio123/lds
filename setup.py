#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import find_packages
from distutils.core import setup

import os
on_rtd = os.environ.get('READTHEDOCS') == 'True'

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

package_name = 'ldds'
version_num = '0.1.0'

def main():
    install_requires = ['m2r'] if on_rtd else []
    setup(
        name = package_name,

        version = version_num,
        
        description = 'Python library for computing Lagrangian descriptors',
        
        long_description = readme,
        
        author = '.',
        
        author_email = '.',
        
        url = 'https://github.com/champsproject/ldds',
        
        license = license,
        
        # packages=find_packages(exclude=('test','data', 'docs')),
        packages=['ldds'],
        package_dir={'ldds': 'ldds'},

        install_requires = requirements
    )


if __name__ == '__main__':
    main()
