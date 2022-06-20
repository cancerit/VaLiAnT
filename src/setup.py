#!/usr/bin/env python3

########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022 Genome Research Ltd
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#############################

from setuptools import setup, find_packages

PACKAGE_NAME = 'valiant'
VERSION = '2.0.0'

config = {
    'name': PACKAGE_NAME,
    'description': 'Variant Library Annotation Tool',
    'author': 'Luca Barbon',
    'author_email': 'cgphelp@sanger.ac.uk',
    'version': VERSION,
    'python_requires': '>= 3.7',
    'install_requires': [
        'chardet',
        'click',
        'cython',
        'pandas>=1.1,<1.2',
        'pyranges',
        'pysam'
    ],
    'tests_require': ['pytest'],
    'packages': find_packages(),
    'package_data': {
        PACKAGE_NAME: ['data/*']
    },
    'entry_points': {
        'console_scripts': [f'{PACKAGE_NAME}={PACKAGE_NAME}.__main__:main']
    }
}

setup(**config)
