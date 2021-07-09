#!/usr/bin/env python3

########## LICENCE ##########
# VaLiAnT, (c) 2020, GRL (the "Software")
#
# The Software remains the property of Genome Research Ltd ("GRL").
#
# The Software is distributed "AS IS" under this Licence solely for non-commercial use in the hope that it will be useful,
# but in order that GRL as a charitable foundation protects its assets for the benefit of its educational and research
# purposes, GRL makes clear that no condition is made or to be implied, nor is any warranty given or to be implied, as to
# the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific
# conditions. Furthermore, GRL disclaims all responsibility for the use which is made of the Software. It further
# disclaims any liability for the outcomes arising from using  the Software.
#
# The Licensee agrees to indemnify GRL and hold GRL harmless from and against any and all claims, damages and liabilities
# asserted by third parties (including claims for negligence) which arise directly or indirectly from the use of the
# Software or the sale of any products based on the Software.
#
# No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic
# or mechanical, without the express permission of GRL. The permission of GRL is not required if the said reproduction,
# modification, transmission or transference is done without financial return, the conditions of this Licence are imposed
# upon the receiver of the product, and all original and amended source code is included in any transmitted product. You
# may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by
# these terms and conditions.
#
# You are not permitted under this Licence to use this Software commercially. Use for which any financial return is
# received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the
# Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software
# or any derivative of it for research with the final aim of developing software products for sale or license to a third
# party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software
# products for sale or license to a third party, or (4) use of the Software to provide any service to an external
# organisation for which payment is received. If you are interested in using the Software commercially, please contact
# legal@sanger.ac.uk. Contact details are: legal@sanger.ac.uk quoting reference Valiant-software.
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
