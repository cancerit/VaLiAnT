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

"""
Functions that may force the program to exit, meant to be called by the CLI commands.
"""

import logging
import sys
from typing import Optional
from .constants import DEFAULT_CODON_TABLE_FILE_NAME
from .models.codon_table import CodonTable
from .utils import is_dna, get_data_file_path


def set_logger(log: str) -> None:
    log_level = log.upper()
    if log_level not in logging._nameToLevel:
        logging.critical("Invalid log level '%s'!" % log)
        sys.exit(1)
    logging.basicConfig(level=logging._nameToLevel[log_level])


def validate_adaptor(adaptor: Optional[str]) -> None:
    if adaptor and not is_dna(adaptor):
        logging.critical("Invalid adaptor sequence '%s'!" % adaptor)
        sys.exit(1)


def load_codon_table(fp: Optional[str]) -> CodonTable:
    if not fp:
        logging.info("Codon table not specified, the default one will be used.")

    codon_table_fp: str = fp or get_data_file_path(
        DEFAULT_CODON_TABLE_FILE_NAME)

    logging.debug("Loading codon table...")
    try:
        return CodonTable.load(codon_table_fp)
    except ValueError as ex:
        logging.critical(ex.args[0])
        logging.critical("Failed to load codon table!")
        sys.exit(1)
