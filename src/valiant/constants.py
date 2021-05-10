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

from .enums import TargetonMutator

DNA_ALPHABET = 'ACGT'
DNA_ALPHABET_SET = set(DNA_ALPHABET)

DEFAULT_CODON_TABLE_FILE_NAME = 'default_codon_table.csv'

REVCOMP_OLIGO_NAME_SUFFIX = '_rc'
CUSTOM_MUTATOR = 'custom'

DEFAULT_OLIGO_MAX_LENGTH = 300

METADATA_FIELDS = [
    'oligo_name',
    'species',
    'assembly',
    'gene_id',
    'transcript_id',
    'src_type',
    'ref_chr',
    'ref_strand',
    'ref_start',
    'ref_end',
    'revc',
    'ref_seq',
    'pam_seq',
    'vcf_alias',
    'vcf_var_id',
    'mut_position',
    'ref',
    'new',
    'ref_aa',
    'alt_aa',
    'mut_type',
    'mutator',
    'oligo_length',
    'mseq'
]

METADATA_FIELDS_SET = set(METADATA_FIELDS)

SRC_TYPES = [
    'ref',
    'cdna'
]

GENERIC_MUTATORS = frozenset([
    TargetonMutator.SNV,
    TargetonMutator.DEL1,
    TargetonMutator.DEL2_0,
    TargetonMutator.DEL2_1
])

CDS_ONLY_MUTATORS = frozenset([
    TargetonMutator.IN_FRAME,
    TargetonMutator.SNV_RE,
    TargetonMutator.STOP,
    TargetonMutator.ALA,
    TargetonMutator.AA
])

MUTATOR_CATEGORIES = sorted([m.value for m in (GENERIC_MUTATORS | CDS_ONLY_MUTATORS)])
