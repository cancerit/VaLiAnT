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

import pandas as pd
from .constants import META_REF_AA
from .mave_hgvs import MAVEPrefix, get_mave_nt


# Metadata table field used to determine whether an entry maps to a CDS
META_CDS_PROBE_FIELD = META_REF_AA


def is_metadata_row_cds(r: pd.Series) -> bool:
    return META_CDS_PROBE_FIELD in r and not pd.isna(r[META_CDS_PROBE_FIELD])


def get_mave_nt_from_row(r: pd.Series) -> str:
    return get_mave_nt(
        MAVEPrefix.LINEAR_GENOMIC,
        r.var_type,
        r.mut_position,
        r.ref,
        r.new)
