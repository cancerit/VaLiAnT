########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020, 2021, 2022, 2023 Genome Research Ltd
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

import csv
from typing import Generator

from .utils import detect_encoding


def load_csv(fp: str, columns: list[str] | None = None, delimiter: str = ',') -> Generator[list[str], None, None]:
    """Validate header and parse all rows of a CSV or TSV file as lists of strings"""

    with open(fp, encoding=detect_encoding(fp)) as fh:
        reader = csv.reader(fh, delimiter=delimiter)

        # Validate header
        if columns is not None and next(reader) != columns:
            t = 'TSV' if delimiter == '\t' else 'CSV'
            raise ValueError(f"Invalid {t} header!")

        # Load rows
        for r in reader:
            yield r
