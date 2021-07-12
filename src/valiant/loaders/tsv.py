########## LICENCE ##########
# VaLiAnT
# Copyright (C) 2020-2021 Genome Research Ltd
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

import logging
from typing import List
import chardet
import csv


def load_tsv(fp: str, columns: List[str]) -> List[List[str]]:

    # Detect encoding
    with open(fp, 'rb') as rfh:
        encoding: str = chardet.detect(rfh.read(10000))['encoding']

    logging.debug("File '%s' encoding: %s." % (fp, encoding))

    # Load rows
    with open(fp, encoding=encoding) as fh:
        reader = csv.reader(fh, delimiter='\t')
        if next(reader) != columns:
            raise ValueError("Invalid TSV header!")
        return [r for r in reader]
