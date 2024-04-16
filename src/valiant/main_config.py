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

import abc
from typing import Type, Union

from pydantic import BaseModel, Field
from typing_extensions import Annotated, Literal

from . import __version__ as APP_VERSION, __package__ as APP_NAME
from .cdna_config import CDNAConfig
from .config import BaseConfig
from .enums import ExecMode
from .errors import InvalidConfig
from .sge_config import SGEConfig


class BaseMainConfig(BaseModel, abc.ABC):
    app_name: Literal['valiant'] = Field(alias='appName')
    app_version: str = Field(alias='appVersion', default=APP_VERSION)
    mode: Literal[ExecMode.SGE, ExecMode.CDNA] = Field()
    params: BaseConfig = Field()

    class Config:
        populate_by_name = True

    def write(self, fp: str) -> None:
        with open(fp, 'w') as fh:
            # TODO: check default delimiters...
            fh.write(self.model_dump_json(by_alias=True))


class SGEMainConfig(BaseMainConfig):
    mode: Literal[ExecMode.SGE]
    params: SGEConfig


class CDNAMainConfig(BaseMainConfig):
    mode: Literal[ExecMode.CDNA]
    params: CDNAConfig


MainConfig = Annotated[Union[SGEMainConfig, CDNAMainConfig], Field(discriminator='mode')]


class MainConfigLoader(BaseModel):
    config: MainConfig


def get_main_config_from_config(config: BaseConfig) -> BaseMainConfig | None:
    t: tuple[Type, ExecMode] | None = (
        (SGEMainConfig, ExecMode.SGE) if isinstance(config, SGEConfig) else
        (CDNAMainConfig, ExecMode.CDNA) if isinstance(config, CDNAConfig) else
        None
    )

    if t is None:
        raise InvalidConfig("Invalid parameters type!")

    cls, mode = t

    return cls(app_name=APP_NAME, params=config, mode=mode)
