import abc
from typing import Optional, Tuple, Type, Union
from typing_extensions import Annotated, Literal

from pydantic import BaseModel, Field

from valiant.errors import InvalidConfig

from .. import __version__ as APP_VERSION, __package__ as APP_NAME
from ..enums import ExecMode
from .cdna_config import CDNAConfig
from .config import BaseConfig
from .sge_config import SGEConfig


class BaseMainConfig(BaseModel, abc.ABC):
    app_name: Literal['valiant']
    app_version: str = APP_VERSION
    mode: Literal[ExecMode.SGE, ExecMode.CDNA]
    params: BaseConfig

    class Config:
        allow_population_by_field_name = True
        fields = {
            'app_name': 'appName',
            'app_version': 'appVersion',
            'mode': 'mode',
            'params': 'params'
        }

    def write(self, fp: str) -> None:
        with open(fp, 'w') as fh:
            fh.write(self.json(by_alias=True, separators=(',', ':')))


class SGEMainConfig(BaseMainConfig):
    mode: Literal[ExecMode.SGE]
    params: SGEConfig


class CDNAMainConfig(BaseMainConfig):
    mode: Literal[ExecMode.CDNA]
    params: CDNAConfig


MainConfig = Annotated[Union[SGEMainConfig, CDNAMainConfig], Field(discriminator='mode')]


class MainConfigLoader(BaseModel):
    config: MainConfig


def get_main_config_from_config(config: BaseConfig) -> Optional[BaseMainConfig]:
    t: Optional[Tuple[Type, ExecMode]] = (
        (SGEMainConfig, ExecMode.SGE) if isinstance(config, SGEConfig) else
        (CDNAMainConfig, ExecMode.CDNA) if isinstance(config, CDNAConfig) else
        None
    )

    if t is None:
        raise InvalidConfig("Invalid parameters type!")

    cls, mode = t

    return cls(app_name=APP_NAME, params=config, mode=mode)
