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

import json
import logging
import click
from typing import Optional
from . import __version__
from .cdna_cli import cdna, run_cdna
from .common_cli import existing_file
from .errors import InvalidConfig
from .models.main_config import BaseMainConfig, CDNAMainConfig, MainConfigLoader, SGEMainConfig
from .sge_cli import run_sge, sge


def load_main_config(fp: str) -> BaseMainConfig:
    with open(fp) as fh:
        try:
            config_dict = json.load(fh)
        except json.JSONDecodeError:
            raise InvalidConfig("not a JSON!")

    return MainConfigLoader.parse_obj({
        'config': config_dict
    }).config


@click.group(invoke_without_command=True)
@click.option('-c', '--config', 'config_fp', type=existing_file, help="Configuration file path")
@click.version_option(__version__)
@click.pass_context
def main(ctx: click.Context, config_fp: Optional[str]):
    if ctx.invoked_subcommand is None:
        if not config_fp:
            raise click.UsageError("Configuration required if no subcommand is specified!")

        # Load configuration
        try:
            config = load_main_config(config_fp)

        except InvalidConfig as ex:
            logging.critical(f"Invalid configuration: {ex.args[0]}")
            ctx.exit(1)

        except (PermissionError, FileNotFoundError) as ex:
            logging.critical(ex)
            ctx.exit(1)

        # Run in the appropriate mode
        if isinstance(config, SGEMainConfig):
            run_sge(config.params, False)

        elif isinstance(config, CDNAMainConfig):
            run_cdna(config.params)

        else:
            logging.critical("Invalid configuration!")
            ctx.exit(1)

        ctx.exit(0)


main.add_command(sge)
main.add_command(cdna)
