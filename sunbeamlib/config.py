import os
import sys
import collections
from pathlib import Path
from pkg_resources import resource_stream

from semantic_version import Version
import ruamel.yaml
from sunbeamlib import __version__

def makepath(path):
    return Path(path).expanduser()


def verify(path):
    path = Path(path)
    if path.exists():
        return path.resolve()
    else:
        raise ValueError(f"Path {path} does not exist")
    
    
def validate_paths(cfg, root):
    """Process paths in config file subsection.

    For each key ending in _fp, the value is:
    - converted to a pathlib.Path
    - ensured to be an absolute path, by appending `root` if needed
    - ensured to exist, if it is not the value from `output_fp`
    - expanded home directory ~
    
    :param cfg: a config file subsection
    :returns: an updated copy of cfg
    """
    new_cfg = {}
    for k, v in cfg.items():
        if k.endswith('_fp'):
            v = makepath(v)
            if not v.is_absolute():
                v = root/v
            if k != 'output_fp':
                try: v = verify(v)
                except ValueError:
                    raise ValueError(
                        "For key '%s': path '%s' does not exist" % (k,v))
        new_cfg[k] = v
    return new_cfg

def check_compatibility(cfg):
    """Returns the major version numbers from the package and config file, respectively"""

    cfg_version = Version(cfg['all'].get('version', '0.0.0'))
    pkg_version = Version(__version__)

    return (pkg_version.major, cfg_version.major)
    
def check_config(cfg):
    """Resolve root in config file, then validate paths."""
    
    root = verify(cfg['all']['root']) if 'root' in cfg['all'] else Path.cwd()
    # Iteratively check paths for each subsection
    new_cfg = {
        section: validate_paths(values, root)
        for section, values in cfg.items()
    }

    new_cfg['all']['root'] = root
    return new_cfg


def output_subdir(cfg, section):
    return cfg['all']['output_fp']/cfg[section]['suffix']


def process_databases(db_dict):
    """Process the list of databases.

    Expands the nucleotide and protein databases specified
    """
    dbs = {'nucl':{}, 'prot':{}}
    root = verify(makepath(db_dict['root_fp']))
    nucl = db_dict.get('nucleotide')
    prot = db_dict.get('protein')
    if nucl:
        dbs['nucl'] = {db: str(root/path) for db, path in nucl.items()}
    if prot:
        dbs['prot'] = {db: str(root/path) for db, path in prot.items()}
    return dbs


def _update_dict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.Mapping):
            # We could use .get() here but ruamel.yaml's weird Mapping
            # subclass outputs errors to stdout if the key doesn't exist
            target[k] = _update_dict(target[k], v) if k in target else _update_dict({}, v)
        else:
            target[k] = v
    return target

def _update_dict_strict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.Mapping) and k in target.keys():
            target[k] = _update_dict_strict(target.get(k, {}), v)
        elif k in target.keys():
            target[k] = v
        else:
            sys.stderr.write("Key '%s' not found in target, skipping\n" % k)
            continue
    return target

def update(config_str, new, strict=False):
    config = ruamel.yaml.round_trip_load(config_str)
    if strict:
        config = _update_dict_strict(config, new)
    else:
        config = _update_dict(config, new)
        if sbx_config := ruamel.yaml.round_trip_load(extension_config()):
            config = _update_dict(config, sbx_config)
    return config

def new(
        project_fp,
        version=__version__,
        template=None):
    if template:
        config = template.read()
    else:
        config = str(resource_stream(
            "sunbeamlib", "data/default_config.yml").read().decode())
        # add config from extensions
        config += extension_config()

    return config.format(
        PROJECT_FP=project_fp,
        SB_VERSION=version)

def extension_config():
    config = ""
    sunbeam_dir = Path(os.getenv("SUNBEAM_DIR", os.getcwd()))
    for sbx in os.listdir(sunbeam_dir/"extensions"):
        try:
            sbx_files = os.listdir(sunbeam_dir/"extensions"/sbx)
        except NotADirectoryError:
            continue
        if 'config.yml' in sbx_files:
            # append it to the existing config
            sbx_config_fp = sunbeam_dir/"extensions"/sbx/"config.yml"
            with open(sbx_config_fp) as sbx_configfile:
                sbx_config = "\n"+sbx_configfile.read()
            config = str(config + sbx_config)
    return config

def load_defaults(default_name):
    return ruamel.yaml.safe_load(
        resource_stream("sunbeamlib", f"data/{default_name}.yml")
        .read()
        .decode()
    )
    
def dump(config, out=sys.stdout):
    if isinstance(config, collections.Mapping):
        ruamel.yaml.round_trip_dump(config, out)
    else:
        out.write(config)
