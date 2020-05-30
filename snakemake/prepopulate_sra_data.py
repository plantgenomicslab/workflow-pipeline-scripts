import argparse
import logging
import os
import pprint
import sys
import tempfile

from snakemake.shell import shell
from parallel_fastq_dump import pfd

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

logger = logging.getLogger('Prepopulate SRA Data')
logger.setLevel(logging.WARNING)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(filename)s[%(lineno)d] %(funcName)s - %(message)s',
                              datefmt='%F %T')
handler.setFormatter(formatter)
logger.addHandler(handler)


def get_config_from_file(fpath):
    if not os.path.exists(fpath):
        logger.error('Configuration file not found')
        return None
    with open(fpath, 'rb') as fhandle:
        data = yaml.load(fhandle, Loader)
    if not data.get('sradata'):
        logger.error('No SRA data found in config')
        return None
    return data


def get_samples_and_units(config):
    samples = []
    sra_ids = []
    units = []
    outdir = config.get('outdir')
    for sra in config.get('sra_list'):
        label = sra.get('label')
        accession = sra.get('id')
        condition = sra.get('condition')
        unit = sra.get('unit')
        sra_ids.append(accession)
        samples.append([label, condition])
        units.append([label, unit, f"{outdir}/{accession}_1.fastq", f"{outdir}/{accession}_2.fastq"])
    return samples, sra_ids, units


def get_sradata_from_config(data):
    if not isinstance(data, dict):
        logger.error('Invalid data object: %s', type(data))
        return None
    if not data.get('sradata'):
        logger.error('No SRA data found in config')
        return None
    return data.get('sradata')


def generate_samples_and_units(samples, units):
    samples_tsv = "sample\tcondition\n"
    for sample in samples:
        samples_tsv += "\t".join(sample) + "\n"
    units_tsv = "sample\tunit\tfq1\tfq2\n"
    for unit in units:
        units_tsv += "\t".join(unit) + "\n"
    return samples_tsv, units_tsv


def populate(config, threads=0):
    if not config:
        logger.error('No configuration provided')
        return None
    sradata = get_sradata_from_config(config)
    if not sradata:
        logger.error('No SRA data found in configuration')
        return None
    samples, sra_ids, units = get_samples_and_units(sradata)
    if not (samples and sra_ids and units):
        logger.error('Failed to get samples and units data from configuration')
        return None
    if not pull_sra_data(sradata, sra_ids, threads):
        logger.error('Failed to pull SRA data')
        return False
    try:
        write_samples_and_units_tsv(samples, units)
    except Exception as err:
        logger.error('Unhandled exception attempting to write samples and units TSV files: %s', err)
        return False
    return True


def pull_sra_data(config, sra_ids, threads=0):
    if not config:
        logger.error('No configuration provided')
        return None
    if not config.get('method'):
        logger.error('No pull method specified')
        return None
    valid_methods = ('fastq', 'fasterq')
    method = config.get('method')
    if config.get('method') not in valid_methods:
        logger.error('Invalid pull method %s. Must be one of: %s', method, ', '.join(valid_methods))
        return None
    if not config.get('outdir'):
        logger.error('No output directory specified')
        return None
    outdir = config.get('outdir')
    extra = ""
    if config.get('extra'):
        extra = config.get('extra')
    if threads:
        if isinstance(threads, (bytes, str)) and threads.isnumeric():
            try:
                threads = int(threads)
            except ValueError:
                threads = ""
        elif not isinstance(threads, int):
            threads = ""
    else:
        threads = ""
    if method == 'fastq':
        args = argparse.Namespace()
        args.minSpotId = 1
        args.maxSpotId = None
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        args.outdir = outdir
        args.threads = threads if threads else 1
        if extra:
            if isinstance(extra, (str, bytes)):
                extra = extra.split()
            elif not isinstance(extra, list):
                extra = []
        else:
            extra = []
        with tempfile.TemporaryDirectory() as tmp:
            args.tmpdir = tmp
            for accession in sra_ids:
                downloaded = True
                for path in [os.path.join(outdir, f"{accession}{suffix}") for suffix in ['.fastq', '_1.fastq', '_2.fastq']]:
                    if not os.path.exists(path):
                        logger.warning('Missing SRA file %s. Re-downloading %s.', path, accession)
                        downloaded = False
                        break
                if not downloaded:
                    logger.info('Pulling SRA accession %s, threads "%s"', accession, threads)
                    pfd(args, accession, extra)
                else:
                    logger.info('All files for %s already downloaded', accession)
        return True
    elif method == 'fasterq':
        outarg = f"--outdir {outdir}"
        threads = f"--threads {threads}" if threads else ""
        with tempfile.TemporaryDirectory() as tmp:
            for accession in sra_ids:
                downloaded = True
                for path in [os.path.join(outdir, f"{accession}{suffix}") for suffix in ['.fastq', '_1.fastq', '_2.fastq']]:
                    if not os.path.exists(path):
                        logger.warning('Missing SRA file %s. Re-downloading %s.', path, accession)
                        downloaded = False
                        break
                if not downloaded:
                    logger.info('Pulling SRA accession %s, threads "%s"', accession, threads)
                    shell(
                        "fasterq-dump --temp {tmp} {threads} "
                        "{extra} {outarg} {accession}"
                    )
                else:
                    logger.info('All files for %s already downloaded', accession)
        return True
    return None


def write_samples_and_units_tsv(samples, units):
    if not (samples and units):
        logger.error('No samples and units data provided')
        return None
    samples_tsv, units_tsv = generate_samples_and_units(samples, units)
    if not (samples_tsv and units_tsv):
        logger.error('Failed to generate TSV data for samples and units')
        return None
    with open('samples.tsv', 'w') as fhandle:
        fhandle.write(samples_tsv)
    with open('units.tsv', 'w') as fhandle:
        fhandle.write(units_tsv)
    return True

def main():
    logger.setLevel(logging.INFO)
    parser = argparse.ArgumentParser(description='Prepopulate SRA Data')
    parser.add_argument('--config', type=str,
                        help='YAML configuration file to load SRA info from')
    parser.add_argument('--threads', type=int,
                        help='Number of threads to use')
    args, extra = parser.parse_known_args()
    logger.info('Args: %s, Extra: %s', args, extra)
    config = get_config_from_file(args.config)
    if not config:
        logger.error('Unable to get config from file %s', args.config)
        return 1
    logger.info('Config:\n%s', pprint.pformat(config))
    if populate(config, args.threads):
        return 0
    logger.error('Failed to populate from config')
    return 1


if __name__ == '__main__':
    sys.exit(main())
