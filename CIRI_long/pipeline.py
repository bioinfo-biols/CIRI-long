import logging
import subprocess
LOGGER = logging.getLogger('CIRI-long')


def run_ccs(in_file, out_dir, prefix, threads, debugging):
    ccs_fa = '{}/tmp/{}.ccs.fa'.format(out_dir, prefix)
    raw_fa = '{}/tmp/{}.raw.fa'.format(out_dir, prefix)

    ccs_cmd = 'ccs -i {} -o {} -r {} -t {}'.format(
        in_file, ccs_fa, raw_fa, threads,
    )
    with open('{}/{}.log'.format(out_dir, prefix), 'a') as log:
        LOGGER.debug(ccs_cmd)
        ret = subprocess.call(ccs_cmd, shell=True, stdout=log)

    return ret
