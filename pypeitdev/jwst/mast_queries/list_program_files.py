"""
List the NIRSpec files MAST has for one JWST program, without downloading anything.

Usage: python3 -m pypeitdev.jwst.mast_queries.list_program_files 4713
"""
import argparse
from collections import Counter

from astroquery.mast import MastMissionsClass


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('program', type=str, help='JWST program ID, e.g. 4713')
    args = parser.parse_args()

    mast = MastMissionsClass(mission='JWST')
    datasets = mast.query_criteria(program=args.program)
    print(f'{len(datasets)} datasets in program {args.program}')
    cols = [c for c in ['fileSetName', 'targprop', 'exp_type', 'opticalElements',
                        'opt_elem', 'opt_elem2', 'date_obs', 'access'] if c in datasets.colnames]
    datasets[cols].pprint(max_lines=-1, max_width=-1)

    products = mast.get_unique_product_list(datasets)
    suffixes = Counter(fn.rsplit('_', 1)[-1] for fn in products['filename'])
    print('\nProduct counts by file suffix:', dict(suffixes))
    msa = [fn for fn in products['filename'] if fn.endswith('_msa.fits')]
    print('Multi-object setup (msa) files:', msa)


if __name__ == '__main__':
    main()
