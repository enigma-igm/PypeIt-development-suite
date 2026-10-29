"""
Download from MAST the exposure files that jwst_targets lists for one target, skipping files already on disk.

Usage: python3 -m pypeitdev.jwst.mast_queries.download_target_files 4713 G140M 4713_6733
"""
import argparse
import os

import numpy as np
from astroquery.mast import MastMissionsClass

from pypeitdev.jwst.jwst_targets import jwst_targets


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('progid', type=str, help='JWST program ID, e.g. 4713')
    parser.add_argument('disperser', type=str, help='Disperser name as used in jwst_targets, e.g. G140M')
    parser.add_argument('target', type=str, help='Target name as used in jwst_targets, e.g. 4713_6733')
    parser.add_argument('--slit', type=str, default=None, help='Slit name, for fixed-slit targets')
    parser.add_argument('--dry_run', action='store_true', help='Only list the files that would be downloaded')
    args = parser.parse_args()

    exp_list, _, rawpath = jwst_targets(args.progid, args.disperser, args.target, slit=args.slit)
    wanted = [f for det_files in exp_list for f in det_files]
    missing = [os.path.basename(f) for f in wanted if not os.path.isfile(f)]
    print(f'{len(wanted)} files listed, {len(missing)} missing from {rawpath}')
    if len(missing) == 0:
        return

    mast = MastMissionsClass(mission='JWST')
    products = mast.get_unique_product_list(mast.query_criteria(program=args.progid))
    products = products[np.isin(products['filename'], missing)]
    not_on_mast = sorted(set(missing) - set(products['filename']))
    if len(not_on_mast) > 0:
        print('Not found on MAST:', not_on_mast)
    for fn in products['filename']:
        print('  ' + fn)
    if args.dry_run:
        return

    os.makedirs(rawpath, exist_ok=True)
    mast.download_products(products, download_dir=rawpath, flat=True)


if __name__ == '__main__':
    main()
