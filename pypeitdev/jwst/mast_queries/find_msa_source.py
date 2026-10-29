"""
Download a NIRSpec multi-object setup file (_msa.fits) and report the shutters assigned to one source.

Usage: python3 -m pypeitdev.jwst.mast_queries.find_msa_source jw04713001002_01_msa.fits 6733 --outdir ~/jwst_redux/Raw/NIRSPEC_MSA/4713
"""
import argparse
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import MastMissionsClass


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('msa_file', type=str, help='Name of the _msa.fits file on MAST')
    parser.add_argument('source_id', type=int, help='Source ID in the multi-object setup catalog')
    parser.add_argument('--outdir', type=str, default='.', help='Directory to save the _msa.fits file in')
    args = parser.parse_args()

    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    local = os.path.join(outdir, args.msa_file)
    if not os.path.isfile(local):
        mast = MastMissionsClass(mission='JWST')
        program = args.msa_file[3:7]
        products = mast.get_unique_product_list(mast.query_criteria(program=program))
        products = products[products['filename'] == args.msa_file]
        mast.download_products(products, download_dir=outdir, flat=True)

    with fits.open(local) as hdul:
        shutters = Table(hdul['SHUTTER_INFO'].data)
        sources = Table(hdul['SOURCE_INFO'].data)

    src = sources[sources['source_id'] == args.source_id]
    print('Source catalog entry:')
    src.pprint(max_width=-1)

    rows = shutters[shutters['source_id'] == args.source_id]
    print(f'\nShutters assigned to source {args.source_id}:')
    rows.pprint(max_lines=-1, max_width=-1)
    print('\nmsa_metadata_id values containing the source:', np.unique(rows['msa_metadata_id']).tolist())
    print('slitlet_id values:', np.unique(rows['slitlet_id']).tolist())


if __name__ == '__main__':
    main()
