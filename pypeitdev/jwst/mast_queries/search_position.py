"""
List JWST datasets on MAST around a sky position, grouped by program, instrument mode and optical elements.

Usage: python3 -m pypeitdev.jwst.mast_queries.search_position 15.0659247 28.0350864 --radius 3
"""
import argparse
from collections import Counter

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import MastMissionsClass


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('ra', type=float, help='Right ascension [deg]')
    parser.add_argument('dec', type=float, help='Declination [deg]')
    parser.add_argument('--radius', type=float, default=3.0, help='Search radius [arcmin]')
    args = parser.parse_args()

    mast = MastMissionsClass(mission='JWST')
    coord = SkyCoord(args.ra, args.dec, unit='deg')
    rows = mast.query_region(coord, radius=args.radius * u.arcmin, limit=5000)
    print('{0} datasets within {1} arcmin'.format(len(rows), args.radius))
    keys = Counter((str(r['program']), str(r['exp_type']), str(r['opticalElements']))
                   for r in rows if 'program' in rows.colnames)
    for (prog, exptype, optics), n in sorted(keys.items()):
        print('  program {0:>5}  {1:<15} {2:<20} {3} datasets'.format(prog, exptype, optics, n))


if __name__ == '__main__':
    main()
