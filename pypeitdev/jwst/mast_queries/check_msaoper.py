"""
Report the entries of NIRSpec MSA shutter-operability (msaoper) reference files near given shutters.

Usage: python3 -m pypeitdev.jwst.mast_queries.check_msaoper ~/crds_cache/references/jwst/nirspec/jwst_nirspec_msaoper_*.json --quadrant 4 --row 68 --columns 67 68 69
"""
import argparse
import json
import os


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('files', type=str, nargs='+', help='msaoper JSON files')
    parser.add_argument('--quadrant', type=int, required=True, help='MSA quadrant')
    parser.add_argument('--row', type=int, required=True, help='Shutter row from the _msa.fits file')
    parser.add_argument('--columns', type=int, nargs='+', required=True, help='Shutter columns from the _msa.fits file')
    parser.add_argument('--radius', type=int, default=1, help='Also report entries this many shutters away')
    args = parser.parse_args()

    for f in args.files:
        with open(os.path.expanduser(f)) as fh:
            ref = json.load(fh)
        entries = [e for e in ref['msaoper'] if e['Q'] == args.quadrant]
        print('{0}: useafter {1}; {2} non-operable shutters listed in quadrant {3}'.format(
            os.path.basename(f), ref.get('useafter'), len(entries), args.quadrant))
        cols = range(min(args.columns) - args.radius, max(args.columns) + args.radius + 1)
        rows = range(args.row - args.radius, args.row + args.radius + 1)
        # The mapping of x, y to row, column is checked both ways
        for label, match in [('x=row, y=column', lambda e: e['x'] in rows and e['y'] in cols),
                             ('x=column, y=row', lambda e: e['x'] in cols and e['y'] in rows)]:
            hits = [e for e in entries if match(e)]
            print('  if {0}: {1}'.format(label, hits if hits else 'no entries'))


if __name__ == '__main__':
    main()
