"""
Run jwst_run_redux on one target from jwst_targets.

Usage: python3 -m pypeitdev.jwst.run_target_redux 4713 G140M 4713_6733 --sources 4713_6733
"""
import argparse

from pypeitdev.jwst.jwst_targets import jwst_targets
from pypeitdev.jwst.jwst_run_redux import jwst_run_redux


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('progid', type=str, help='JWST program ID, e.g. 4713')
    parser.add_argument('disperser', type=str, help='Disperser name as used in jwst_targets, e.g. G140M')
    parser.add_argument('target', type=str, help='Target name as used in jwst_targets, e.g. 4713_6733')
    parser.add_argument('--slit', type=str, default=None, help='Slit name for jwst_targets, for fixed-slit targets')
    parser.add_argument('--sources', type=str, nargs='+', default=None, help='MSA source names to reduce; default is all')
    parser.add_argument('--slits', type=str, nargs='+', default=None, help='Slit names to reduce; default is all')
    parser.add_argument('--source_type', type=str, default='EXTENDED', choices=['EXTENDED', 'POINT'],
                        help='Source type for the spec2 pipeline: EXTENDED for MSA, POINT for fixed slit')
    parser.add_argument('--bkg_redux', action='store_true',
                        help='Subtract the other nod positions in the same exposure set instead of fitting the sky')
    parser.add_argument('--bkg_redux_skysub', action='store_true',
                        help='With --bkg_redux, also fit and remove residual sky in the difference image (fixed slits)')
    parser.add_argument('--pypeit_subdir', type=str, default=None,
                        help='Subdirectory of the reduction directory for PypeIt output; '
                             'default is pypeit, or pypeit_bkg_redux with --bkg_redux')
    parser.add_argument('--redux_suffix', type=str, default='',
                        help='Appended to the reduction directory from jwst_targets, e.g. _point, so a run with '
                             'different calwebb settings does not overwrite the existing stage 2 output')
    parser.add_argument('--snr_thresh', type=float, default=10.0, help='Object-finding signal-to-noise threshold')
    parser.add_argument('--boxcar_radius', type=float, default=None,
                        help='Boxcar extraction radius [arcsec]; 0.1 arcsec per pixel, so 0.65 gives a 13-pixel box')
    parser.add_argument('--show', action='store_true', help='Show images and QA in ginga')
    parser.add_argument('--overwrite_stage2', action='store_true', help='Rerun the spec2 pipeline even if cal files exist')
    args = parser.parse_args()

    exp_list, redux_dir, _ = jwst_targets(args.progid, args.disperser, args.target, slit=args.slit)
    redux_dir = redux_dir.rstrip('/') + args.redux_suffix
    is_rate = all(f.endswith('_rate.fits') for det_files in exp_list for f in det_files)
    file_kwargs = {'rate_list': exp_list} if is_rate else {'uncal_list': exp_list}
    pypeit_subdir = args.pypeit_subdir
    if pypeit_subdir is None:
        pypeit_subdir = 'pypeit'
        if args.bkg_redux:
            pypeit_subdir = 'pypeit_bkg_redux_skysub' if args.bkg_redux_skysub else 'pypeit_bkg_redux'

    jwst_run_redux(redux_dir, args.source_type, reduce_slits=args.slits, reduce_sources=args.sources,
                   show=args.show, overwrite_stage2=args.overwrite_stage2, snr_thresh=args.snr_thresh,
                   bkg_redux=args.bkg_redux, pypeit_subdir=pypeit_subdir, boxcar_radius=args.boxcar_radius,
                   bkg_redux_skysub=args.bkg_redux_skysub,
                   **file_kwargs)


if __name__ == '__main__':
    main()
