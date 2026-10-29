"""
Rerun PypeIt object finding on the sky-subtracted image stored in a spec2d file, to see why no object was found.

Usage: python3 -m pypeitdev.jwst.debug_objfind /path/to/spec2d_file.fits --snr_thresh 10 --fwhm 2 --trim 5
"""
import argparse

import numpy as np
from pypeit import spec2dobj
from pypeit.core import findobj_skymask


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('spec2d_file', type=str, help='spec2d file')
    parser.add_argument('--snr_thresh', type=float, default=10.0, help='Object-finding signal-to-noise threshold')
    parser.add_argument('--fwhm', type=float, default=2.0, help='Assumed object full width at half maximum [pixels]')
    parser.add_argument('--trim', type=int, default=5, help='Pixels trimmed from each slit edge before finding')
    parser.add_argument('--find_window', type=str, default='slit', choices=['omit', 'slit', 'none_none'],
                        help='find_min_max passed to objs_in_slit: omit it, use the slit specmin/specmax '
                             '(as FindObjects does by default), or [None, None]')
    parser.add_argument('--trace_window', type=str, default='slit', choices=['omit', 'slit', 'none_none'],
                        help='spec_min_max (trace window) passed to objs_in_slit, same choices as --find_window. '
                             'FindObjects always passes one, so omit does not match a real run')
    args = parser.parse_args()

    s = spec2dobj.Spec2DObj.from_file(args.spec2d_file, 'DET01')
    left, right, _ = s.slits.select_edges()
    slitid_img = s.slits.slit_img()
    thismask = slitid_img == s.slits.spat_id[0]
    inmask = (s.bpmmask.mask == 0) & thismask
    image = s.sciimg - s.skymodel

    windows = {'omit': None, 'slit': [s.slits.specmin[0], s.slits.specmax[0]], 'none_none': [None, None]}
    find_min_max = windows[args.find_window]
    spec_min_max = windows[args.trace_window]

    print('pixels in slit: {0}, good pixels: {1}, find_min_max: {2}, spec_min_max: {3}'.format(
        thismask.sum(), inmask.sum(), find_min_max, spec_min_max))
    sobjs = findobj_skymask.objs_in_slit(image, s.ivarraw, thismask, left[:, 0], right[:, 0], inmask=inmask,
                                         fwhm=args.fwhm, snr_thresh=args.snr_thresh,
                                         trim_edg=(args.trim, args.trim), find_min_max=find_min_max,
                                         spec_min_max=spec_min_max)
    print('objects found: {0}'.format(sobjs.nobj))
    for sobj in sobjs:
        print('  spatial pixel {0:.1f}, signal-to-noise {1:.1f}'.format(
            np.median(sobj.TRACE_SPAT), sobj.smash_snr if sobj.smash_snr is not None else np.nan))


if __name__ == '__main__':
    main()
