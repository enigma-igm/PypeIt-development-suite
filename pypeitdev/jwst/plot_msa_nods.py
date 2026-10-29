"""
Plot the calwebb _cal.fits 2D cutouts of one MSA source for several nod exposures, stacked on one set of axes,
with the expected source position from the slit WCS marked on each.

Usage: python3 -m pypeitdev.jwst.plot_msa_nods 4713_6733 /path/to/level2/*_nrs1_cal.fits --outfile /path/to/plots/nods.png
"""
import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
from stdatamodels.jwst import datamodels


def setup_rcparams():
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def load_cutout(cal_file, source_name, ref_wave_um):
    with datamodels.open(cal_file) as model:
        for slit in model.slits:
            if str(slit.source_name).strip() != source_name:
                continue
            data = np.array(slit.data, dtype=float)
            wave = np.array(slit.wavelength, dtype=float)
            # Expected source position in subimage pixels from the WCS (world -> detector)
            _, ysrc = slit.meta.wcs.invert(slit.source_ra, slit.source_dec, ref_wave_um)
            return data, wave, float(ysrc), str(slit.shutter_state)
    raise ValueError('{0} not found in {1}'.format(source_name, cal_file))


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('source_name', type=str, help='Source name, e.g. 4713_6733')
    parser.add_argument('cal_files', type=str, nargs='+', help='calwebb _cal.fits files, one per nod')
    parser.add_argument('--outfile', type=str, required=True, help='Output PNG file')
    parser.add_argument('--lya_wave', type=float, default=9440.0, help='Wavelength to mark [Angstrom]')
    parser.add_argument('--ref_wave', type=float, default=1.1, help='Wavelength for the WCS source position [micron]')
    args = parser.parse_args()

    cutouts = [load_cutout(f, args.source_name, args.ref_wave) for f in args.cal_files]
    nspat, nspec = cutouts[0][0].shape
    gap = 3

    # Stack the nods vertically with a NaN gap; remove a rough sky with the median of each column
    stack, offsets = [], []
    for data, _, _, _ in cutouts:
        offsets.append(sum(s.shape[0] for s in stack))
        sky = np.nanmedian(data, axis=0)
        stack.append(data - sky[None, :])
        stack.append(np.full((gap, nspec), np.nan))
    image = np.vstack(stack[:-1])

    # Approximate wavelength axis: median wavelength of each column in the first exposure
    wave_axis = 1e4 * np.nanmedian(cutouts[0][1], axis=0)
    good = np.isfinite(wave_axis)
    spec_pix = np.arange(nspec)
    finite = image[np.isfinite(image)]
    vmin, vmax = np.percentile(finite, [5, 99])

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(16, 7))
    ax.imshow(image, origin='lower', aspect='auto', cmap='gray_r', vmin=vmin, vmax=vmax,
              extent=[-0.5, nspec - 0.5, -0.5, image.shape[0] - 0.5], interpolation='nearest')

    for (data, _, ysrc, state), off, f in zip(cutouts, offsets, args.cal_files):
        ax.axhline(off + ysrc, color='tab:orange', lw=1.5, ls='--')
        expname = os.path.basename(f).split('_nrs')[0].replace('jw04713001002_', '').replace('_', r'\_')
        ax.text(5, off + 1.5, r'$\rm {0}~~shutters~{1}$'.format(expname, state.replace('x', r'\times')),
                color='tab:blue', fontsize=14)

    lya_pix = np.interp(args.lya_wave, wave_axis[good], spec_pix[good])
    ax.axvline(lya_pix, color='tab:red', lw=1.5, ls=':')
    ax.text(lya_pix + 5, image.shape[0] - 4, r'$\rm Ly\alpha~at~{0:.0f}~\AA$'.format(args.lya_wave),
            color='tab:red', fontsize=14)

    tick_waves = np.arange(7500, 12600, 500)
    ax.set_xticks(np.interp(tick_waves, wave_axis[good], spec_pix[good]))
    ax.set_xticklabels([r'${0:d}$'.format(int(w)) for w in tick_waves])
    ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
    ax.set_ylabel(r'$\rm Spatial~pixel~(nods~stacked)$')
    ax.tick_params(which='major', direction='in', top=True, right=True, length=10, width=1.5)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=5, width=1)
    ax.minorticks_on()

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    fig.savefig(args.outfile, dpi=200, bbox_inches='tight')
    print('Wrote', args.outfile)
    for (data, _, ysrc, state), f in zip(cutouts, args.cal_files):
        print('{0}: shutters {1}, expected source spatial pixel {2:.1f} of {3}'.format(
            os.path.basename(f), state, ysrc, nspat))


if __name__ == '__main__':
    main()
