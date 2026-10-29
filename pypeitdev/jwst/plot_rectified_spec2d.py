"""
Plot the inverse-variance weighted mean of background-subtracted spec2d images from several exposures on the same
pixel grid, rectified along the object trace from the matching spec1d files.

Usage: python3 -m pypeitdev.jwst.plot_rectified_spec2d --spec2d a.fits b.fits --spec1d a1.fits b1.fits --outfile spec2d.png
"""
import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage
from pypeit import spec2dobj, specobjs

LYA_REST = 1215.67


def setup_rcparams():
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def rectify(image, trace_spat, half_width):
    """Sample image (nspec, nspat) at trace_spat + offsets for each spectral row; linear interpolation."""
    nspec = image.shape[0]
    offsets = np.arange(-half_width, half_width + 1)
    spec_idx = np.repeat(np.arange(nspec)[:, None], offsets.size, axis=1)
    spat_idx = trace_spat[:, None] + offsets[None, :]
    return ndimage.map_coordinates(image, [spec_idx, spat_idx], order=1, mode='constant', cval=np.nan), offsets


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('--spec2d', type=str, nargs='+', required=True, help='spec2d files')
    parser.add_argument('--spec1d', type=str, nargs='+', required=True, help='spec1d files, same order as --spec2d')
    parser.add_argument('--half_width', type=int, default=10, help='Spatial pixels to show either side of the trace')
    parser.add_argument('--z_qso', type=float, default=6.327, help='Quasar redshift')
    parser.add_argument('--z_gal', type=float, default=6.77, help='Galaxy redshift')
    parser.add_argument('--wave_min', type=float, default=8000.0, help='Minimum wavelength [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=13000.0, help='Maximum wavelength [Angstrom]')
    parser.add_argument('--outfile', type=str, required=True, help='Output PNG file')
    args = parser.parse_args()

    wave_trace, offsets, mean, box_radius = rectified_mean(args.spec2d, args.spec1d, args.half_width)

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(16, 3.2))
    draw_spec2d(ax, wave_trace, offsets, mean, box_radius, args.z_qso, args.z_gal, args.wave_min, args.wave_max)
    ax.set_xlabel(r'$\rm Wavelength~[\AA]$')

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    fig.savefig(args.outfile, dpi=200, bbox_inches='tight')
    print('Wrote', args.outfile)


def regrid(wave_new, wave_old, rect, rect_ivar):
    """Linearly interpolate a rectified image (nspec, noffset) and its ivar from wave_old onto wave_new."""
    good = np.isfinite(wave_old) & (wave_old > 1.0)
    out = np.full((wave_new.size, rect.shape[1]), np.nan)
    out_ivar = np.zeros((wave_new.size, rect.shape[1]))
    for j in range(rect.shape[1]):
        ok = good & np.isfinite(rect[:, j]) & (rect_ivar[:, j] > 0)
        if ok.sum() < 2:
            continue
        out[:, j] = np.interp(wave_new, wave_old[ok], rect[ok, j], left=np.nan, right=np.nan)
        out_ivar[:, j] = np.interp(wave_new, wave_old[ok], rect_ivar[ok, j], left=0.0, right=0.0)
    out_ivar = np.where(np.isfinite(out), out_ivar, 0.0)
    return out, out_ivar


def rectified_mean(spec2d_files, spec1d_files, half_width):
    """Inverse-variance weighted mean of the trace-rectified, background-subtracted spec2d images.

    Each rectified image is linearly interpolated in wavelength onto the wavelength grid along the trace of the
    first file, so exposures from different slits (different wavelength grids) can be combined; for display only.

    Returns the wavelength along the trace, the spatial offsets, the mean image (nspec, noffset),
    and the boxcar radius in pixels.
    """
    num, den, wave_trace, box_radius = 0.0, 0.0, None, None
    for f2d, f1d in zip(spec2d_files, spec1d_files):
        s = spec2dobj.Spec2DObj.from_file(f2d, 'DET01')
        sobj = specobjs.SpecObjs.from_fitsfile(f1d, chk_version=False)[0]
        trace = np.asarray(sobj.TRACE_SPAT, dtype=float)
        gpm = (s.bpmmask.mask == 0) & (s.ivarraw > 0)
        image = np.where(gpm, s.sciimg - s.skymodel, np.nan)
        ivar = np.where(gpm, s.ivarraw, 0.0)
        rect, offsets = rectify(image, trace, half_width)
        rect_ivar, _ = rectify(ivar, trace, half_width)
        rect_ivar = np.where(np.isfinite(rect) & np.isfinite(rect_ivar), rect_ivar, 0.0)
        wave_this = ndimage.map_coordinates(np.asarray(s.waveimg, dtype=float),
                                            [np.arange(trace.size), trace], order=1)
        if wave_trace is None:
            wave_trace = wave_this
            # The boxcar radius is not stored; half the number of pixels summed gives it
            box_radius = 0.5 * np.median(sobj.BOX_NPIX[sobj.BOX_NPIX > 0])
        else:
            rect, rect_ivar = regrid(wave_trace, wave_this, rect, rect_ivar)
        num = num + np.nan_to_num(rect) * rect_ivar
        den = den + rect_ivar
    mean = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return wave_trace, offsets, mean, box_radius


def draw_spec2d(ax, wave_trace, offsets, mean, box_radius, z_qso, z_gal, wave_min, wave_max):
    """Draw the rectified 2D spectrum, the boxcar window, and the Lyman-alpha lines on ax (no x label)."""
    good = np.isfinite(wave_trace) & (wave_trace > 1.0)
    use = good & (wave_trace > wave_min) & (wave_trace < wave_max)
    shown = mean[use]
    vmin, vmax = np.nanpercentile(shown, [2, 98])
    # Rows of the rectified image are spectral pixels: plot wavelength along x, offset from trace along y
    ax.pcolormesh(wave_trace[use], offsets, shown.T, cmap='gray_r', vmin=vmin, vmax=vmax, shading='nearest')
    for sign in [-1, 1]:
        ax.axhline(sign * box_radius, color='tab:blue', ls='-', lw=1.5,
                   label=r'$\rm boxcar~extraction~window~(\pm{0:g}~pix)$'.format(box_radius) if sign > 0 else None)
    ax.legend(fontsize=14, loc='upper right', frameon=True, framealpha=0.9)
    ax.axvline(LYA_REST * (1 + z_qso), color='tab:red', ls='--', lw=2)
    ax.axvline(LYA_REST * (1 + z_gal), color='orange', ls='--', lw=2)
    ax.set_xlim(wave_min, wave_max)
    ax.set_ylabel(r'$\rm Offset~[pix]$')
    ax.tick_params(which='major', direction='in', top=True, right=True, length=8, width=1.5)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=4, width=1)
    ax.minorticks_on()


if __name__ == '__main__':
    main()
