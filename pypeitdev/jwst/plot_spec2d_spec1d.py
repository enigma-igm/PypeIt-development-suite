"""
Plot the rectified 2D spectrum (top) above the combined 1D spectrum (bottom) on a shared wavelength axis.

Usage: python3 -m pypeitdev.jwst.plot_spec2d_spec1d combined.fits --spec2d a.fits b.fits --spec1d a1.fits b1.fits --outfile spec.png
"""
import argparse
import os

from matplotlib import pyplot as plt

from pypeitdev.jwst.plot_combined_spec1d import load_combined, draw_spec1d, setup_rcparams
from pypeitdev.jwst.plot_rectified_spec2d import rectified_mean, draw_spec2d


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('combined_file', type=str, help='Output of combine_spec1d')
    parser.add_argument('--spec2d', type=str, nargs='+', required=True, help='spec2d files')
    parser.add_argument('--spec1d', type=str, nargs='+', required=True, help='spec1d files, same order as --spec2d')
    parser.add_argument('--half_width', type=int, default=10, help='Spatial pixels to show either side of the trace')
    parser.add_argument('--z_qso', type=float, default=6.327, help='Quasar redshift')
    parser.add_argument('--z_gal', type=float, default=6.77, help='Galaxy redshift')
    parser.add_argument('--wave_min', type=float, default=8000.0, help='Minimum wavelength [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=13000.0, help='Maximum wavelength [Angstrom]')
    parser.add_argument('--ymin', type=float, default=-0.35, help='Minimum flux to plot [microjansky]')
    parser.add_argument('--ymax', type=float, default=1.1, help='Maximum flux to plot [microjansky]')
    parser.add_argument('--label', type=str, default='PypeIt nod subtracted, 2 exposures', help='Legend label for the 1D spectrum')
    parser.add_argument('--label_lines', action='store_true', help='Mark Si II, Si IV and C IV absorption at z_gal')
    parser.add_argument('--outfile', type=str, required=True, help='Output PNG file')
    args = parser.parse_args()

    wave_trace, offsets, mean, box_radius = rectified_mean(args.spec2d, args.spec1d, args.half_width)
    wave, flux, sigma = load_combined(args.combined_file)

    setup_rcparams()
    fig, (ax2d, ax1d) = plt.subplots(2, 1, figsize=(16, 8), sharex=True,
                                     gridspec_kw={'height_ratios': [1, 4.5], 'hspace': 0.06})
    draw_spec2d(ax2d, wave_trace, offsets, mean, box_radius, args.z_qso, args.z_gal, args.wave_min, args.wave_max)
    # The 2D panel is short: fewer ticks and a smaller label so the label stays within the panel
    ax2d.set_yticks([-5, 0, 5])
    ax2d.set_ylabel(r'$\rm Offset~[pix]$', fontsize=15)
    ax2d.get_legend().remove()
    draw_spec1d(ax1d, wave, flux, sigma, args.z_qso, args.z_gal, args.wave_min, args.wave_max, args.ymin, args.ymax,
                label=args.label, label_lines=args.label_lines)
    ax1d.set_xlabel(r'$\rm Wavelength~[\AA]$')

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    fig.savefig(args.outfile, dpi=200, bbox_inches='tight')
    print('Wrote', args.outfile)


if __name__ == '__main__':
    main()
