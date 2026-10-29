"""
Plot a combined 1D spectrum from combine_spec1d (flux density in microjansky) in the style of the
Eilers et al. 2025 galaxy spectra, marking the quasar and galaxy Lyman-alpha wavelengths.

Usage: python3 -m pypeitdev.jwst.plot_combined_spec1d combined.fits --z_qso 6.327 --z_gal 6.77 --outfile spec.png
"""
import argparse
import os

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

LYA_REST = 1215.67
# Interstellar absorption lines [Angstrom, vacuum rest frame]; doublets are marked at both components
ABS_LINES = [(r'Si\,II', [1260.42]), (r'Si\,IV', [1393.76, 1402.77]), (r'C\,IV', [1548.20, 1550.78])]


def setup_rcparams():
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('combined_file', type=str, help='Output of combine_spec1d')
    parser.add_argument('--z_qso', type=float, default=6.327, help='Quasar redshift')
    parser.add_argument('--z_gal', type=float, default=6.77, help='Galaxy redshift')
    parser.add_argument('--wave_min', type=float, default=8000.0, help='Minimum wavelength to plot [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=13000.0, help='Maximum wavelength to plot [Angstrom]')
    parser.add_argument('--ymin', type=float, default=-0.4, help='Minimum flux to plot [microjansky]')
    parser.add_argument('--ymax', type=float, default=2.1, help='Maximum flux to plot [microjansky]')
    parser.add_argument('--label', type=str, default='PypeIt nod subtracted, 2 exposures', help='Legend label')
    parser.add_argument('--label_lines', action='store_true', help='Mark Si II, Si IV and C IV absorption at z_gal')
    parser.add_argument('--outfile', type=str, required=True, help='Output PNG file')
    args = parser.parse_args()

    wave, flux, sigma = load_combined(args.combined_file)

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(16, 6))
    draw_spec1d(ax, wave, flux, sigma, args.z_qso, args.z_gal, args.wave_min, args.wave_max, args.ymin, args.ymax,
                label=args.label, label_lines=args.label_lines)
    ax.set_xlabel(r'$\rm Wavelength~[\AA]$')

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    fig.savefig(args.outfile, dpi=200, bbox_inches='tight')
    print('Wrote', args.outfile)


def load_combined(combined_file):
    """Return wavelength, flux and 1 sigma error of the good pixels of a combine_spec1d output."""
    t = Table.read(combined_file)
    good = np.asarray(t['mask'], dtype=bool)
    wave = np.asarray(t['wave'])[good]
    flux = np.asarray(t['flux'])[good]
    ivar = np.asarray(t['ivar'])[good]
    sigma = np.where(ivar > 0, 1.0 / np.sqrt(np.where(ivar > 0, ivar, 1.0)), np.nan)
    return wave, flux, sigma


def draw_absorption_lines(ax, z_gal, y_tick, ymax):
    """Mark interstellar absorption lines at the galaxy redshift with short ticks and labels above y_tick."""
    tick_len = 0.04 * ymax
    for name, rest_waves in ABS_LINES:
        obs = [w * (1 + z_gal) for w in rest_waves]
        for w in obs:
            ax.plot([w, w], [y_tick, y_tick + tick_len], color='tab:green', lw=1.5)
        ax.text(np.mean(obs), y_tick + 1.4 * tick_len, r'$\rm {0}$'.format(name), color='tab:green',
                ha='center', va='bottom', fontsize=14)


def draw_spec1d(ax, wave, flux, sigma, z_qso, z_gal, wave_min, wave_max, ymin, ymax,
                label='PypeIt nod subtracted, 2 exposures', label_lines=False):
    """Draw the spectrum, its error, the Lyman-alpha lines and optionally absorption lines on ax (no x label)."""
    if label_lines:
        draw_absorption_lines(ax, z_gal, 0.36 * ymax, ymax)
    ax.plot(wave, flux, color='black', lw=1.5, drawstyle='steps-mid', label=r'$\rm {0}$'.format(label.replace(' ', '~')))
    ax.plot(wave, sigma, color='0.6', lw=1.5, drawstyle='steps-mid', label=r'$\rm 1\sigma~error$')
    ax.axhline(0.0, color='0.4', ls='--', lw=1.5)

    lya_qso = LYA_REST * (1 + z_qso)
    lya_gal = LYA_REST * (1 + z_gal)
    ax.axvline(lya_qso, color='tab:red', ls='--', lw=2)
    ax.axvline(lya_gal, color='orange', ls='--', lw=2)
    ax.text(lya_qso - 20, 0.97 * ymax, r'$z_{\rm QSO}$', color='tab:red', rotation=90, ha='right', va='top')
    ax.text(lya_gal + 20, 0.97 * ymax, r'$z_{\rm gal}$', color='orange', rotation=90, ha='left', va='top')

    ax.set_xlim(wave_min, wave_max)
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(r'$f_\nu~[10^{-29}~\rm erg~s^{-1}~cm^{-2}~Hz^{-1}]$')
    ax.legend(fontsize=14, loc='upper right', frameon=False)
    ax.tick_params(which='major', direction='in', top=True, right=True, length=10, width=1.5)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=5, width=1)
    ax.minorticks_on()


if __name__ == '__main__':
    main()
