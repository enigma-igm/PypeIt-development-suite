"""
Overlay several combined 1D spectra (outputs of combine_spec1d, in microjansky) on one set of axes.

Usage: python3 -m pypeitdev.jwst.plot_overlay_spec1d a.fits b.fits --labels "visit 002" "observation 002" --outfile overlay.png
"""
import argparse
import os

from matplotlib import pyplot as plt

from pypeitdev.jwst.plot_combined_spec1d import load_combined, setup_rcparams, LYA_REST

COLORS = ['black', 'tab:blue', 'tab:orange', 'tab:green']


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().split('\n')[0])
    parser.add_argument('combined_files', type=str, nargs='+', help='Outputs of combine_spec1d')
    parser.add_argument('--labels', type=str, nargs='+', required=True, help='Legend label for each file')
    parser.add_argument('--z_qso', type=float, default=6.327, help='Quasar redshift')
    parser.add_argument('--z_gal', type=float, default=6.77, help='Galaxy redshift')
    parser.add_argument('--wave_min', type=float, default=8000.0, help='Minimum wavelength [Angstrom]')
    parser.add_argument('--wave_max', type=float, default=13000.0, help='Maximum wavelength [Angstrom]')
    parser.add_argument('--ymin', type=float, default=-0.35, help='Minimum flux to plot [microjansky]')
    parser.add_argument('--ymax', type=float, default=1.1, help='Maximum flux to plot [microjansky]')
    parser.add_argument('--outfile', type=str, required=True, help='Output PNG file')
    args = parser.parse_args()

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(16, 6))
    for f, label, color in zip(args.combined_files, args.labels, COLORS):
        wave, flux, sigma = load_combined(f)
        ax.plot(wave, flux, color=color, lw=1.2, drawstyle='steps-mid', alpha=0.85,
                label=r'$\rm {0}$'.format(label.replace(' ', '~')))
        ax.plot(wave, sigma, color=color, lw=1.0, ls=':', drawstyle='steps-mid', alpha=0.7)
    ax.axhline(0.0, color='0.4', ls='--', lw=1.5)

    lya_qso = LYA_REST * (1 + args.z_qso)
    lya_gal = LYA_REST * (1 + args.z_gal)
    ax.axvline(lya_qso, color='tab:red', ls='--', lw=2)
    ax.axvline(lya_gal, color='orange', ls='--', lw=2)
    ax.text(lya_qso - 20, 0.97 * args.ymax, r'$z_{\rm QSO}$', color='tab:red', rotation=90, ha='right', va='top')
    ax.text(lya_gal + 20, 0.97 * args.ymax, r'$z_{\rm gal}$', color='orange', rotation=90, ha='left', va='top')

    ax.set_xlim(args.wave_min, args.wave_max)
    ax.set_ylim(args.ymin, args.ymax)
    ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
    ax.set_ylabel(r'$f_\nu~[10^{-29}~\rm erg~s^{-1}~cm^{-2}~Hz^{-1}]$')
    ax.legend(fontsize=14, loc='upper right', frameon=False, title=r'$\rm dotted:~1\sigma~error$', title_fontsize=14)
    ax.tick_params(which='major', direction='in', top=True, right=True, length=10, width=1.5)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=5, width=1)
    ax.minorticks_on()

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    fig.savefig(args.outfile, dpi=200, bbox_inches='tight')
    print('Wrote', args.outfile)


if __name__ == '__main__':
    main()
