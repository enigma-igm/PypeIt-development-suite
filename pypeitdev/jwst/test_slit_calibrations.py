"""
Check that NIRSpecSlitCalibrations finds a slit by name in the cal and interpolatedflat models separately.

Usage: python3 -m pypeitdev.jwst.test_slit_calibrations
"""
import os

from stdatamodels.jwst import datamodels
from pypeit.spectrographs.util import load_spectrograph

from pypeitdev.jwst.jwst_utils import NIRSpecSlitCalibrations

L2 = os.path.expanduser('~/jwst_redux/redux/NIRSPEC_MSA/4713/G140M/4713_6733_obs2/output/level2')


def check(det, slit_name):
    base = os.path.join(L2, 'jw04713002001_03101_00002_{0}'.format(det))
    spectrograph = load_spectrograph('jwst_nirspec')
    det_par = spectrograph.get_detector_par(1 if det == 'nrs1' else 2)
    with datamodels.open(base + '_cal.fits') as cal, datamodels.open(base + '_interpolatedflat.fits') as flat:
        calib = NIRSpecSlitCalibrations(det_par, cal, flat, slit_name)
        print('{0} {1}: on_detector {2}, cal index {3}, source {4}, waveimg shape {5}'.format(
            det, slit_name, calib.on_detector, calib.slit_indx, calib.source_name,
            None if not calib.on_detector else calib.waveimg.shape))
        return calib


def main():
    c = check('nrs1', 'S200A1')
    assert c.on_detector and c.source_name.strip() == '4713_6733'
    c = check('nrs2', '5')
    assert c.on_detector
    print('OK')


if __name__ == '__main__':
    main()
