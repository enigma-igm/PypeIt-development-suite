""" For the development and testing of TraceSlits
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from astropy.io import fits
import glob

from pypit import msgs
from pypit import ardebug as debugger
from pypit import ginga
from pypit import arload
from pypit import arproc
from pypit import arcomb
from pypit import ardeimos
from pypit import arlris
from pypit import arpixels
from pypit import arsave
from pypit import traceslits

debug = debugger.init()
debug['develop'] = True
msgs.reset(debug=debug, verbosity=2)

from pypit import artrace

def_settings = traceslits.default_settings
#def_settings=dict(trace={'slits': {'single': [],
#                                   'function': 'legendre',
#                                   'polyorder': 3,
#                                   'diffpolyorder': 2,
#                                   'fracignore': 0.01,
#                                   'medrep': 0,
#                                   'pad': 0,
#                                   'number': -1,
#                                   'maxgap': None,
#                                   'sigdetect': 20.,
#                                   'pca': {'params': [3,2,1,0,0,0], 'type': 'pixel', 'extrapolate': {'pos': 0, 'neg':0}},
#                                   'sobel': {'mode': 'nearest'}},
def_settings['trace']['combine'] = {'match': -1.,
                          'satpix': 'reject',
                          'method': 'weightmean',
                          'reject': {'cosmics': 20., 'lowhigh': [0,0], 'level': [3.,3.],
                                     'replace': 'maxnonsat'}}

def pypit_trace_slits():
    pass

def combine_frames(spectrograph, files, det, settings, saturation=None, numamplifiers=None):
    # Grab data info
    datasec, oscansec, naxis0, naxis1 = arproc.get_datasec(spectrograph, files[0],
                                                           numamplifiers=numamplifiers, det=det)

    # Load + process images
    frames = []
    for ifile in files:
        rawframe, head0 = arload.load_raw_frame(spectrograph, ifile, pargs.det, disp_dir=0)
        # Bias subtract
        newframe = arproc.sub_overscan(rawframe.copy(), numamplifiers, datasec, oscansec)
        # Trim
        frame = arproc.trim(newframe, numamplifiers, datasec)
        # Append
        frames.append(frame)

    # Convert to array
    frames_arr = np.zeros((frames[0].shape[0], frames[0].shape[1], len(frames)))
    for ii in range(len(frames)):
        frames_arr[:,:,ii] = frames[ii]

    # Combine
    mstrace = arcomb.core_comb_frames(frames_arr, frametype='Unknown',
                                 method=settings['trace']['combine']['method'],
                                 reject=settings['trace']['combine']['reject'],
                                 satpix=settings['trace']['combine']['satpix'],
                                 saturation=saturation)
    # Return
    return mstrace

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Developing/testing/checking trace_slits [v1.1]')
    parser.add_argument("spectrograph", type=str, help="Instrument [keck_deimos, keck_lris_red, keck_lris_blue]")
    parser.add_argument("--files", type=str, help="File(s) of trace flats")
    parser.add_argument("--det", default=1, type=int, help="Detector")
    parser.add_argument("--show", default=False, action="store_true", help="Show the image with traces")
    parser.add_argument("--outfile", type=str, help="Output to a MasterFrame formatted FITS file")
    parser.add_argument("--pclass", default=False, action="store_true", help="Use the ProcessImages class")
    #parser.add_argument("--driver", default=False, action="store_true", help="Show the image with traces")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):

    xgap = 0.0   # Gap between the square detector pixels (expressed as a fraction of the x pixel size)
    ygap = 0.0   # Gap between the square detector pixels (expressed as a fraction of the x pixel size)
    ysize = 1.0  # The size of a pixel in the y-direction as a multiple of the x pixel size
    binbpx = None
    numamplifiers = None
    saturation = None
    add_user_slits = None
    user_settings = None

    settings = def_settings.copy()

    # Read files
    if pargs.files is not None:
        files = glob.glob(pargs.files+'*')
    else:
        files = None

    # Instrument specific
    if pargs.spectrograph == 'keck_deimos':

        saturation = 65535.0              # The detector Saturation level
        numamplifiers=1

        if files is None:
            #files = glob.glob('../RAW_DATA/Keck_DEIMOS/830G_M/DE.20100913.57*')  # Mask (57006, 57161)
            #files = glob.glob('data/DEIMOS/DE.20100913.57*')  # Mask (57006, 57161)
            #  The following are with sigdetect=20;  sigdetect=50 gets rid of the junk (I think)
            #      det=1 :: 25 slits including star boxes
            #      det=2 :: 26 slits including stars + short first one
            #      det=3 :: 27 slits with a fake, short slit at the start
            #      det=4 :: 27 slits with several junk slits in saturated star boxes
            #      det=5 :: 25 slits with one junk slit on a saturated box
            #      det=6 :: 26 slits with a short, legit first slit
            #      det=7 :: 26 slits with stars
            #      det=8 :: 25 slits well done slits including a short first one
            #
            #files = ['../RAW_DATA/Keck_DEIMOS/830G_L/'+ifile for ifile in [  # Longslit in dets 3,7
            #    'd0914_0014.fits', 'd0914_0015.fits']]
            # Fred's latest
            #files = glob.glob('data/DEIMOS/Trace_flats/d0526_0*')
            # Longslit
            files = glob.glob('data/DEIMOS/Trace_flats/d0423_0*')

            if len(files) == 0:
                print("No files!!")
                debugger.set_trace()


        # Bad pixel mask (important!!)
        binbpx = ardeimos.bpm(pargs.det)

        #hdul = fits.open('trace_slit.fits')
        settings['trace']['slits']['sigdetect'] = 50.0
        settings['trace']['slits']['fracignore'] = 0.02  # 0.02 removes star boxes
        settings['trace']['slits']['pca']['params'] = [3,2,1,0]
        # For combine
        user_settings = dict(run={'spectrograph': 'keck_deimos'})
    elif pargs.spectrograph == 'keck_lris_red':
        saturation = 65535.0              # The detector Saturation level
        numamplifiers=2

        if files is None:
            #python dev_trace_slits.py keck_lris_red --outfile=../Cooked/Trace/MasterTrace_KeckLRISr_150420_402 --det=1
            #files = glob.glob('data/LRIS/Trace_flats/r150420_402*')
            #add_user_slits = [[489,563,1024]] # Goes with r150420_402*  ; and it works
            #    det1 : Missing a slit between two standard stars
            #    det2 : 12 solid slits

            #python dev_trace_slits.py keck_lris_red --outfile=../Cooked/Trace/MasterTrace_KeckLRISr_20160110_A --det=1
            files = ['data/LRIS/Trace_flats/LR.20160110.10103.fits.gz',  # det=1: finds a ghost slit;  crazy edge case..
                     'data/LRIS/Trace_flats/LR.20160110.10273.fits.gz']  # det=2: solid

            #files = ['data/LRIS/Trace_flats/LR.20160110.10644.fits.gz',  # det=1:  Well done! including an overlapping slit
            #         'data/LRIS/Trace_flats/LR.20160110.10717.fits.gz']  # det=2: 21 solid slits including stars

        # Read
        head0 = fits.open(files[0])[0].header
        xbin, ybin = [int(ii) for ii in head0['BINNING'].split(',')]
        binbpx = arlris.core_bpm(xbin, ybin, 'red', pargs.det)

        settings['trace']['slits']['sigdetect'] = 50.0
        settings['trace']['slits']['pca']['params'] = [3,2,1,0]
    elif pargs.spectrograph == 'keck_lris_blue':
        saturation = 65535.0              # The detector Saturation level
        numamplifiers=2
        if files is None:
            files = glob.glob('../RAW_DATA/Keck_LRIS_blue/long_600_4000_d560/b150910_2051*') # Single Twilight
            #files = glob.glob('data/LRIS/Trace_flats/LB.20160109.*')  # det=1 : solid; det=2 solid [sigdetect=30]
            #files = glob.glob('data/LRIS/Trace_flats/LB.20160406.*')  # det=1 : solid;

        settings['trace']['slits']['pca']['params'] = [3,2,1,0]
        settings['trace']['slits']['sigdetect'] = 30.0
    else:
        debugger.set_trace()

    # Combine
    if pargs.pclass:
        from pypit import processimages
        tflats = processimages.ProcessImages(files, user_settings=user_settings)
        mstrace = tflats.combine(bias_subtract='overscan', trim=True)
    else:
        mstrace = combine_frames(pargs.spectrograph, files, pargs.det, settings,
                             saturation=saturation, numamplifiers=numamplifiers)

    # binpx
    if binbpx is None:
        binbpx = np.zeros_like(mstrace)

    # pixlocn
    pixlocn = arpixels.core_gen_pixloc(mstrace)

    # Trace
    tslits = traceslits.TraceSlits(mstrace, pixlocn, binbpx=binbpx, settings=settings)
    tslit_dict = tslits.run(armlsd=True)#, add_user_slits=add_user_slits)

    lordloc = tslit_dict['lcen']
    rordloc = tslit_dict['rcen']

    # Show in Ginga?
    nslit = lordloc.shape[1]
    print("Found {:d} slits".format(nslit))
    if pargs.show:
        viewer, ch = ginga.show_image(mstrace)
        ginga.show_slits(viewer, ch, lordloc, rordloc,
                         np.arange(nslit) + 1, pstep=50)

    # Output to a MasterFrame?
    if pargs.outfile is not None:
        tslits.save_master(pargs.outfile)

if __name__ == '__main__':
    pargs = parser()
    main(pargs)
