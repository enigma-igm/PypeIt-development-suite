import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt

from IPython import embed

# set environment variables
os.environ['CRDS_PATH'] = '/Users/joe/crds_cache/'
os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu/'
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import stats
#from gwcs import wcstools

## JWST imports
# The calwebb_spec and spec3 pipelines
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Spec2Pipeline
from jwst import datamodels

DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']

# PypeIt imports
from pypeitdev.jwst.jwst_utils import NIRSpecSlitCalibrations, jwst_mosaic, jwst_reduce
from pypeitdev.jwst import jwst_targets
from pypeit.metadata import PypeItMetaData
from pypeit.display import display
from pypeit.images import combineimage
from pypeit import specobjs
from pypeit.utils import inverse, fast_running_median, nan_mad_std

from pypeit.spectrographs.util import load_spectrograph
from pypeit import msgs
from pypeit import spec2dobj


# JFH: This is the main up to date routine. Ignore the others.
DO_NOT_USE = datamodels.dqflags.pixel['DO_NOT_USE']


#disperser = 'PRISM_02073'
# Targets
#target = 'J0020'
#target = 'J0411'
#target = 'J0410'
#target = 'J0313'
#target = 'J0252'
#target = 'J0313'
#target = 'J1007'
#target = 'J1342'
# Dispersers
#disperser = '140H_bogus_F100LP'
#disperser = '140H'
#disperser = '235H'
#disperser = '140H'
#disperser = '395H'
#disperser = 'PRISM_02073'

#progid = '1764'
# J0313-1806
# Slit S200A1
#target, slits, disperser = 'J0313-1806', 'S200A1', '140H'
#target, slits, disperser = 'J0313-1806', 'S200A1', '235H'
#target, slits, disperser = 'J0313-1806', 'S200A1', '395H'
# Slit S200A2
#target, slits, disperser = 'J0313-1806', 'S200A2', '140H'
#target, slits, disperser = 'J0313-1806', 'S200A2', '235H'
#target, slits, disperser = 'J0313-1806', 'S200A2', '395H'


# Define the information to get the exp_list
#progid = '3526'
# Slit S200A1
#target, slits, disperser = 'J1007+2115', 'S200A1', '235H'
#target, slits, disperser = 'J0313-1806', 'S200A1', '235H'
# Slit S200A2
#target, slits, disperser = 'J0038-0653', 'S200A2', '235H'
#target, slits, disperser= 'J0038-1527', 'S200A2', '235H'
#target, slits, disperser = 'J0252-0503', 'S200A2', '235H'
#target, disperser, slits = 'J0410-0139', 'S200A2', '235H'

def validate_redux_input(reduce_input, name):
    """
    Utility method to validate the reduction input keyword parameters and guarantee they are lists. 
    
    Parameters
    ----------
    reduce_input : str or list
        The input to validate.
    name : str
        The name of the parameter.
        
    Returns
    -------
    _reduce_input : list
    """

    if isinstance(reduce_input, str):
        _reduce_input = [reduce_input] 
    elif isinstance(reduce_input, list):
        _reduce_input = reduce_input
    elif reduce_input is None:
        _reduce_input = None
    else: 
        msgs.error(f'{name} must be a string or a list of strings.')
        
    return _reduce_input
    

def validate_interpolatedflat_files(intflat_fs_output_files, intflat_output_files, detector):
    """
    Utility method to validate the interpolated flat files for FS slits.
    
    Parameters
    ----------
    intflat_fs_output_files : list
        List of interpolated flat files for FS slits.
    intflat_output_files : list
        List of interpolated flat files for MOS slits.
    detector : str
        The detector name, i.e. 'nrs1' or 'nrs2'.
    
    Returns
    -------
    merge_fs : bool
        True if there are FS slits to merge, False otherwise
    
    Raises
    ------
    ValueError 
        If the length of the _interpolatedflat_fs.fits files does not match the length of the _interpolatedflat.fits files.
    
    """

    #  intflat for FS slits for NRS1
    if (len(intflat_fs_output_files) > 0) & (len(intflat_fs_output_files) == len(intflat_output_files)):
        merge_fs = True
        msgs.info('Found _nrs1_interpolatedflat_fs.fits files. There are FS slits and MOS slits on nrs1')
    elif len(intflat_fs_output_files) == 0:
        merge_fs = False
    else: 
        raise ValueError('The length of the _interpolatedflat_fs.fits files does not match the length of the _interpolatedflat.fits files for {:s}'.format(detector))
    return merge_fs

class DitherOffsets:
    """
    Compute and hold **detector pixel offsets** per (slit, source) from per-exposure
    **SI-ideal** dither commands, using each slit’s GWCS.

    This class converts the APT/telemetry dither offsets (``XOFFSET``, ``YOFFSET``;
    SI-ideal arcsec) into detector-frame pixel shifts for each requested
    (sl it, source) pair and exposure. For each slit, it samples the per-slit GWCS
    at the **center of the extracted subimage** to measure the local plate scale
    (arcsec/pixel) along detector +X and +Y, then divides the SI-ideal offset by
    those scales:

        dx_pix = XOFFSET_arcsec / (arcsec_per_pixel_along +X)
        dy_pix = YOFFSET_arcsec / (arcsec_per_pixel_along +Y)

    By design here, **detector Y is treated as the spatial direction** and detector
    X as the spectral direction. No projection perpendicular to dispersion is
    applied (i.e. “spatial” ≡ +Y in the detector frame).

    Parameters
    ----------
    cal_data : np.ndarray of jwst.datamodels.MultiSlitModel, shape (2, nexp)
        The calibrated spec2 models for all exposures. Axis 0 is detector
        (0=NRS1, 1=NRS2); axis 1 is exposure index. Each element must have a
        ``.slits`` container with per-slit GWCS (``slit.meta.wcs``) and a
        ``meta.dither`` node carrying ``x_offset``/``y_offset`` (arcsec).
    slits_sources : list[tuple[str, str]]
        The list of **(slit_name, source_name)** pairs to compute. Slit names are
        matched to ``slit.name`` and source names to ``slit.source_name`` (after
        ``str(...).strip()``). Each pair is treated as a unique entry.

    Notes
    -----
    - **SI / SI-ideal**: distortion-removed instrument frame (defined by the SIAF)
    used by JWST/APT to express small commanded offsets; keywords
    ``meta.dither.x_offset`` and ``meta.dither.y_offset`` are in arcsec.
    - The SI-ideal origin and axes are **aperture/slit dependent**. In nod patterns
    that move between slits (e.g., S200A1 ↔ S200A2), the reported ``YOFFSET``
    can differ between those slits even for the same “2-point” dither.
    - The local plate scale is measured from sky separations between
    WCS(x0, y0) and WCS(x0±1, y0) / WCS(x0, y0±1) at the **subimage center**.
    - Coverage is tracked per (detector, exposure) and per (slit, source); some
    entries may exist on only one detector or only a subset of exposures.

    After construction
    ------------------
    - ``offsets_pix_x[slit]``, ``offsets_pix_y[slit]`` :
        dict[str -> ndarray] of shape (2, nexp) — per-detector pixel offsets for
        that slit (0=NRS1, 1=NRS2). NaN where the (slit, source) is absent.
        (Dictionaries are keyed by **slit name**; the corresponding ``source_name``
        is available via ``self.source_names`` aligned with ``self.slit_names``.)
    - ``xoff_arcsec``, ``yoff_arcsec`` :
        ndarray, shape (2, nexp) — global SI-ideal offsets (arcsec) per
        detector/exposure, copied from the data models.
    - ``filters``, ``gratings`` :
        list[str], length ``nexp`` — per-exposure instrument configuration
        (taken from detector 0 for convenience).
    - ``on_detector[slit]`` :
        ndarray(bool), shape (2, nexp) — True where that (slit, source) is present.
    - ``ndet``, ``nexp``, ``det_names``, ``slit_names``, ``source_names`` :
        basic metadata.

    Examples
    --------
    Build once after opening all spec2 ``_cal.fits`` models:

    >>> # cal_data shape: (2, nexp)  [axis 0: 0=NRS1, 1=NRS2]
    >>> pairs = [('S200A1', 'QSO1'), ('S200A2', 'QSO1')]
    >>> offsets = DitherOffsets(cal_data, slits_sources=pairs)

    1) Absolute pixel offsets per (slit, source)
    -----------------------------------------
    Return detector-frame pixel offsets (dx_pix, dy_pix).

    >>> # Both detectors (NaN where absent): shapes (2, nexp)
    >>> dx_all, dy_all = offsets.get_pixel_offsets(slit_name='S200A1', best=False)
    >>> dx_all.shape, dy_all.shape
    ((2, nexp), (2, nexp))

    >>> # Or select by source name (unique within input pairs)
    >>> dx_best, dy_best = offsets.get_pixel_offsets(source_name='QSO1', best=True)
    >>> dx_best.shape, dy_best.shape
    ((nexp,), (nexp,))

    2) Global SI-ideal arcsec offsets (not slit/source-specific)
    ---------------------------------------------------------
    Return the commanded offsets in arcsec for each detector/exposure.

    >>> x_arcsec, y_arcsec = offsets.get_arcsec_offsets()     # shapes (2, nexp)
    >>> x1, y1 = offsets.get_arcsec_offsets(det_name='NRS1')  # shapes (nexp,)

    3) Which detector should be used for a slit/source?
    ------------------------------------------------
    Uses the “good across **all** exposures” rule; ties resolved by ``prefer_det``.

    >>> offsets.get_best_detector(slit_name='S200A1', prefer_det=0)
    (0, 'NRS1')   # example

    4) Relative SI-ideal arcsec offsets (global)
    -----------------------------------------
    Differences every exposure from a reference exposure or external reference.

    >>> rel_x, rel_y = offsets.get_rel_arcsec_offsets(ref_exp=0)  # shapes (2, nexp)
    >>> rel_x0, rel_y0 = offsets.get_rel_arcsec_offsets(external_ref_x=0.0, external_ref_y=0.0)

    5) Relative pixel offsets for a slit/source
    ----------------------------------------
    Differences every exposure from a reference in **pixels** for that entry.

    >>> rdx, rdy = offsets.get_rel_pixel_offsets(slit_name='S200A1', ref_exp=0, best=True)  # (nexp,)
    >>> rdx2, rdy2 = offsets.get_rel_pixel_offsets(source_name='QSO1', ref_exp=0, best=False)  # (2, nexp)

    6) Offsets in the format expected by coadd2d (spatial-only, best detector)
    -----------------------------------------------------------------------
    Convenience wrapper that returns **dy_pix** relative to the reference.

    >>> y_offsets = offsets.get_coadd2d_offsets(slit_name='S200A1', ref_exp=0)  # (nexp,)
    >>> y_offsets_ext = offsets.get_coadd2d_offsets(source_name='QSO1', external_ref=0.0)

    Printing
    --------
    The object’s ``__repr__`` prints a compact table with the global SI-ideal
    X/Y (arcsec) per detector and sample per-pair spatial offsets (dy[pix]) on a
    detector with full coverage when possible.

    Caveats
    -------
    - If the same ``source_name`` appears with multiple slits in ``slits_sources``,
    calls that select by ``source_name`` will pick the **first** match.
    - If the per-slit GWCS is invalid at the subimage center or ±1-pixel steps,
    the computed pixel offsets may be NaN or unreliable.
    """
    def __init__(self, cal_data, slits_sources, verbose=True):
        """
        Parameters
        ----------
        cal_data : np.ndarray of jwst.datamodels.MultiSlitModel, shape (2, nexp)
            The calibrated spec2 models for all exposures. Axis 0 is detector
            (0=NRS1, 1=NRS2); axis 1 is exposure index. Each element must have
            a ``.slits`` container with per-slit WCS (``slit.meta.wcs``) and a
            ``meta.dither`` node carrying ``x_offset``/``y_offset`` (arcsec).
        slits_sources : list of tuples (str, str)
            Slit/Sources names to compute offsets for (e.g., ``[('S200A1', 'QSO1'), ('S200A2', 'QSO1')]``).
            Slit names are matched exactly (after ``str(...).strip()``) to the
            ``slit.name`` entries present in each ``MultiSlitModel``, source names
            are matched exactly after ``str(...).strip()`` to ``slit.source_name``.

        Behavior
        --------
        For each ``slit_source`` and for each (detector, exposure):
          1) Find the slit/source by name in ``cal_data[det, exp].slits``.
          2) Take the subimage center (x0, y0). If ``slit.data`` is absent,
             fall back to the parent model's subarray size (or 2048×2048).
          3) Evaluate the per-slit GWCS at (x0, y0), (x0+1, y0), (x0, y0+1)
             and compute **arcsec/pixel** along +X and +Y via sky separations.
          4) Divide the global SI-ideal offsets for that (det, exp) by those
             arcsec/pixel values to obtain ``dx_pix`` and ``dy_pix``.
          5) Record coverage in ``on_detector[slit][det, exp]``.

        Assumptions & Caveats
        ---------------------
        - “Spatial” ≡ **detector Y**. No projection orthogonal to the dispersion
          is applied.
        - The SI-ideal offsets are small (|X|≪|Y| for standard NIRSpec nods).
        - The method relies on the slit GWCS being valid at (x0, y0) and at
          ±1-pixel steps; if the WCS is invalid there, the pixel offsets will be
          NaN or unreliable.

        Raises
        ------
        ValueError
            If a slit WCS returns outputs that cannot be coerced to
            ``astropy.coordinates.SkyCoord`` for separation calculations.

        Side Effects
        ------------
        Populates:
          ``offsets_pix_x``, ``offsets_pix_y``, ``on_detector``,
          ``xoff_arcsec``, ``yoff_arcsec``, ``filters``, ``gratings``,
          ``filenames``, ``ndet``, ``nexp``, ``det_names``, ``slit_names``.
        """        
        self.cal_data = cal_data
        # Validate slit_source input
        if not isinstance(slits_sources, (list, tuple)) or len(slits_sources) == 0:
            msgs.error('slits_sources must be a non-empty list of (slit_name, source_name) tuples.')
        self.slit_names, self.source_names = zip(*slits_sources)
        self.ndet, self.nexp = cal_data.shape
        self.det_names = ['NRS1', 'NRS2']
    

        # Global per-(det,exp) SI-ideal offsets + per-exposure config
        # SI is a distortion-removed, instrument-centric frame (units: arcsec) defined in the 
        # Science Instrument Aperture File (SIAF). Since in the standard dither pattherns, the 
        # xoff_arcsec is extremely small ~ 1e-6 asec, whereas the yoff_arcsec is significant ~ 1 asec, 
        # we simply treat the x-axis as the spectral direction and the y-axis as the spatial direction. 
        self.xoff_arcsec = np.full((self.ndet, self.nexp), np.nan, float)
        self.yoff_arcsec = np.full((self.ndet, self.nexp), np.nan, float)
        self.filters = []
        self.gratings = []
        self.filenames = np.full((self.ndet, self.nexp), '', str)

        for iexp in range(self.nexp):
            self.filters.append(str(cal_data[0, iexp].meta.instrument.filter))
            self.gratings.append(str(cal_data[0, iexp].meta.instrument.grating))
            for idet in range(self.ndet):
                self.xoff_arcsec[idet, iexp] = float(cal_data[idet, iexp].meta.dither.x_offset)
                self.yoff_arcsec[idet, iexp] = float(cal_data[idet, iexp].meta.dither.y_offset)
                self.filenames[idet, iexp] = str(cal_data[idet, iexp].meta.filename)

        # Per-slit pixel offsets and coverage mask
        self.offsets_pix_x = {s: np.full((self.ndet, self.nexp), np.nan, float) for s in self.slit_names}
        self.offsets_pix_y = {s: np.full((self.ndet, self.nexp), np.nan, float) for s in self.slit_names}
        self.on_detector   = {s: np.zeros((self.ndet, self.nexp), bool)         for s in self.slit_names}

        # Compute per-slit arcsec/pix (x,y) from GWCS at slit center; convert global SI-ideal to pixels
        for slit_name, source_name in zip(self.slit_names, self.source_names):
            for idet in range(self.ndet):
                for iexp in range(self.nexp):
                    model = cal_data[idet, iexp]
                    names_here = [sl.name for sl in model.slits]
                    try:
                        islit = next(i for i, nm in enumerate(names_here) if str(nm).strip() == slit_name.strip())
                    except StopIteration:
                        msgs.warn(f"Slit: {slit_name} not on detector NRS{idet+1}, exp={iexp}.")
                        continue  # slit not present on this (det,exp)

                    slit = model.slits[islit]
                    wcs  = slit.meta.wcs

                    if getattr(slit, "data", None) is not None and slit.data is not None:
                        ny, nx = slit.data.shape
                    else:
                        msgs.warn(f"Slit {slit_name} on det=NRS{idet+1}, exp={iexp} has no data array; assuming 2048^2.")
                        ny = getattr(model.meta.subarray, "ysize", 2048)
                        nx = getattr(model.meta.subarray, "xsize", 2048)
                    x0, y0 = (nx - 1)/2.0, (ny - 1)/2.0

                    # evaluate WCS at subimage center and +1 pix in x and y
                    c0, lam0  = wcs(x0,       y0,       with_units=True)
                    cx, lamx1 = wcs(x0 + 1.0, y0,       with_units=True)
                    cy, lamy1 = wcs(x0,       y0 + 1.0, with_units=True)

                    asec_per_pix_x = c0.separation(cx).to(u.arcsec).value
                    asec_per_pix_y = c0.separation(cy).to(u.arcsec).value

                    self.offsets_pix_x[slit_name][idet, iexp] = self.xoff_arcsec[idet, iexp] / asec_per_pix_x
                    self.offsets_pix_y[slit_name][idet, iexp] = self.yoff_arcsec[idet, iexp] / asec_per_pix_y
                    self.on_detector[slit_name][idet, iexp]   = True
                    
        if verbose: 
            text = repr(self)
            msgs.info("\n" + text)

    def _check_name(self, slit_name=None, source_name=None):
        """ 
        Utility method to check that the slit name exists in the list of slit names.
        
        Parameters
        ----------
        slit_name : str, optional
            The slit name. Either slit_name or source_name must be provided.
        source_name : str, optional
            The source name. Either slit_name or source_name must be provided.
    
        Returns
        -------
        idx : int
            The index of the slit/source in the list of slit/source names.
        _slit_name : str
            The slit name.
        _source_name : str
            The source name.
        """
        if (slit_name is None) and (source_name is None):
            msgs.error('Either slit or source name must be provided.')
        if slit_name is not None:
            if slit_name not in self.slit_names:
                msgs.error(f'Slit name {slit_name} not found. Available slits are: {self.slit_names}')
            else: 
                _slit_name = str(slit_name)
                idx = self.slit_names.index(_slit_name)
                _source_name = self.source_names[idx]
                return idx, _slit_name, _source_name
        if source_name is not None:
            if source_name not in self.source_names:
                msgs.error(f'Source name {source_name} not found. Available sources are: {self.source_names}')
            else: 
                _source_name = str(source_name)
                idx = self.source_names.index(_source_name)
                _slit_name = self.slit_names[idx]
                return idx, _slit_name, _source_name
            

    def get_pixel_offsets(self, slit_name=None, source_name=None, best=False):
        """
        Return the (dx_pix, dy_pix) pixel offsets for a given slit. 

        Parameters
        ----------
        slit_name : str
            The slit name, optional. Either slit_name or source_name must be provided. 
        source_name : str
            The source name, optional. Either slit_name or source_name must be provided.
        best : bool, optional
            If True, return only the offsets for the best detector (see `get_best_detector`).
            If False (default), return offsets for both detectors, with NaNs where the slit is not present on 
            that detector.
            
        Returns
        -------
        dx_pix : np.ndarray, shape (ndet, nexp) or (nexp,) if best=True
            The pixel offsets in the x direction.
        dy_pix : np.ndarray, shape (ndet, nexp) or (nexp,) if best=True
            The pixel offsets in the y direction.
        """
        idx, _slit_name, _source_name = self._check_name(slit_name=slit_name, source_name=source_name)
        if best: 
            idet, det_name = self.get_best_detector(_slit_name)
            if np.isnan(idet):
                msgs.error(f'Slit {slit_name} has no coverage on either detector for any exposure.')
            return self.offsets_pix_x[_slit_name][idet, :], self.offsets_pix_y[_slit_name][idet, :]
        else: 
            return self.offsets_pix_x[_slit_name].copy(), self.offsets_pix_y[_slit_name].copy()
        

    def get_arcsec_offsets(self, idet=None, det_name=None):
        """
        Return the global offsets (arcsec) for both detectors or the best detector.
        
        Parameters
        ----------
        idet : int, optional
            The detector index (0 or 1) to return offsets for. If both idet and det_name are None, return offsets for both detectors. If both idet and det_name are specified, raise an error.
        det_name : str, optional
            The detector name ('NRS1' or 'NRS2') to return offsets for. If both idet and det_name are None, return offsets for both detectors. If both idet and det_name are specified, raise an error.
        
        Returns
        -------
        dx_arcsec : np.ndarray, shape (ndet, nexp) or (nexp,)
            The offsets in the x (spectral) direction (arcsec).
        dy_arcsec : np.ndarray, shape (ndet, nexp) or (nexp,)
            The offsets in the y (spatial)  direction (arcsec).
        """
        if (idet is not None) and (det_name is not None):
            msgs.error('Specify only one of idet or det_name.')
        if idet is not None:
            return self.xoff_arcsec[idet, :].copy(), self.yoff_arcsec[idet, :].copy()
        if det_name is not None:
            det_name_u = det_name.upper()
            if det_name_u == 'NRS1':
                idet = 0
            elif det_name_u == 'NRS2':
                idet = 1
            else:
                msgs.error('det_name must be NRS1 or NRS2.')
            return self.xoff_arcsec[idet, :].copy(), self.yoff_arcsec[idet, :].copy()
        return self.xoff_arcsec.copy(), self.yoff_arcsec.copy()

    def get_best_detector(self, slit_name=None, source_name=None, prefer_det=0):
        """
        Determine the best detector for a slit/source for computing the offsets, using the detector-selection rule:
            1. If the slit/source is on only one detector => use that detector.
            2. If the slit/source is on both detectors => choose `prefer_det` (0=NRS1, 1=NRS2).

        Parameters
        ----------
        slit_name : str
            The slit name. Optional, either slit_name or source_name must be provided.
        source_name : str
            The source name. Optional, either slit_name or source_name must be provided.
        prefer_det : int, optional
            The preferred detector if the slit is on both detectors (0=NRS1, 1=NRS2). Default is 0 (NRS1).

        Returns
        -------
        idet : int
            The chosen detector index (0 or 1), or np.nan if the slit is not on either detector. If the slit is not on either detector, idet is np.nan. 
        det_name : str or None
            The chosen detector name ('NRS1' or 'NRS2'), or None if the slit is not on either detector. If the slit is not on either detector, det_name is None.    
        """
        
        if prefer_det not in (0, 1):
            msgs.error('prefer_det must be 0 (NRS1) or 1 (NRS2).')
        
        idx, _slit_name, _source_name = self._check_name(slit_name=slit_name, source_name=source_name)

        det_good = np.all(self.on_detector[_slit_name], axis=1)
        # If only one detector has any coverage, use that one, otherwise just use NRS1
        if not det_good.any():
            msgs.warn(f'Slit/Source {_slit_name}/{_source_name} has no coverage on either detector for any exposure.')
            return np.nan, None
        
        if det_good.sum() == 1:
            idet = int(np.argmax(det_good))
            return idet, self.det_names[idet]
        
        # both True -> honor prefer_det
        idet = int(prefer_det)
        return idet, self.det_names[idet]
    
    def get_rel_arcsec_offsets(self, ref_exp=0, external_ref_x=None, external_ref_y=None):
        """
        Relative SI-ideal offsets (arcsec) over exposures, using the same detector rule,
        but these are global (not slit-specific).

        Parameters
        ----------
        ref_exp : int, optional
            The reference exposure index. Default is 0, i.e. the first exposure in the series. If external_reference is provided, this is ignored.
        external_ref_x : float or np.ndarray of shape (2,) 
            If provided, use this external reference offsets (arcsec) instead of a set of offsets from one of the exposures 
            (ref_exp) in this sequence of exposures as the reference. If a float is provided, it is used for both detectors.
        external_ref_y : float or np.ndarray of shape (2,) 
            If provided, use this external reference offsets (arcsec) instead of a set of offsets from one of the exposures 
            (ref_exp) in this sequence of exposures as the reference. If a float is provided, it is used for both detectors.
        
        Returns
        -------
        rel_arcsec_offsets_x : np.ndarray, shape (2, nexp)
            The relative x arcsec offsets for both detectors (NRS1, NRS2).
        rel_arcsec_offsets_y : np.ndarray, shape (2, nexp)            
            The relative y arcsec offsets for both detectors (NRS1, NRS2).
        """
        if (external_ref_x is None) and (external_ref_y is None):
            if (ref_exp < 0) or (ref_exp >= self.nexp):
                msgs.error(f'ref_exp must be between 0 and {self.nexp-1}.')
            ref_val_x = self.xoff_arcsec[:, ref_exp]
            ref_val_y = self.yoff_arcsec[:, ref_exp]
        else:
            if (external_ref_x is None) or (external_ref_y is None):
                msgs.error('If providing external references, provide both external_ref_x and external_ref_y.')
            # normalize to (2,) arrays
            def _to_det2(val, label):
                if isinstance(val, (float, int)):
                    return np.full(2, float(val), float)
                val = np.asarray(val, dtype=float)
                if val.shape != (2,):
                    msgs.error(f'{label} must be a scalar or shape (2,) array.')
                return val
            ref_val_x = _to_det2(external_ref_x, 'external_ref_x')
            ref_val_y = _to_det2(external_ref_y, 'external_ref_y')

        return ref_val_x[:, np.newaxis] - self.xoff_arcsec, ref_val_y[:, np.newaxis] - self.yoff_arcsec 


    def get_rel_pixel_offsets(self, slit_name=None, source_name=None, ref_exp=0, external_ref_x=None, external_ref_y=None, prefer_det=0, best=True):
        """
        Relative SI-ideal offsets (arcsec) over exposures, using the same detector rule,
        but these are global (not slit-specific).

        Parameters
        ----------
        slit_name : str
            The slit name. Optional, either slit_name or source_name must be provided.
        source_name : str
            The source name. Optional, either slit_name or source_name must be provided.
        ref_exp : int, optional
            The reference exposure index. Default is 0, i.e. the first exposure in the series. If external_reference is provided, this is ignored.
        external_ref_x : float or np.ndarray of shape (2,) 
            If provided, use this external reference offsets (arcsec) instead of a set of offsets from one of the exposures 
            (ref_exp) in this sequence of exposures as the reference. If a float is provided, it is used for both detectors.
        external_ref_y : float or np.ndarray of shape (2,) 
            If provided, use this external reference offsets (arcsec) instead of a set of offsets from one of the exposures 
            (ref_exp) in this sequence of exposures as the reference. If a float is provided, it is used for both detectors.
        prefer_det : int, optional
            The preferred detector if both detectors have coverage (0=NRS1, 1=NRS2). Default is 0 (NRS1).
            This is passed to `get_best_detector
        best : bool, optional
            If True (default), return one set of offsets for both detectors using the best detector (see `get_best_detector`). 
            If False, return offsets for both detectors, with NaNs where the slit is not present on that detector. 

        Returns
        -------
        rel_arcsec_offsets (nexp,) float
        """
        idx, _slit_name, _source_name = self._check_name(slit_name=slit_name, source_name=source_name)
   
        # Check external_reference is provided, if best=True it is a float, and if best=False, it is a 2d array
        if external_ref_x is None and external_ref_y is None:
            if (ref_exp < 0) | (ref_exp >= self.nexp):
                msgs.error(f'ref_exp must be between 0 and {self.nexp-1}.')
            ref_val_x = self.offsets_pix_x[_slit_name][:, ref_exp]
            ref_val_y = self.offsets_pix_y[_slit_name][:, ref_exp]
        else: 
            if (external_ref_x is None) or (external_ref_y is None):
                msgs.error('If providing external references, provide both external_ref_x and external_ref_y.')
            def _to_det2(val, label):
                if isinstance(val, (float, int)):
                    return np.full(2, float(val), float)
                val = np.asarray(val, dtype=float)
                if val.shape != (2,):
                    msgs.error(f'{label} must be a scalar or shape (2,) array.')
                return val
            ref_val_x = _to_det2(external_ref_x, 'external_ref_x')
            ref_val_y = _to_det2(external_ref_y, 'external_ref_y')

        rel_pix_offsets_x = ref_val_x[:, np.newaxis] - self.offsets_pix_x[_slit_name]
        rel_pix_offsets_y = ref_val_y[:, np.newaxis] - self.offsets_pix_y[_slit_name]

        if best: 
            idet, det_name = self.get_best_detector(slit_name=_slit_name, prefer_det=prefer_det)
            if np.isnan(idet):
                msgs.error(f'Slit/Source {_slit_name}/{_source_name} has no coverage on either detector for any exposure.')
            return rel_pix_offsets_x[idet, :], rel_pix_offsets_y[idet, :]

        return rel_pix_offsets_x, rel_pix_offsets_y
        
    def get_coadd2d_offsets(self, slit_name=None, source_name=None, ref_exp=0, external_ref=None, prefer_det=0):
        """
        Wrapper for `get_rel_pixel_offsets` to return the offsets in the format required by `core.coadd2d`.
        
        Parameters
        ----------
        slit_name : str
            The slit name. Optional, either slit_name or source_name must be provided.
        source_name : str
            The source name. Optional, either slit_name or source_name must be provided.
        ref_exp : int, optional
            The reference exposure index. Default is 0, i.e. the first exposure in the series. If external_reference is provided, this is ignored.
        external_ref : float 
            If provided, use this external reference offset, which is the y_offset (spatial in pixels) instead of an offset from 
            one of the exposures (ref_exp) in this sequence of exposures as the reference.
        prefer_det : int, optional
            The preferred detector if both detectors have coverage (0=NRS1, 1=NRS2). Default is 0 (NRS1).
            This is passed to `get_best_detector

        """
        if external_ref is not None:
            external_ref_x = 0.0
            external_ref_y = float(external_ref)
        else:
            external_ref_x = None
            external_ref_y = None
            
        return self.get_rel_pixel_offsets(slit_name=slit_name, source_name=source_name, ref_exp=ref_exp, external_ref_x=external_ref_x, 
                                          external_ref_y=external_ref_y, prefer_det=prefer_det, best=True)[1]        
    
    def __repr__(self):
        # ---------- helpers ----------
        def _fmt_row(arr, width=8, prec=3):
            # arr is 1D (nexp,); show NaN as '·'
            out = []
            for v in np.asarray(arr, float):
                if np.isfinite(v):
                    out.append(f"{v:{width}.{prec}f}")
                else:
                    out.append(" " * (width-1) + "·")
            return " ".join(out)

        # layout
        sep = "-" * max(72, 20 + 10*self.nexp)
        nitems = len(self.slit_names)
        header = (
            "\n" + "=" * 50 + "\n"
            + f"DitherOffsets  ndet={self.ndet}  nexp={self.nexp}  nslit/source pairs={nitems}\n"
            + "=" * 50
        )

        # single config for all (as constructed)
        cfg = "Config: "
        if len(self.filters) and len(self.gratings):
            cfg += f"{self.filters[0]}/{self.gratings[0]}"
        else:
            cfg += "(unknown)"

        exp_line = "Exp idx:           " + " ".join(f"{i:^8d}" for i in range(self.nexp))

        # global SI-ideal X/Y (arcsec), per detector
        x1 = _fmt_row(self.xoff_arcsec[0, :])
        y1 = _fmt_row(self.yoff_arcsec[0, :])
        x2 = _fmt_row(self.xoff_arcsec[1, :])
        y2 = _fmt_row(self.yoff_arcsec[1, :])

        lines = [
            header, cfg, sep, exp_line,
            f"NRS1  X[arcsec]:   {x1}",
            f"NRS1  Y[arcsec]:   {y1}",
            f"NRS2  X[arcsec]:   {x2}",
            f"NRS2  Y[arcsec]:   {y2}",
            sep
        ]

        # show first few slit/source pairs: dy[pix] only, pick a detector with full coverage if possible
        show_n = min(8, nitems)
        for slit, src in list(zip(self.slit_names, self.source_names))[:show_n]:
            arr = self.offsets_pix_y[slit]  # shape (2, nexp) — keyed by slit
            full = np.all(np.isfinite(arr), axis=1)   # (2,) fully covered across all exps?
            any_  = np.any(np.isfinite(arr), axis=1)  # (2,) any coverage?

            if full.any():
                det = 0 if full[0] else 1
            elif any_.any():
                det = 0 if any_[0] else 1
            else:
                lines.append(f"{slit}/{src:<12}: (no coverage)")
                continue

            dy = _fmt_row(arr[det, :])
            label = f"{slit}/{src}"
            lines.append(f"{label:<24}: {self.det_names[det]} dy[pix]: {dy}")

        if nitems > show_n:
            lines.append(f"... ({nitems - show_n} more slit/source pairs)")

        return "\n".join(lines)


# def compute_spatial_positions(cal_data, slit_names):
#     """
#     Routine to compute the spatial offsets in pixels for a given slit from the dither pattern offsets. 
    
#     Parameters
#     ----------
#     cal_data : np.narray of datamodels.SpecModel objects
#         Calibrated data models from the Spec2 pipeline, shape = (2, nexp), i.e. the first
#         dimension is the detector (0 = NRS1, 1 = NRS2) and the second dimension is the exposure number.
#     slit_names : list of str
#         List of slit names to compute the offsets for.

#     Returns
#     -------
#     offsets_dict : dict
#         Dictionary with the slit names as keys. The values are np.ndarrays of shape (2, nexp,) with the spatial offsets 
#         in pixels for each exposure
    
#     Simplest possible offsets:
#       - per-slit WCS (from cal_data[det,exp].slits[index])
#       - arcsec_per_pix_y via sky separation between (x0,y0) and (x0,y0+1)
#       - dy_pix = meta.dither.y_offset / arcsec_per_pix_y
#       - relative to (ref_det, ref_exp) -> NRS1, exposure 0 by default

#     Returns
#     -------
#     offsets_nrs1 : (nexp,) float
#         Detector-Y spatial offsets in pixels for NRS1, relative to NRS1 exp0.
#     """

#     nslits = len(slit_names)    
#     ndet, nexp = cal_data.shape
#     dx_pix = np.full((nslits, ndet, nexp), np.nan, dtype=float)
#     dy_pix = np.full((nslits, ndet, nexp), np.nan, dtype=float)
#     on_detector = np.zeros((nslits, ndet, nexp), dtype=bool)
#     # Setup parameters 
#     filter = cal_model.meta.instrument.filter
#     grating = cal_model.meta.instrument.grating

#     # These offsets are global, i.e. the are the same for all slits and correspond to the dither patern. 
#     xoff_arcsec = float(cal_model.meta.dither.x_offset) 
#     yoff_arcsec = float(cal_model.meta.dither.y_offset)

#     for islit, slit_name in enumerate(slit_names):   
#         for idet in range(ndet):
#             for iexp in range(nexp):
#                 cal_model = cal_data[idet, iexp]
#                 slit_names = np.array([slit.name for slit in cal_model.slits])
#                 _indx = np.where((slit_names == slit_name))[0]
#                 on_detector = _indx.size != 0
#                 if on_detector:                 
#                     slit = cal_model.slits[_indx[0]]
#                     wcs_obj = slit.meta.wcs  # per-slit 2-input GWCS (x,y)->(ra,dec,lam)
#                     ny, nx = slit.data.shape
#                     # compute everything relative to the middle of the subimage
#                     x0, y0 = (nx - 1)/2.0, (ny - 1)/2.0
#                     # arcsec per pixel along detector-Y at the slit center
#                     c0,  lam0  = wcs_obj(x0, y0, with_units=True)
#                     yc1, lamx1 = wcs_obj(x0, y0 + 1.0, with_units=True) 
#                     xc1, lamy1 = wcs_obj(x0 + 1.0, y0, with_units=True)  
#                     asec_per_pix_y = c0.separation(yc1).to(u.arcsec).value
#                     asec_per_pix_x = c0.separation(xc1).to(u.arcsec).value
#                     # convert to detector pixels along Y
#                     dy_pix[idet, iexp] = yoff_arcsec/asec_per_pix_y
#                     dx_pix[idet, iexp] = xoff_arcsec/asec_per_pix_x
#         # relative to NRS1, exposure 0 (your original convention: ref - current)
#         ref_y = dy_pix[ref_det, ref_exp]
#         #offsets_nrs1 = ref_y - dy_pix[0, :]   # shape (nexp,)
#         offsets_nrs1 = ref_y - dy_pix # shape (ndet, nexp)
#         #print(f"Slit_name = {slit_name}, xoff (arcsec) = {xoff_arcsec}, yoff (arcsec) = {yoff_arcsec}, yoffsets (pixels): {offsets_nrs1[0, :]}")   
#     return offsets_nrs1[0, :]



def jwst_run_redux(redux_dir, source_type, uncal_list=None, rate_list=None, 
                   reduce_slits=None, reduce_sources=None, 
                   show=False, overwrite_stage1=False, overwrite_stage2=False, 
                   kludge_err=1.5, snr_thresh=10.0, find_trim_edge=[5,5], bkg_redux=False, run_bogus_f100lp=False):
    """
    Main routine to reduce JWST NIRSpec data
    
    Parameters
    ----------
    redux_dir : str
        Path to the directory where the data will be reduced.
    source_type = str
        source_type for the spec2d pipeline. Options are 'POINT' or 'EXTENDED'. Set this parameter to be 'POINT' for fixed slit data reduction, and 
        'EXTENDED' for MSA reductions. This impacts the flux calibration, I believe via slit loss corrections and what is assumed.  
    uncal_list : list
        List of lists of uncalibrated files for each exposure. exp_list[0] is for nrs1 and exp_list[1] is for nrs2.  Optional, default is None. Either uncal_list or rate_list must be provided.
    rate_list : list
        List of lists of rate files for each exposure. rate_list[0] is for nrs1 and rate_list[1] is for nrs2. Optional, default is None. Either uncal_list or rate_list must be provided.
    reduce_slits : str or list, optional
        List of slits to reduce. If None, reduce all. 
    reduce_sources : str or list, optional
        List of sources to reduce. If None, reduce all. 
    show : bool, optional
        Show the images and QA for the reduction in the ginga window. 
    overwrite_stage1 : bool, optional
        Rerun the jwst pipeline overwriting existing stage1 reductions 
    overwrite_stage2 : bool, optional
        Rerun the jwst pipeline overwriting existing stage2 reductions
    kludge_err : float, optional
        Factor to scale the sigma error maps up by to account for the incorrect error propagation in the JWST pipeline. 
    snr_thresh : float, optional
        SNR threshold for the spec2d extraction. Default is 10.0. 
    find_trim_edge : list, optional
        Trim the slit by this number of pixels left/right before finding objects. 
    bkg_redux : bool, optional
        If True, perform a background redux using image differencing. If False, model the background with a bspline. bkg_redux should typically be set to True for MSA reductions, 
        and False for FS reductions.
    run_bogus_f100lp : bool
        If True, run bogus redux for the 140H/F100LP grating to accomodate data taken with the F070LP filter? Default is False. 
        
    Returns
    -------
    DitherOffsets : DitherOffsets object
        Object containing the dither offsets and slit/source information.
    """

    _reduce_slits = validate_redux_input(reduce_slits, 'reduce_slits')
    _reduce_sources = validate_redux_input(reduce_sources, 'reduce_sources')
    
    # Define the path
    output_dir = os.path.join(redux_dir, 'output')
    output_dir_level1 = os.path.join(output_dir, 'level1')
    output_dir_level2 = os.path.join(output_dir, 'level2')
    pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

    # Make the new Science dir
    # TODO: This needs to be defined by the user
    scipath = os.path.join(pypeit_output_dir, 'Science')
    if not os.path.isdir(scipath):
        msgs.info('Creating directory for Science output: {0}'.format(scipath))
        os.makedirs(scipath)
    if not os.path.isdir(output_dir):
        msgs.info('Creating directory for calwebb output: {0}'.format(output_dir))
        os.makedirs(output_dir)
    if not os.path.isdir(output_dir_level1):
        msgs.info('Creating directory for calwebb level 1 output: {0}'.format(output_dir_level1))
        os.makedirs(output_dir_level1)
    if not os.path.isdir(output_dir_level2):
        msgs.info('Creating directory for calwebb level 2 output: {0}'.format(output_dir_level2))
        os.makedirs(output_dir_level2)
    
    # Did the user pass in an uncal_list? 
    if uncal_list is None and rate_list is None:
        msgs.error('Either uncal_list or rate_list must be provided.')
    elif uncal_list is not None and rate_list is not None:
        msgs.error('Only one of uncal_list or rate_list can be provided.')
    elif uncal_list is not None and rate_list is None: 
        uncalfiles_1 = uncal_list[0]
        uncalfiles_2 = uncal_list[1] if len(uncal_list) > 1 else []
        uncalfiles_all = uncalfiles_1 + uncalfiles_2
        nexp = len(uncalfiles_1)

        basenames, basenames_1, basenames_2, rate_files_1, rate_files_2 = [], [], [], [], []
        for sci1, sci2 in zip(uncalfiles_1, uncalfiles_2):
            b1 = os.path.basename(sci1).replace('_uncal.fits', '')
            b2 = os.path.basename(sci2).replace('_uncal.fits', '')
            basenames_1.append(b1)
            basenames_2.append(b2)
            basenames.append(b2.replace('_nrs2', ''))
            rate_files_1.append(os.path.join(output_dir_level1, b1 + '_rate.fits'))
            rate_files_2.append(os.path.join(output_dir_level1, b2 + '_rate.fits'))



        #parameter_dict_det1 = {"jump": {"maximum_cores": 'quarter'},}
        # These clean_flicker_noise parameters for NIRSpec are based on the recommended values here: 
        # https://jwst-docs.stsci.edu/known-issues-with-jwst-data/1-f-noise#gsc.tab=0 
        # which were created after JWST moved the 1/f noise correction to the level1 correction stage
        parameter_dict_det1 = {'clean_flicker_noise': 
            {'skip': False, 'fit_method': 'fft', 'n_sigma': 2, 'mask_science_regions': False, 'background_method': None}}
        for uncal in uncalfiles_all:
            ratefile = os.path.join(output_dir_level1,  os.path.basename(uncal).replace('_uncal', '_rate'))
            if os.path.isfile(ratefile) and not overwrite_stage1:
                msgs.info('Using existing rate file: {0}'.format(ratefile))
                continue
            Detector1Pipeline.call(uncal, save_results=True, output_dir=output_dir_level1, steps=parameter_dict_det1)

    elif uncal_list is None and rate_list is not None:
        rate_files_1 = rate_list[0]
        rate_files_2 = rate_list[1] if len(rate_list) > 1 else []
        nexp = len(rate_files_1)
        
        basenames, basenames_1, basenames_2 = [], [], []
        for rate1, rate2 in zip(rate_files_1, rate_files_2):
            b1 = os.path.basename(rate1).replace('_rate.fits', '')
            b2 = os.path.basename(rate2).replace('_rate.fits', '')
            basenames_1.append(b1)
            basenames_2.append(b2)
            basenames.append(b2.replace('_nrs2', ''))


    rate_files_all = rate_files_1 + rate_files_2
    bkg_indices = [(1,2), (0,2), (0,1)]

    # TODO Should we flat field. The flat field and flat field error are wonky and probably nonsense
    param_dict_spec2 = {
        'assign_wcs': {'save_results': True}, # This output now has the full 2d image 
        'extract_2d': {'save_results': True},
        'bkg_subtract': {'skip': True},
        'imprint_subtract': {'save_results': True}, # TODO Check up on whether imprint subtraction is being done by us???
        'master_background_mos': {'skip': True},
        'srctype': {'source_type': source_type},        
        # Default to setting the source type to extended for MSA data and point for FS data. This impacts flux calibration.
        #'srctype': {'source_type': 'POINT'} if 'FS' in mode else {'source_type': 'EXTENDED'},
        # 'flat_field': {'skip': True},
        'resample_spec': {'skip': True},
        'extract_1d': {'skip': True},
        'flat_field': {'save_interpolated_flat': True}, 
         # Forces to always run the barshadow step. Default is to apply  barshadow only for extended sources, which means slit won't be flat
        'barshadow': {'source_type': 'EXTENDED'},  
#        'nsclean': {'skip': True, 'save_results': False},
    }
    # So the nsclean is now being done via clean_flicker_noise in the Det1 pipeline. 

    # TODO I'm rather unclear on what to do with the src_type since I'm not following what the JWST ipeline is doing there very
    # well. I think they are trying to model slit losses for point sources but not for extended sources. Changing src_type
    # changes the output units from MJy/sr to MJy. It seems more intuititive that to have the 2d images fluxed in SB
    # units. Right now I'm leaning towards use srctype='EXTENDED' for MSA but srctype='POINT' for FS.

    # For MSA data, we need to run the MSA flagging step and this generates the full 2d frames we operate on, for
    # FS data, the MSA flagging is not performed and so the last step is the assign_wcs step. 
    # This is no longer needed now that nslcean is after both. 
    #if mode =='MSA':
    #    param_dict_spec2['msa_flagging'] = {'save_results': True}
    #elif mode == 'FS':
    #    param_dict_spec2['assign_wcs'] = {'save_results': True}
    #param_dict_spec2['nsclean'] = {'save_results': True, 'skip': False}
    # Run the spec2 pipeline

    for sci in rate_files_all:
        assign_wcs_file = os.path.join(output_dir_level2, os.path.basename(sci).replace('_rate', '_assign_wcs'))
        if os.path.isfile(assign_wcs_file) and not overwrite_stage2:
            msgs.info('Found existing assign_wcs file: {0}; not running Spec2'.format(assign_wcs_file))
            continue
        Spec2Pipeline.call(sci, save_results=True, output_dir=output_dir_level2, steps=param_dict_spec2)


    # Some pypeit things
    spectrograph = load_spectrograph('jwst_nirspec')
    par = spectrograph.default_pypeit_par()
    det_container_list = [spectrograph.get_detector_par(1), spectrograph.get_detector_par(2)]
    fitstbl_1 = PypeItMetaData(spectrograph, par=par,files=rate_files_1, strict=True)
    fitstbl_2 = PypeItMetaData(spectrograph, par=par,files=rate_files_2, strict=True)
    fitstbls = [fitstbl_1, fitstbl_2]

    pypeline = 'MultiSlit'
    par['rdx']['redux_path'] = pypeit_output_dir
    qa_dir = os.path.join(pypeit_output_dir, 'QA')
    par['rdx']['qadir'] = 'QA'
    png_dir = os.path.join(qa_dir, 'PNGs')
    if not os.path.isdir(qa_dir):
        msgs.info('Creating directory for QA output: {0}'.format(qa_dir))
        os.makedirs(qa_dir)
    if not os.path.isdir(png_dir):
        os.makedirs(png_dir)

    # Set some parameters for difference imaging
    if bkg_redux:
        par['reduce']['findobj']['skip_skysub'] = True # Do not sky-subtract when object finding
        par['reduce']['extraction']['skip_optimal'] = True # Skip local_skysubtraction and profile fitting

    # Turn off 2d model masking
    par['reduce']['extraction']['use_2dmodel_mask'] = False
    # Set the SNR threshold for the spec2d object finding
    par['reduce']['findobj']['snr_thresh'] = snr_thresh
    par['reduce']['findobj']['find_fwhm'] = 2.0 # JWST data are often undersampled
    par['reduce']['findobj']['find_trim_edge'] = find_trim_edge


    # Output file names
    intflat_output_files_1 = []
    msa_output_files_1 = []
    cal_output_files_1 = []

    intflat_output_files_2 = []
    msa_output_files_2 = []
    cal_output_files_2 = []

    for base1, base2 in zip(basenames_1, basenames_2):
        msa_output_files_1.append(os.path.join(output_dir_level2, base1 + '_assign_wcs.fits'))
        msa_output_files_2.append(os.path.join(output_dir_level2, base2 + '_assign_wcs.fits'))

        intflat_output_files_1.append(os.path.join(output_dir_level2, base1 + '_interpolatedflat.fits'))
        cal_output_files_1.append(os.path.join(output_dir_level2, base1 + '_cal.fits'))
        intflat_output_files_2.append(os.path.join(output_dir_level2, base2 + '_interpolatedflat.fits'))
        cal_output_files_2.append(os.path.join(output_dir_level2, base2 + '_cal.fits'))

    # Check to see if there are also _interpolatedflat_fs.fits files, since the Spec2 pipeline writes these out separately for MOS + FS data
    # If they exist, we need to append these slits to the  MOS slits in the _interpolatedflat.fits files that we read in
    intflat_fs_output_files_1, intflat_fs_output_files_2 = [], []
    for file_1, file_2 in zip(intflat_output_files_1, intflat_output_files_2):
        if os.path.isfile(file_1.replace('_interpolatedflat.fits', '_interpolatedflat_fs.fits')):
            intflat_fs_output_files_1.append(file_1.replace('_interpolatedflat.fits', '_interpolatedflat_fs.fits'))
        if os.path.isfile(file_2.replace('_interpolatedflat.fits', '_interpolatedflat_fs.fits')):
            intflat_fs_output_files_2.append(file_2.replace('_interpolatedflat.fits', '_interpolatedflat_fs.fits'))
    # Validate 
    merge_fs_nrs1 = validate_interpolatedflat_files(intflat_fs_output_files_1, intflat_output_files_1, 'nrs1')
    merge_fs_nrs2 = validate_interpolatedflat_files(intflat_fs_output_files_2, intflat_output_files_2, 'nrs2')


    # Read in calwebb outputs for everytihng

    # Read in multi exposure calwebb outputs
    msa_multi_list_1 = []
    intflat_multi_list_1 = []
    final_multi_list_1 = []
    msa_multi_list_2 = []
    intflat_multi_list_2 = []
    final_multi_list_2 = []
    #nslits_1 = np.zeros(nexp, dtype=int)
    #nslits_2 = np.zeros(nexp, dtype=int)
    #t_eff = np.zeros(nexp, dtype=float)

    ndetectors = 2
    # Create arrays to hold JWST spec2, but only load the files when they're needed
    msa_data = np.empty((ndetectors, nexp), dtype=object)
    flat_data = np.empty((ndetectors, nexp), dtype=object)
    cal_data = np.empty((ndetectors, nexp), dtype=object)

    ## OLD OFFSETS CODE

    # Get the dither offsets in pixels
    #dither_offsets = np.zeros((ndetectors,nexp), dtype=float)
    # TODO: This probably isn't correct.  I.e., need to know offsets and slit position angle.
    #for iexp in range(nexp):
    #    with fits.open(rate_files_1[iexp]) as hdu:
    #        dither_offsets[0,iexp] = hdu[0].header['YOFFSET']
    #for idet in range(1,ndetectors):
    #    dither_offsets[idet] = dither_offsets[0]
    #dither_offsets_pixels_old = dither_offsets.copy()
    #for idet in range(ndetectors):
    #    dither_offsets_pixels_old[idet] /= det_container_list[idet].platescale
    ## NOTE: Sign convention requires this calculation of the offset
    #dither_offsets_pixels_old = dither_offsets_pixels_old[:,0,None] - dither_offsets_pixels_old
    #print(dither_offsets_pixels_old[0,:])
        

    print('Reading in calwebb outputs. This may take a while...')

    # TODO Figure out why this is so damn slow! I suspect it is calwebb1
    for iexp in range(nexp):
        # Open some JWST data models
        msa_data[0, iexp] = datamodels.open(msa_output_files_1[iexp])
        cal_data[0, iexp] = datamodels.open(cal_output_files_1[iexp])
        flat_data_1 = datamodels.open(intflat_output_files_1[iexp])
        # Merge the FS slits into the MOS slits for NRS1      
        if merge_fs_nrs1: 
            msgs.info('Appending interpolatedflat FS slits into MOS output for {:s}'.format(basenames_1[iexp]))
            flat_data_fs_1 = datamodels.open(intflat_fs_output_files_1[iexp])
            for slit in flat_data_fs_1.slits:
                flat_data_1.slits.append(slit)
        flat_data[0, iexp] = flat_data_1


        msa_data[1, iexp] = datamodels.open(msa_output_files_2[iexp])
        cal_data[1, iexp] = datamodels.open(cal_output_files_2[iexp])
        flat_data_2 = datamodels.open(intflat_output_files_2[iexp])
        # Merge the FS slits into the MOS slits for NRS2
        if merge_fs_nrs2:
            msgs.info('Appending interpolatedflat FS slits into MOS output for {:s}'.format(basenames_2[iexp]))
            flat_data_fs_2 = datamodels.open(intflat_fs_output_files_2[iexp])
            for slit in flat_data_fs_2.slits:
                flat_data_2.slits.append(slit)
        flat_data[1, iexp] = flat_data_2


    # Use the first exposure to se the slit names
    # (ndet, nslit)
    slit_names_1 = [slit.name for slit in cal_data[0,0].slits]
    slit_names_2 = [slit.name for slit in cal_data[1,0].slits]
    slit_names_tot = np.hstack([slit_names_1, slit_names_2])
    source_names_1 = [slit.source_name for slit in cal_data[0,0].slits]
    source_names_2 = [slit.source_name for slit in cal_data[1,0].slits]
    source_names_tot = np.hstack([source_names_1, source_names_2])

    # Find the unique slit names and the unique sources aligned with those slits
    slit_names_uni, uni_indx = np.unique(slit_names_tot, return_index=True)
    source_names_uni = source_names_tot[uni_indx]
    slit_sources_uni = [(slit, source) for slit, source in zip(slit_names_uni, source_names_uni)]

    iexp_ref = 0
    if not os.path.isdir(scipath):
        msgs.info('Creating directory for Science output: {0}'.format(scipath))


    if _reduce_slits is not None:
        gd_slits_sources = [(slt, src) 
                            for slt, src in slit_sources_uni 
                            for slit in _reduce_slits if slt.strip() == slit.strip()]
    elif _reduce_sources is not None:
        gd_slits_sources = [(slt, src) 
                            for slt, src in slit_sources_uni 
                            for source in _reduce_sources if src == source.strip()]
    else:
        gd_slits_sources = slit_sources_uni


    # Construct the offsets object
    SlitOffsets = DitherOffsets(cal_data, gd_slits_sources)

    # Loop over all slits. For each exposure create a mosaic and save them to individual PypeIt spec2d files.
    for ii, (islit, isource) in enumerate(gd_slits_sources):

        # Clear the ginga canvas for each new source/slit
        if show:
            display.clear_all(allow_new=True)

        # TODO This step is only performed with a reference exposure because calwebb has an annoying property that
        # it does not always extract the same subimage spectral pixels for the different dithers in the dither pattern.
        # This seems to be a bug in calwebb, since it is unclear why the subimage calibrations should change.
        # This is problem for 2d coadding, since then the offsets in the detector frame will be not allow one to register
        # the frames. It is possible to fix this by using the RA/DEC images provided by calwebb to determine the
        # actual locations on the sky, which would be preferable. However, this does not appear to be working correctly
        # in calwebb. So for now, we just use the first exposure as the reference exposure for the calibrations.
        CalibrationsNRS1 = NIRSpecSlitCalibrations(det_container_list[0], cal_data[0, iexp_ref], flat_data[0, iexp_ref],
                                                islit, f070_f100_rescale=run_bogus_f100lp)
        CalibrationsNRS2 = NIRSpecSlitCalibrations(det_container_list[1], cal_data[1, iexp_ref], flat_data[1, iexp_ref],
                                                islit, f070_f100_rescale=run_bogus_f100lp)

        for iexp in range(nexp):
            # Container for all the Spec2DObj, different spec2dobj and specobjs for each slit
            all_spec2d = spec2dobj.AllSpec2DObj()
            all_spec2d['meta']['bkg_redux'] = bkg_redux
            all_spec2d['meta']['find_negative'] = bkg_redux
            # Container for the specobjs
            all_specobjs = specobjs.SpecObjs()

            # Create the image mosaic
            sciImg, slits, waveimg, tilts, ndet = jwst_mosaic(msa_data[:, iexp], [CalibrationsNRS1, CalibrationsNRS2], kludge_err=kludge_err,
                noise_floor=par['scienceframe']['process']['noise_floor'], show=show)

            # If this is a bkg_redux, perform background subtraction
            if bkg_redux:
                bkgImg_list = []
                for iexp_bkg in bkg_indices[iexp]:
                    bkgImg_i, _, _, _, _= jwst_mosaic(msa_data[:, iexp_bkg], [CalibrationsNRS1, CalibrationsNRS2], kludge_err=kludge_err,
                                            noise_floor=par['scienceframe']['process']['noise_floor'])
                    bkgImg_list.append(bkgImg_i)

                # TODO the parset label here may change in Pypeit to bkgframe
                combineImage = combineimage.CombineImage(bkgImg_list, par['scienceframe']['process'])
                bkgImg = combineImage.run(ignore_saturation=True)
                sciImg = sciImg.sub(bkgImg)

            # Run the reduction
            all_spec2d[sciImg.detector.name], tmp_sobjs = jwst_reduce(sciImg, slits, waveimg, tilts, spectrograph, par,
                                                                    show=show, find_negative=bkg_redux, bkg_redux=bkg_redux,
                                                                    clear_ginga=False, show_peaks=show, show_skysub_fit=show,
                                                                    basename=basenames[iexp])
            # Hold em
            if tmp_sobjs.nobj > 0:
                all_specobjs.add_sobj(tmp_sobjs)

                if show:
                    # Plot boxcar
                    wv_gpm_box = tmp_sobjs[0].BOX_WAVE > 1.0
                    gpm_temp = tmp_sobjs[0].BOX_MASK[wv_gpm_box]
                    flux_temp = tmp_sobjs[0].BOX_COUNTS[wv_gpm_box] * gpm_temp

                    data = np.ma.MaskedArray(flux_temp, mask=np.logical_not(gpm_temp))
                    sigclip = stats.SigmaClip(sigma=5.0, maxiters=15, cenfunc='median', stdfunc=nan_mad_std)
                    data_clipped, lower, upper = sigclip(data, masked=True, return_bounds=True)
                    gpm_clip = np.logical_not(data_clipped.mask)  # mask_stack = True are good values
                    flux_sm = fast_running_median(flux_temp*gpm_clip, 10)
                    sigma_sm = fast_running_median(tmp_sobjs[0].BOX_COUNTS_SIG[wv_gpm_box], 100)
                    ymax = 1.2*np.max(flux_sm)
                    ymin = -np.max(sigma_sm)
                    plt.plot(tmp_sobjs[0].BOX_WAVE[wv_gpm_box],tmp_sobjs[0].BOX_COUNTS[wv_gpm_box] * tmp_sobjs[0].BOX_MASK[wv_gpm_box],
                    color='green', drawstyle='steps-mid', label='Boxcar Counts')
                    plt.plot(tmp_sobjs[0].BOX_WAVE[wv_gpm_box], tmp_sobjs[0].BOX_COUNTS_SIG[wv_gpm_box] * tmp_sobjs[0].BOX_MASK[wv_gpm_box],
                        color='cyan', drawstyle='steps-mid', label='Boxcar Counts Error')

                    # plot optimal
                    if tmp_sobjs[0].OPT_WAVE is not None:
                        wv_gpm_opt = tmp_sobjs[0].OPT_WAVE > 1.0
                        plt.plot(tmp_sobjs[0].OPT_WAVE[wv_gpm_opt], tmp_sobjs[0].OPT_COUNTS[wv_gpm_opt] * tmp_sobjs[0].OPT_MASK[wv_gpm_opt],
                                color='black', drawstyle='steps-mid', label='Optimal Counts')
                        plt.plot(tmp_sobjs[0].OPT_WAVE[wv_gpm_opt], tmp_sobjs[0].OPT_COUNTS_SIG[wv_gpm_opt] * tmp_sobjs[0].OPT_MASK[wv_gpm_opt],
                                color='red', drawstyle='steps-mid', label='Optimal Counts Error')

                    plt.ylim([ymin, ymax])
                    plt.legend()
                    plt.show()

            # THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()
            base_suffix = 'slit_{:s}'.format(islit) if _reduce_slits is not None else 'source_{:s}'.format(isource)
            #base_suffix = 'source_{:s}'.format(isource) if isource is not None else 'slit_{:s}'.format(islit)
            basename = '{:s}_{:s}'.format(basenames[iexp], base_suffix)

            # TODO Populate the header with metadata relevant to this source?

            # Write out specobjs
            # Build header for spec2d
            head2d = fits.getheader(rate_files_1[iexp])
            subheader = spectrograph.subheader_for_spec(fitstbl_1[iexp], head2d, allow_missing=False)
            # Overload the target name with the source name
            subheader['target'] = isource
            if all_specobjs.nobj > 0:
                outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
                all_specobjs.write_to_fits(subheader, outfile1d)

            # Info
            outfiletxt = os.path.join(scipath, 'spec1d_{:s}.txt'.format(basename))
            all_specobjs.write_info(outfiletxt, spectrograph.pypeline)

            # Build header for spec2d
            outfile2d = os.path.join(scipath, 'spec2d_{:s}.fits'.format(basename))
            # TODO For the moment hack so that we can write this out
            pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph, subheader=subheader,
                                                redux_path=None, calib_dir=None)
            # Write spec2d
            all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr, overwrite=True)

    return SlitOffsets

#slits = 'S200A2' #'S200A1' #'S200A2' #'S200A1' # 'S200A2'
#disperser = '235H' #'235H' #'235H_bogus_FS' #'140H_bogus_FS_F100LP'
#exp_list, redux_dir = jwst_targets(progid, disperser, target, slits=slits)

# Redux parameters 
#mode = 'MSA'
#mode = 'FS'
# This only impacts the srctype
#bkg_redux = False #False #False #False #False
#run_stage1 = True
#run_stage2 = True
#overwrite_stage1 = True
#overwrite_stage2 = True 
#show=False
#reduce_slits = [slits] # None
#reduce_sources = None

#kludge_err = 1.5 # More recent reductions suggest this number ought to be more like 1.2
# Is this correct?

# if __name__ == '__main__':
    
        
#     #progid = '1764'
#     # J0313-1806
#     progid = '1764'
#     target = 'J0313-1806'
#     slits = ['S200A1', 'S200A2']
#     dispersers = ['140H', '235H', '395H']
#     # Slit S200A1
#     #target, slits, disperser = 'J0313-1806', 'S200A1', '140H'
#     #target, slits, disperser = 'J0313-1806', 'S200A1', '235H'
#     #target, slits, disperser = 'J0313-1806', 'S200A1', '395H'
#     # Slit S200A2
#     #target, slits, disperser = 'J0313-1806', 'S200A2', '140H'
#     #target, slits, disperser = 'J0313-1806', 'S200A2', '235H'
#     #target, slits, disperser = 'J0313-1806', 'S200A2', '395H'

#     for disp in dispersers:
#         for slt in slits:
#             print('Running NIRSpec redux for program {0}, target {1}, slit {2}, disperser {3}'.format(progid, target, slt, disp))            
#             #jwst_mosaic(progid, disp, target, slits=slt, show=False, overwrite_stage1=False, overwrite_stage2=False, reduce_slits=None, reduce_sources=None, kludge_err=1.5, bkg_redux=False)