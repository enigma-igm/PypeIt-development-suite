

# This is my re-definition of the specobj class. I started from arspecobj and am making changes that I need
from __future__ import absolute_import, division, print_function

import copy
from collections import OrderedDict

import numpy as np

from astropy import units
from astropy.table import Table

from pypit import msgs
from pypit import arparse
from pypit.core import artraceslits
from pypit import ardebug as debugger

class SpecObj(object):
    """Class to handle object spectra from a single exposure
    One generates one of these Objects for each spectrum in the exposure. They are instantiated by the object
    finding routine, and then all spectral extraction information for the object are assigned as attributes

    Parameters:
    ----------
    shape: tuple (nspec, nspat)
       dimensions of the spectral image that the object is identified on
    slit_spat_pos: tuple of floats (spat_left,spat_right)
        The spatial pixel location of the left and right slit trace arrays evaluated at slit_spec_pos (see below). These
        will be in the range (0,nspat)
    slit_spec_pos: float
        The midpoint of the slit location in the spectral direction. This will typically be nspec/2, but must be in the
        range (0,nspec)

    Optional Parameters:
    -------------------
    det:   int
        Detector number. (default = 1, max = 99)
    config: str
       Instrument configuration (default = None)
    scidx: int
       Exposure index (deafult = 1, max=9999)
    objtype: str, optional
       Type of object ('unknown', 'standard', 'science')

    Attributes:
    ----------
    slitcen: float
       Center of slit in fraction of total (trimmed) detector size at ypos
    slitid: int
       Identifier for the slit (max=9999)
    objid: int
       Identifier for the object (max=999)
    """
    # Attributes
    # Init

    def __init__(self, shape, slit_spat_pos, slit_spec_pos, det = 1, config = None, slitid = None, scidx = 1, objtype='unknown'):

        #Assign from init parameters
        self.shape = shape
        self.slit_spat_pos = slit_spat_pos
        self.slit_spec_pos = slit_spec_pos
        self.config = config
        self.slitid = slitid
        self.scidx = copy.deepcopy(scidx)
        self.det = det
        self.objtype = objtype

        # ToDo add all attributes here and to the documentaiton

        # Object finding attributes
        self.objid = None
        self.idx = None
        self.spat_fracpos = None
        self.smash_peakflux = None
        self.trace_spat = None
        self.trace_spec = None
        self.fwhm = None
        self.spat_medpos = None

        # Attributes for HAND apertures, which are object added to the extraction by hand
        self.HAND_SPEC = None
        self.HAND_SPAT = None
        self.HAND_DET = None
        self.HAND_FWHM = None
        self.HAND_FLAG = False


        # Dictionaries holding boxcar and optimal extraction parameters
        self.boxcar = {}   # Boxcar extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'
        self.optimal = {}  # Optimal extraction 'wave', 'counts', 'var', 'sky', 'mask', 'flam', 'flam_var'


        # Generate IDs
        #self.slitid = int(np.round(self.slitcen*1e4))
        #self.objid = int(np.round(xobj*1e3))

        # Set index
        #self.set_idx()

        #

    def set_idx(self):
        # Generate a unique index for this exposure
        #self.idx = '{:02d}'.format(self.setup)
        self.idx = 'O{:03d}'.format(self.objid)
        self.idx += '-S{:04d}'.format(self.slitid)
        sdet = arparse.get_dnum(self.det, prefix=False)
        self.idx += '-D{:s}'.format(sdet)
        self.idx += '-I{:04d}'.format(self.scidx)

    def check_trace(self, trace, toler=1.):
        """Check that the input trace matches the defined specobjexp

        Parameters:
        ----------
        trace: ndarray
          Trace of the object
        toler: float, optional
          Tolerance for matching, in pixels
        """
        # Trace
        yidx = int(np.round(self.ypos*trace.size))
        obj_trc = trace[yidx]
        # self
        nslit = self.shape[1]*(self.xslit[1]-self.xslit[0])
        xobj_pix = self.shape[1]*self.xslit[0] + nslit*self.xobj
        # Check
        if np.abs(obj_trc-xobj_pix) < toler:
            return True
        else:
            return False

    def copy(self):
        slf = SpecObjExp(self.shape, self.config, self.scidx, self.det, self.xslit, self.ypos, self.xobj,
                       objtype=self.objtype)
        slf.boxcar = self.boxcar.copy()
        slf.optimal = self.optimal.copy()
        return slf

    # Printing
    def __repr__(self):
        # Generate sets string
        sdet = arparse.get_dnum(self.det, prefix=False)
        return ('<SpecObjExp: Setup = {:}, Slit = {:} at spec = {:7.2f} & spat = ({:7.2f},{:7.2f}) on det={:s}, scidx={:}, objid = {:} and objtype={:s}>'.format(
            self.config, self.slitid, self.slit_spec_pos, self.slit_spat_pos[0], self.slit_spat_pos[1], sdet, self.scidx, self.objid, self.objtype))


def init_exp(lordloc, rordloc, shape, maskslits,
             det, scidx, fitstbl, tracelist, settings, ypos=0.5, **kwargs):
    """ Generate a list of SpecObjExp objects for a given exposure

    Parameters
    ----------
    self
       Instrument "setup" (min=10,max=99)
    scidx : int
       Index of file
    det : int
       Detector index
    tracelist : list of dict
       Contains trace info
    ypos : float, optional [0.5]
       Row on trimmed detector (fractional) to define slit (and object)

    Returns
    -------
    specobjs : list
      List of SpecObjExp objects
    """

    # Init
    specobjs = []
    if fitstbl is None:
        fitsrow = None
    else:
        fitsrow = fitstbl[scidx]
    config = instconfig(fitsrow=fitsrow, binning=settings['detector']['binning'])
    slits = range(len(tracelist))
    gdslits = np.where(~maskslits)[0]

    # Loop on slits
    for sl in slits:
        specobjs.append([])
        # Analyze the slit?
        if sl not in gdslits:
            specobjs[sl].append(None)
            continue
        # Object traces
        if tracelist[sl]['nobj'] != 0:
            # Loop on objects
            #for qq in range(trc_img[sl]['nobj']):
            for qq in range(tracelist[sl]['traces'].shape[1]):
                slitid, slitcen, xslit = artraceslits.get_slitid(shape, lordloc, rordloc,
                                                                 sl, ypos=ypos)
                # xobj
                _, xobj = get_objid(lordloc, rordloc, sl, qq, tracelist, ypos=ypos)
                # Generate
                if tracelist[sl]['object'] is None:
                    specobj = SpecObjExp((tracelist[0]['object'].shape[:2]), config, scidx, det, xslit, ypos, xobj, **kwargs)
                else:
                    specobj = SpecObjExp((tracelist[sl]['object'].shape[:2]), config, scidx, det, xslit, ypos, xobj,
                                         **kwargs)
                # Add traces
                specobj.trace = tracelist[sl]['traces'][:, qq]
                # Append
                specobjs[sl].append(copy.deepcopy(specobj))
        else:
            msgs.warn("No objects for slit {0:d}".format(sl+1))
            specobjs[sl].append(None)
    # Return
    return specobjs


def objnm_to_dict(objnm):
    """ Convert an object name or list of them into a dict

    Parameters
    ----------
    objnm : str or list of str

    Returns
    -------
    odict : dict
      Object value or list of object values
    """
    if isinstance(objnm, list):
        tdict = {}
        for kk,iobj in enumerate(objnm):
            idict = objnm_to_dict(iobj)
            if kk == 0:
                for key in idict.keys():
                    tdict[key] = []
            # Fill
            for key in idict.keys():
                tdict[key].append(idict[key])
        # Generate the Table
        return tdict
    # Generate the dict
    prs = objnm.split('-')
    odict = {}
    for iprs in prs:
        odict[iprs[0]] = int(iprs[1:])
    # Return
    return odict


def mtch_obj_to_objects(iobj, objects, stol=50, otol=10, **kwargs):
    """
    Parameters
    ----------
    iobj : str
      Object identifier in format O###-S####-D##
    objects : list
      List of object identifiers
    stol : int
      Tolerance in slit matching
    otol : int
      Tolerance in object matching

    Returns
    -------
    matches : list
      indices of matches in objects
      None if none
    indcies : list

    """
    # Parse input object
    odict = objnm_to_dict(iobj)
    # Generate a Table of the objects
    tbl = Table(objnm_to_dict(objects))

    # Logic on object, slit and detector [ignoring sciidx for now]
    gdrow = (np.abs(tbl['O']-odict['O']) < otol) & (np.abs(tbl['S']-odict['S']) < stol) & (tbl['D'] == odict['D'])
    if np.sum(gdrow) == 0:
        return None
    else:
        return np.array(objects)[gdrow].tolist(), np.where(gdrow)[0].tolist()



def get_objid(lordloc, rordloc, islit, iobj, trc_img, ypos=0.5):
    """ Convert slit position to a slitid
    Parameters
    ----------
    det : int
    islit : int
    iobj : int
    trc_img : list of dict
    ypos : float, optional

    Returns
    -------
    objid : int
    xobj : float
    """
    yidx = int(np.round(ypos*lordloc.shape[0]))
    pixl_slit = lordloc[yidx, islit]
    pixr_slit = rordloc[yidx, islit]
    #
    xobj = (trc_img[islit]['traces'][yidx,iobj]-pixl_slit) / (pixr_slit-pixl_slit)
    objid= int(np.round(xobj*1e3))
    # Return
    return objid, xobj


def instconfig(fitsrow=None, binning=None):
    """ Returns a unique config string

    Parameters
    ----------
    fitsrow : Row
    binnings : str, optional

    Returns
    -------
    config : str
    """

    config_dict = OrderedDict()
    config_dict['S'] = 'slitwid'
    config_dict['D'] = 'dichroic'
    config_dict['G'] = 'dispname'
    config_dict['T'] = 'dispangle'
    #
    config = ''
    for key in config_dict.keys():
        try:
            comp = str(fitsrow[config_dict[key]])
        except (KeyError, TypeError):
            comp = '0'
        #
        val = ''
        for s in comp:
            if s.isdigit():
                val += s
        config = config + key+'{:s}-'.format(val)
    # Binning
    if binning is None:
        msgs.warn("Assuming 1x1 binning for your detector")
        binning = '1x1'
    val = ''
    for s in binning:
        if s.isdigit():
            val = val + s
    config += 'B{:s}'.format(val)
    # Return
    return config


def dummy_specobj(fitstbl, det=1, extraction=True):
    """ Generate dummy specobj classes
    Parameters
    ----------
    fitstbl : Table
      Expecting the fitsdict from dummy_fitsdict
    Returns
    -------

    """
    shape = fitstbl['naxis1'][0], fitstbl['naxis0'][0]
    config = 'AA'
    scidx = 5 # Could be wrong
    xslit = (0.3,0.7) # Center of the detector
    ypos = 0.5
    xobjs = [0.4, 0.6]
    specobjs = []
    for xobj in xobjs:
        specobj = SpecObjExp(shape, config, scidx, det, xslit, ypos, xobj)
        # Dummy extraction?
        if extraction:
            npix = 2001
            specobj.boxcar['wave'] = np.linspace(4000., 6000., npix)*units.AA
            specobj.boxcar['counts'] = 50.*(specobj.boxcar['wave'].value/5000.)**-1.
            specobj.boxcar['var']  = specobj.boxcar['counts'].copy()
        # Append
        specobjs.append(specobj)
    # Return
    return specobjs

#TODO We need a method to write these objects to a fits file