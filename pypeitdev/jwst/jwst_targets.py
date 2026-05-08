import os
from IPython import embed 

def jwst_spec1d_files(progid, disperser, target, slit=None, source=None):
    """
    Routine to return a list of spec1d filenames for JWST NIRSpec exposures for a given target and disperser and possibly slit.
    
    Parameters
    ----------
    progid : str
        Program ID for the observations. 
    disperser : str
        Name of the disperser.
    target : str
        Name of the target.
    slit : str
        Slit requested, optional default is None.
    source : str
        Source name, optional default is None.
    
    Return 
    ------
    spec1d_list : list
        List of spec1d filenames.
    """

    if isinstance(slit, (tuple, list)) and len(slit) == 1:
        slit = slit[0]
    if isinstance(source, (tuple, list)) and len(source) == 1:
        source = source[0]

    if slit is None and source is None:
        raise ValueError("Either 'slit' or 'source' must be specified.")
    elif slit is not None and source is not None:
        raise ValueError("Only one of 'slit' or 'source' can be specified.")


    uncal_list, redux_path, rawpath_level2 = jwst_targets(progid, disperser, target, slit=slit)
    suffix = f'_slit_{slit}.fits' if slit is not None else f'_source_{source}.fits'
    spec1d_pre = ['spec1d_' + os.path.basename(file).split('_nrs1')[0] + suffix for file in uncal_list[0]]
    spec1d_filenames = [os.path.join(redux_path, 'pypeit', 'Science', spec) for spec in spec1d_pre]
    
    return redux_path, spec1d_filenames


def jwst_targets(progid, disperser, target, slit=None):
    """
    Routine to return a list of JWST NIRSpec exposures for a given target and disperser and possibly slit.
    
    Parameters
    ----------
    progid : str
        Program ID for the observations. 
    disperser : str
        Name of the disperser.
    target : str
        Name of the target.
    slit : str
        Slit requested, optional default is None.
        
    Returns
    -------
    exp_list : list
        List of lists of uncalibrated files for each exposure. exp_list[0] contains 
        all of the exposure files for nrs1 and exp_list[1] contains all of the exposure files for nrs2.
    redux_dir : str
        Path to the directory where the data will be reduced.
    rawpath_level2 : str
        Path to the directory where the raw level 2 output files are.
    
    """

    exp_list = []
    detectors = ['nrs1', 'nrs2']

    # If bkg_redux is False, the code will model the sky and the object profile and perform optimal extraction.
    # If bkg_redux is True, the code will difference image and simply boxcar extract (optimal not implemented yet)
    for detname in detectors:

        if '3543' in progid:
            if 'G395M' == disperser:
                ## BHstar object
                if target == 'BHstar':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO3543/G395M/stage1_rate'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/3543/G395M/', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw03543001001_07101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw03543001001_07101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw03543001001_07101_00004_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3]) 
        
        if '1181' in progid:
            if target == 'JADES-GN-28074':

                base_raw = "/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO1181"
                base_redux = "/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/1181"

                if disperser == 'G140M':
                    subdir = 'G140M'
                    visit_id = '04101'
                    #visit_id = '18101'

                elif disperser == 'G235M':
                    subdir = 'G235M'
                    #visit_id = '06101'
                    visit_id = '20101'

                elif disperser == 'G395M':
                    subdir = 'G395M'
                    visit_id = '08101'

                elif disperser == 'PRISM':
                    subdir = 'PRISM'
                    #visit_id = '10101'
                    #visit_id = '16101'
                    #visit_id = '12101'
                    visit_id = '14101'

                else:
                    raise ValueError(f"Disperser '{disperser}' not recognized")

                rawpath_level2 = os.path.join(base_raw, subdir)
                redux_dir = os.path.join(base_redux, subdir, target)

                # Construct uncal file paths
                uncalfile1 = os.path.join(rawpath_level2, f'jw01181004001_{visit_id}_00001_{detname}_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, f'jw01181004001_{visit_id}_00002_{detname}_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, f'jw01181004001_{visit_id}_00003_{detname}_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        
            elif target == 'JADES-GN-72127':
                base_raw = "/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO1181"
                base_redux = "/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/1181"

                if disperser == 'G140M':
                    subdir = 'G140M'
                    visit_id = '09101'

                elif disperser == 'G235M':
                    subdir = 'G235M'
                    #visit_id = '17101'
                    visit_id = '07101'

                elif disperser == 'G395M':
                    subdir = 'G395M'
                    visit_id = '25101'

                elif disperser == 'G395H':
                    subdir = 'G395H'
                    visit_id = '03101'

                elif disperser == 'PRISM':
                    subdir = 'PRISM'
                    visit_id = '11101'
                    #visit_id = '13101'
                    #visit_id = '31101'

                else:
                    raise ValueError(f"Disperser '{disperser}' not recognized")

                rawpath_level2 = os.path.join(base_raw, subdir)
                redux_dir = os.path.join(base_redux, subdir, target)

                uncalfile1 = os.path.join(rawpath_level2, f'jw01181009001_{visit_id}_00001_{detname}_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, f'jw01181009001_{visit_id}_00002_{detname}_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, f'jw01181009001_{visit_id}_00003_{detname}_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        # RUBIES LRD program
        if '4233' in progid:
            if disperser == 'PRISM':
                if target == 'RUBIES-UDS-154183':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/4233/PRISM'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/PRISM', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233003001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233003001_03101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233003001_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw04233003002_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw04233003002_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw04233003003_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile7 = os.path.join(rawpath_level2, 'jw04233003003_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile8 = os.path.join(rawpath_level2, 'jw04233003003_03101_00004_' + detname + '_uncal.fits')

                    exp_list.append([uncalfile1, uncalfile2, uncalfile3])

                if target == 'RUBIES-EGS-42046':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/4233/PRISM'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/PRISM', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233005001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233005001_03101_00004_' + detname + '_uncal.fits')

                    exp_list.append([uncalfile1, uncalfile2])

                if target == 'RUBIES-EGS-40579':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO4233/40579/PRISM'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/PRISM', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233001001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233001001_03101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233001001_03101_00004_' + detname + '_uncal.fits')

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw04233001002_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw04233001002_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw04233001002_03101_00004_' + detname + '_uncal.fits')

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw04233001003_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw04233001003_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw04233001003_03101_00004_' + detname + '_uncal.fits')


                    exp_list.append([uncalfile1, uncalfile2, uncalfile3])


                if target == 'RUBIES-UDS-QG-z7':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO4233/977881/PRISM'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/PRISM', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233003002_03101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233003002_03101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233003002_03101_00004_' + detname + '_uncal.fits')


                    exp_list.append([uncalfile1, uncalfile2, uncalfile3])

            elif disperser == 'G395M':
                if target == 'RUBIES-UDS-154183':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/4233/G395M'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/G395M', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233003001_05101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233003001_05101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233003001_05101_00004_' + detname + '_uncal.fits')

                    exp_list.append([uncalfile1, uncalfile2, uncalfile3])

                if target == 'RUBIES-EGS-42046':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/4233/G395M'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/G395M', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233005001_05101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233005001_05101_00004_' + detname + '_uncal.fits')

                    exp_list.append([uncalfile1, uncalfile2])

                if target == 'RUBIES-EGS-49140':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/4233/G395M'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/G395M', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233006001_05101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233006001_05101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233006002_05101_00002_' + detname + '_uncal.fits')
                    uncalfile4 = os.path.join(rawpath_level2, 'jw04233006002_05101_00004_' + detname + '_uncal.fits')
                    uncalfile5 = os.path.join(rawpath_level2, 'jw04233006003_05101_00002_' + detname + '_uncal.fits')
                    uncalfile6 = os.path.join(rawpath_level2, 'jw04233006003_05101_00004_' + detname + '_uncal.fits')

                    exp_list.append([uncalfile1, uncalfile2, uncalfile3, 
                                    uncalfile4, uncalfile5, uncalfile6])
                    

                if target == 'RUBIES-UDS-QG-z7':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO4233/977881/G395M'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/4233/G395M', target)

                    uncalfile1 = os.path.join(rawpath_level2, 'jw04233003001_05101_00002_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw04233003001_05101_00003_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw04233003001_05101_00004_' + detname + '_uncal.fits')


                    exp_list.append([uncalfile1, uncalfile2, uncalfile3])


        # UNCOVER LRD
        if '2561' in progid:
            if 'PRISM' == disperser:
                
                if target == 'A2744-45924':
                    rawpath_level2 = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/2561/PRISM'
                    redux_dir = os.path.join('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/2561/PRISM', target)

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002001_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002001_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002001_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002001_03101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002001_03101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002001_03101_00007_' + detname + '_uncal.fits')

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002001_05101_00001_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002001_05101_00002_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002001_05101_00003_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002001_05101_00004_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002001_05101_00005_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002001_05101_00006_' + detname + '_uncal.fits')

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002002_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002002_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002002_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002002_03101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002002_03101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002002_03101_00007_' + detname + '_uncal.fits')

                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002003_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002003_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002003_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002003_03101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002003_03101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002003_03101_00007_' + detname + '_uncal.fits')

                    # THIS WORKS FOR A2744-45924, need to check others
                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002004_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002004_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002004_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002004_03101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002004_03101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002004_03101_00007_' + detname + '_uncal.fits')

                    # THIS WORKS FOR A2744-45924, need to check others
                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002005_03101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002005_03101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002005_03101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002005_03101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002005_03101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002005_03101_00007_' + detname + '_uncal.fits')

                    uncalfile1 = os.path.join(rawpath_level2, 'jw02561002005_11101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02561002005_11101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02561002005_11101_00003_' + detname + '_uncal.fits')
                    uncalfile4 = os.path.join(rawpath_level2, 'jw02561002005_11101_00004_' + detname + '_uncal.fits')
                    uncalfile5 = os.path.join(rawpath_level2, 'jw02561002005_11101_00005_' + detname + '_uncal.fits')
                    uncalfile6 = os.path.join(rawpath_level2, 'jw02561002005_11101_00006_' + detname + '_uncal.fits')

                    # THIS WORKS FOR A2744-45924, need to check others
                    # uncalfile1 = os.path.join(rawpath_level2, 'jw02561002005_15101_00002_' + detname + '_uncal.fits')
                    # uncalfile2 = os.path.join(rawpath_level2, 'jw02561002005_15101_00003_' + detname + '_uncal.fits')
                    # uncalfile3 = os.path.join(rawpath_level2, 'jw02561002005_15101_00004_' + detname + '_uncal.fits')
                    # uncalfile4 = os.path.join(rawpath_level2, 'jw02561002005_15101_00005_' + detname + '_uncal.fits')
                    # uncalfile5 = os.path.join(rawpath_level2, 'jw02561002005_15101_00006_' + detname + '_uncal.fits')
                    # uncalfile6 = os.path.join(rawpath_level2, 'jw02561002005_15101_00007_' + detname + '_uncal.fits')

                #exp_list.append([uncalfile1, uncalfile2, uncalfile3, uncalfile4, uncalfile5, uncalfile6])
                #exp_list.append([uncalfile1, uncalfile3, uncalfile5])
                #exp_list.append([uncalfile2, uncalfile4, uncalfile6])
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        # SPURS - Abell2744-QSO1 (image B, Furtak et al. 2024); Cycle 4, PID 9214.
        # Three medium gratings: G140M/F100LP, G235M/F170LP, G395M/F290LP.
        # Data staged by mgii_forest_jwst.utils.stage_msa_for_target into
        #   /Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO9214/Abell2744-QSO1/<disperser>/
        # MSA source table source_id = 41, source_name = 9214_41.
        if '9214' in progid:
            if target == 'Abell2744-QSO1':
                base_raw = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO9214/Abell2744-QSO1'
                base_redux = '/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/9214'

                # (obs_visit, activity) pairs that actually cover image B per the
                # stage_9214_qso1.py matcher. In program 9214, obs=009 carries
                # image B; within obs=009 the activities below hold each grating.
                # Each entry in visit_blocks is one 3-dither activity block.
                # mgii_forest_jwst.redux.redux_utils.calwebb_pypeit chunks the
                # exposure list by 3 and runs jwst_run_redux per chunk, so the
                # hardcoded bkg_indices=[(1,2),(0,2),(0,1)] nod pattern still
                # applies to each block.
                if disperser == 'G140M':
                    subdir = 'G140M'
                    visit_blocks = [
                        ('09214009001', '03101'),
                        ('09214009001', '03103'),
                        ('09214009001', '03105'),
                        ('09214009002', '03101'),
                    ]
                elif disperser == 'G235M':
                    subdir = 'G235M'
                    visit_blocks = [
                        ('09214009002', '03103'),
                    ]
                elif disperser == 'G395M':
                    subdir = 'G395M'
                    visit_blocks = [
                        ('09214009002', '03105'),
                    ]
                else:
                    raise ValueError(f"Disperser '{disperser}' not supported for Abell2744-QSO1 in GO9214")

                rawpath_level2 = os.path.join(base_raw, subdir)
                redux_dir = os.path.join(base_redux, subdir, target)
                exp_nums = ['00001', '00002', '00003']

                uncalfiles = []
                for obs_visit, act_id in visit_blocks:
                    for en in exp_nums:
                        uncalfiles.append(os.path.join(
                            rawpath_level2,
                            f'jw{obs_visit}_{act_id}_{en}_{detname}_uncal.fits'))
                exp_list.append(uncalfiles)

            # Program-level "all sources" reductions (one shared redux dir per
            # pointing × grating; rate files come from the GO9214 pool dirs
            # that stage_9214_qso1.flatten_into_pool() already populates).
            #
            # `_pypdef` suffix: same input cal files / visit blocks, but
            # outputs land in 9214_pypdef/... so a stock-PypeIt-defaults
            # extraction can run side-by-side with the standard 9214/... tree.
            _go9214_pointings = ('GO9214_obs008', 'GO9214_obs009',
                                  'GO9214_obs008_pypdef', 'GO9214_obs009_pypdef')
            if target in _go9214_pointings:
                obs = '008' if 'obs008' in target else '009'
                is_pypdef = target.endswith('_pypdef')
                base_raw = '/Users/jiamuh/jwst_redux/Raw/NIRSPEC_MSA/GO9214/pool'
                base_redux = ('/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/9214_pypdef'
                               if is_pypdef
                               else '/Users/jiamuh/jwst_redux/redux/NIRSPEC_MSA/9214')

                # Each entry is (jw_obs_visit_prefix, 5-digit activity id).
                # G140M is identical between obs 008 / 009 (both use 03101/3/5
                # in vis001 and 03101 in vis002), but G235M and G395M activity
                # IDs differ between the two pointings (vis002 only).
                if disperser == 'G140M':
                    subdir = 'G140M'
                    visit_blocks = [
                        (f'09214{obs}001', '03101'),
                        (f'09214{obs}001', '03103'),
                        (f'09214{obs}001', '03105'),
                        (f'09214{obs}002', '03101'),
                    ]
                elif disperser == 'G235M':
                    subdir = 'G235M'
                    g235m_act = '05101' if obs == '008' else '03103'
                    visit_blocks = [(f'09214{obs}002', g235m_act)]
                elif disperser == 'G395M':
                    subdir = 'G395M'
                    g395m_act = '07101' if obs == '008' else '03105'
                    visit_blocks = [(f'09214{obs}002', g395m_act)]
                else:
                    raise ValueError(f"Disperser '{disperser}' not supported for {target}")

                rawpath_level2 = os.path.join(base_raw, subdir)
                redux_dir = os.path.join(base_redux, subdir, f'all_obs{obs}')
                exp_nums = ['00001', '00002', '00003']

                uncalfiles = []
                for obs_visit, act_id in visit_blocks:
                    for en in exp_nums:
                        uncalfiles.append(os.path.join(
                            rawpath_level2,
                            f'jw{obs_visit}_{act_id}_{en}_{detname}_uncal.fits'))
                exp_list.append(uncalfiles)

        if '2073' in progid:
            if 'PRISM' == disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2073/level_12/02073/'
                #redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02073_CLEAR_PRISM', target)
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/2073/', target)

                #J0252
                if target == 'J0252-0503':
                    uncalfile1 = os.path.join(rawpath_level2, 'jw02073007001_03101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02073007001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02073007001_03101_00003_' + detname + '_uncal.fits')

                elif target == 'J1007+2115':
                    # J1007
                    # NIRSPEC 3-point dither
                    uncalfile1 = os.path.join(rawpath_level2, 'jw02073008001_03101_00001_' + detname + '_uncal.fits')
                    uncalfile2 = os.path.join(rawpath_level2, 'jw02073008001_03101_00002_' + detname + '_uncal.fits')
                    uncalfile3 = os.path.join(rawpath_level2, 'jw02073008001_03101_00003_' + detname + '_uncal.fits')

                    #uncalfile1 = os.path.join(rawpath_level2, 'jw02073006001_03101_00001_' + detname + '_uncal.fits')
                    #uncalfile2 = os.path.join(rawpath_level2, 'jw02073006001_03101_00002_' + detname + '_uncal.fits')
                    #uncalfile3 = os.path.join(rawpath_level2, 'jw02073006001_03101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3]) 

        elif '2756' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_2756/level_12/02756/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/02756_CLEAR_PRISM/calwebb'
                output_dir = os.path.join(redux_dir, 'output')
                pypeit_output_dir = os.path.join(redux_dir, 'pypeit')

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw02756001001_03101_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02756001001_03101_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02756001001_03101_00003_' + detname + '_uncal.fits')

                # dither offset
                # uncalfile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')
                # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '1133' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/Raw'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01133_COM_CLEAR_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')

                # dither offset
                # uncalfile  = os.path.join(rawpath_level2, 'jw01133003001_0310x_00003_' + detname + '_uncal.fits')
                # bkgfile1 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00001_' + detname + '_uncal.fits')
                # bkgfile2 = os.path.join(rawpath_level2, 'jw01133003001_0310x_00002_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        elif '2072' in progid:
            if 'PRISM' in disperser:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/2072/level_12'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/2072/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/02027_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw02072002001_05101_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02072002001_05101_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02072002001_05101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        elif '1219' in progid:
            if 'J1342+0928' in target:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_1219/1219/level_12/01219/'
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/1219/J1342+0928/')

                if disperser == '140H':
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_04101_00001_' + detname + '_uncal.fits')  # msa_metadata_id  = 1
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_04101_00002_' + detname + '_uncal.fits')  
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_04101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_10101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 29
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_10101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_10101_00003_' + detname + '_uncal.fits')
                elif disperser == '235H':
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_06101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 15
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_06101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_06101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01219006001_08101_00001_' + detname + '_uncal.fits') # msa_metadata_id  = 16
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01219006001_08101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01219006001_08101_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])           

        elif '1222' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1222/level_12/01222/'
            if 'J0411-0907' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J0411-0907/'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03108_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03106_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222002001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222002001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222002001_03104_00003_' + detname + '_uncal.fits')

            elif 'J0020-3653' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J0020-3653/'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03108_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03106_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222012001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222012001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222012001_03104_00003_' + detname + '_uncal.fits')
                        
            elif 'J1120+0641' in target:
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1222/J1120+0641/'
                if disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_03101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_03101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_03101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_09101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_09101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_09101_00003_' + detname + '_uncal.fits')
                elif disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_05101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_05101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_05101_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01222005001_07101_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01222005001_07101_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01222005001_07101_00003_' + detname + '_uncal.fits')                        

            exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        elif '4713' in progid:
            if 'J0100+2802' in target:
                ## Prorgram for Slit Loss Characterization for MSA shutters
                # PRISM data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/4713/'
                redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_MSA/4713/J0100+2802/')

                if disperser == '140M':
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw04713002001_03101_00002_' + detname + '_uncal.fits')  # msa_metadata_id  = 1
                        uncalfile2 = os.path.join(rawpath_level2, 'jw04713002001_03101_00003_' + detname + '_uncal.fits')  
                exp_list.append([uncalfile1, uncalfile2])           


        if '1764' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1764/level_12/01764/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/', target)
            if 'J0313-1806' in target:
                #redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb'
                #output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb/output'
                #pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J0313-1806/calwebb/pypeit'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03104_00003_' + detname + '_uncal.fits')
                elif disperser == '395H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03106_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_03108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_03108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_03108_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_0310a_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764014001_0310c_00003_' + detname + '_uncal.fits')
                
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

            elif 'J1007+2115' in target:
                #redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb'
                #output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb/output'
                #pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_FS/1764/J1007+2115/calwebb/pypeit'
                if disperser == '235H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04102_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04102_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04102_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04104_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04104_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04104_00003_' + detname + '_uncal.fits')
                elif disperser == '395H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04106_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04106_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04106_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_04108_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_04108_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_04108_00003_' + detname + '_uncal.fits')
                elif disperser == '140H':
                    # NIRSPEC 3-point dither dither center
                    if slit == 'S200A1':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_0410a_00003_' + detname + '_uncal.fits')
                    elif slit == 'S200A2':
                        uncalfile1 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00001_' + detname + '_uncal.fits')
                        uncalfile2 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00002_' + detname + '_uncal.fits')
                        uncalfile3 = os.path.join(rawpath_level2, 'jw01764006001_0410c_00003_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '3526' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/3526/level_12/03526/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/3526/', target)
            #output_dir = os.path.join(redux_dir, 'calwebb', 'output')        
            #pypeit_output_dir = os.path.join(redux_dir, 'calwebb', 'pypeit')
            # All of these are with slit S200A2
            if 'J0410-0139' in target:
                file_list = []
                for ii in range(1,11): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526012001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J0038-1527' in target:
                file_list = []
                for ii in range(1,4): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526007001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)  
            elif 'J0252-0503' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526008001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                 
            elif 'J0038-0653' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526001001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            # All of these are with slit S200A1
            elif 'J0313-1806' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526002001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                       
            elif 'J1007+2115' in target:
                file_list = []
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03526005001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                

        if '9180' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/9180/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/9180/', target)
            file_list = []
            if 'J2356+0017' in target:
                if disperser == '140H': 
                    prefix = 'jw09180041001_03102_000' if slit == 'S200A1' else 'jw09180041001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180041001_03107_000' if slit == 'S200A1' else 'jw09180041001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J1732+6531' in target:
                if disperser == '140H': 
                    prefix = 'jw09180012001_05101_000' if slit == 'S200A1' else 'jw09180012001_07101_000'
                    indx_range = range(1,6)                    
                elif disperser == '235H':
                    prefix = 'jw09180012001_11101_000' if slit == 'S200A1' else 'jw09180012001_09101_000'
                    indx_range = range(1,6)                    
                elif disperser == '395M':
                    prefix = 'jw09180049001_03102_000'
                    indx_range = range(1,4)
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
            if 'J1429-0104' in target:
                if disperser == '140H': 
                    prefix = 'jw09180031001_03102_000' if slit == 'S200A1' else 'jw09180031001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180031001_03108_000' if slit == 'S200A1' else 'jw09180031001_03106_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,3)           
            if 'J1428+0454' in target:        
                if disperser == '140H': 
                    prefix = 'jw09180029001_03102_000' if slit == 'S200A1' else 'jw09180029001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180029001_03107_000' if slit == 'S200A1' else 'jw09180029001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J1450-0144' in target:        
                if disperser == '140H': 
                    prefix = 'jw09180022001_03102_000' if slit == 'S200A1' else 'jw09180022001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180022001_03107_000' if slit == 'S200A1' else 'jw09180022001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)            
            if 'J1609+5328' in target:    
                if disperser == '140H': 
                    prefix = 'jw09180009001_04101_000' if slit == 'S200A1' else 'jw09180009001_06101_000'
                    indx_range = range(1,4)
                elif disperser == '235H':
                    prefix = 'jw09180009001_10101_000' if slit == 'S200A1' else 'jw09180009001_08101_000'
                    indx_range = range(1,4)                    
                elif disperser == '395M':
                    prefix = 'jw09180050001_04101_000'
                    indx_range = range(1,3)
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
            if 'J1440+0019' in target: 
                # Missing images in mast 
                if disperser == '140H': 
                    prefix = 'jw09180037001_03102_000' if slit == 'S200A1' else 'jw09180037001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180037001_03108_000' if slit == 'S200A1' else 'jw09180037001_03106_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,3)
            if 'J0443-5332' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180040001_03102_000' if slit == 'S200A1' else 'jw09180040001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180040001_03107_000' if slit == 'S200A1' else 'jw09180040001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J0412-5638' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180043001_03102_000' if slit == 'S200A1' else 'jw09180043001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180043001_03107_000' if slit == 'S200A1' else 'jw09180043001_03105_000'
                elif disperser == '395M':
                    prefix = 'jw09180055001_03102_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)            
            if 'J0522-5127' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180046001_03102_000' if slit == 'S200A1' else 'jw09180046001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180046001_03107_000' if slit == 'S200A1' else 'jw09180046001_03105_000'
                elif disperser == '395M':
                    prefix = 'jw09180058001_03102_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)            
            if 'J0446-5700' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180061001_03102_000' if slit == 'S200A1' else 'jw09180061001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180061001_03107_000' if slit == 'S200A1' else 'jw09180061001_03105_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J0451-3426' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180060001_04101_000' if slit == 'S200A1' else 'jw09180060001_06101_000'
                elif disperser == '235H':
                    prefix = 'jw09180060001_09101_000' if slit == 'S200A1' else 'jw09180060001_07101_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)
            if 'J0933+7427' in target: 
                if disperser == '140H': 
                    prefix = 'jw09180036001_03102_000' if slit == 'S200A1' else 'jw09180036001_03104_000'
                elif disperser == '235H':
                    prefix = 'jw09180036001_03107_000' if slit == 'S200A1' else 'jw09180036001_03105_000'
                elif disperser == '395M':
                    prefix = 'jw09180053001_03102_000'
                else: 
                    raise ValueError("Disperser not recognized: {}".format(disperser))
                indx_range = range(1,2)     

            else: 
                raise ValueError("Target not recognized: {}".format(target))



            for ii in indx_range: 
                file_list.append(os.path.join(rawpath_level2, prefix + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
            exp_list.append(file_list)


        if '1967' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/1967/level_12/01967/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/1967/', target)
            # All of these are with slit S200A2
            file_list = []
            if 'J2255+0251' in target:
                for ii in range(1,4): 
                    file_list.append(os.path.join(rawpath_level2, 'jw01967012001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J2236+0032' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967011001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            elif 'J0844-0052' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967002001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)            
            elif 'J0844-0132' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967003001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)   
            elif 'J0918+0139' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967005001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list) 
            elif 'J1425-0015' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967008001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                                   
            elif 'J1512+4422' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967009001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1525+4303' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967010001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1146+0124' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967006001_04102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      
            elif 'J1146-0005' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967007001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)     
            elif 'J0217-0208' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967001001_04101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)                     
            elif 'J0911+0152' in target:
                for ii in range(1,4):
                    file_list.append(os.path.join(rawpath_level2, 'jw01967004001_05101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)      


        if '3417' in progid:
            rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/3417/'
            redux_dir = os.path.join('/Users/joe/jwst_redux/redux/NIRSPEC_FS/3417/', target)
            # All of these are with slit S200A2
            file_list = []
            if 'J0844+0226' in target:
                for ii in range(1,6): 
                    file_list.append(os.path.join(rawpath_level2, 'jw03417002001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J0905+0300' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417008001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J0935-0110' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417003001_04102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J0853+0139' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417010001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)               
            if 'J1423+0206' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417006001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)               
            if 'J1254-0014' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417004001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J0207+0238' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417005001_06101_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J1416+0015' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417009001_04102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
            if 'J1423-0018' in target:
                for ii in range(1,6):
                    file_list.append(os.path.join(rawpath_level2, 'jw03417001001_03102_000' + "{:02d}".format(ii) + '_' + detname + '_uncal.fits'))
                exp_list.append(file_list)
                                                                                              

        if '1117' in progid:
            if 'PRISM' in disperser:
                # PRISM data
                rawpath_level2 = '//Users/joe/jwst_redux/Raw/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/level_12/01117'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_PRISM/01117_COM_CLEAR_PRISM/calwebb/pypeit'

                # NIRSPEC 3-point dither
                # dither center
                uncalfile1 = os.path.join(rawpath_level2, 'jw01117007001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01117007001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01117007001_03101_00004_' + detname + '_uncal.fits')

                exp_list.append([uncalfile1, uncalfile2, uncalfile3])
        if '2736' in progid:
            if 'G395M' in disperser:
                # Use islit = 37 for nrs1
                # G395M data
                rawpath_level2 = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/Raw'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/NIRSPEC_ERO/02736_ERO_SMACS0723_G395M/calwebb/pypeit'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw02736007001_03103_00001_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02736007001_03103_00002_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02736007001_03103_00003_' + detname + '_uncal.fits')
                
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

            elif 'G235M' in disperser:
                # Use islit = 38 for nrs1
                # G235M data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_ERO/02736_ERO_SMACS0723_G395MG235M/level_2/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_ERO/02736_ERO_SMACS0723_G235M/calwebb/pypeit/'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw02736007001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw02736007001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw02736007001_03101_00004_' + detname + '_uncal.fits')
            
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

        if '1671' in progid:
            if 'G395M' in disperser:
                # Use islit = 37 for nrs1
                # G395M data
                rawpath_level2 = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/Maseda/'
                redux_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb'
                output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/output'
                pypeit_output_dir = '/Users/joe/jwst_redux/redux/NIRSPEC_MSA/Maseda/395M/calwebb/pypeit'

                # NIRSPEC 3-point dither
                uncalfile1 = os.path.join(rawpath_level2, 'jw01671001001_03101_00002_' + detname + '_uncal.fits')
                uncalfile2 = os.path.join(rawpath_level2, 'jw01671001001_03101_00003_' + detname + '_uncal.fits')
                uncalfile3 = os.path.join(rawpath_level2, 'jw01671001001_03101_00004_' + detname + '_uncal.fits')
                exp_list.append([uncalfile1, uncalfile2, uncalfile3])

    return exp_list, redux_dir, rawpath_level2
    
    