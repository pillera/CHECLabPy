import warnings
import numpy as np


def get_clp_mapping_from_tc_mapping(tc_mapping):
    """
    Get a CHECLabPy mapping dataframe from a TargetCalib Mapping class,
    along with the metadata

    Parameters
    ----------
    tc_mapping : `target_calib.Mapping`

    Returns
    -------
    df : `pandas.DataFrame`

    """
    df = tc_mapping.as_dataframe()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        df.metadata = dict(
            cfgfile=tc_mapping.GetCfgPath(),
            is_single_module=tc_mapping.IsSingleModule(),
            n_pixels=tc_mapping.GetNPixels(),
            n_modules=tc_mapping.GetNModules(),
            n_tmpix=tc_mapping.GetNTMPix(),
            n_superpixels=tc_mapping.GetNSuperPixels(),
            n_rows=tc_mapping.GetNRows(),
            n_columns=tc_mapping.GetNColumns(),
            size=tc_mapping.GetSize(),
            fOTUpRow_l=tc_mapping.fOTUpRow_l,
            fOTUpRow_u=tc_mapping.fOTUpRow_u,
            fOTUpCol_l=tc_mapping.fOTUpCol_l,
            fOTUpCol_u=tc_mapping.fOTUpCol_u,
            fOTUpX_l=tc_mapping.fOTUpX_l,
            fOTUpX_u=tc_mapping.fOTUpX_u,
            fOTUpY_l=tc_mapping.fOTUpY_l,
            fOTUpY_u=tc_mapping.fOTUpY_u
        )
    return df


def get_clp_mapping_from_version(version, single_module=False):
    from target_calib import CameraConfiguration
    tc_mapping = CameraConfiguration(version).GetMapping(
        singleModule=single_module
    )
    return get_clp_mapping_from_tc_mapping(tc_mapping)


def get_ref_path_from_version(version):
    from target_calib import CameraConfiguration
    ref_path = CameraConfiguration(version).GetReferencePulsePath()
    return ref_path


def get_tc_mapping_from_clp_mapping(mapping):
    """
    Get a TargetCalib Mapping class from a CHECLabPy mapping dataframe

    Parameters
    ----------
    mapping : `pandas.DataFrame`
        The CHECLabPy dataframe representation of the mapping

    Returns
    -------
    tc_mapping : `target_calib.Mapping`

    """
    from target_calib import Mapping
    cfgfile = mapping.metadata.cfgfile
    is_single_module = mapping.metadata.is_single_module
    return Mapping(cfgfile, is_single_module)


def get_superpixel_mapping(mapping):
    """
    Reorganise a CHECLabPy mapping dataframe to represent the positions of the
    superpixels on the camera

    Parameters
    ----------
    mapping : `pandas.DataFrame`
        The CHECLabPy dataframe representation of the mapping

    Returns
    -------
    `pandas.DataFrame`

    """
    df = mapping[['superpixel', 'slot', 'asic', 'row', 'col', 'xpix', 'ypix']]
    f_rowcol = lambda v: v.values[0] // 2
    f = dict(slot='first', asic='first', row=f_rowcol, col=f_rowcol,
             xpix='mean', ypix='mean')
    df = df.groupby('superpixel').agg(f).reset_index()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        df.metadata = mapping.metadata.copy()
    df.metadata['n_rows'] = df['row'].max() + 1
    df.metadata['n_columns'] = df['col'].max() + 1
    df.metadata['size'] *= 2
    return df


def get_tm_mapping(mapping):
    """
    Reorganise a CHECLabPy mapping dataframe to represent the positions of the
    TMs on the camera

    Parameters
    ----------
    mapping : `pandas.DataFrame`
        The CHECLabPy dataframe representation of the mapping

    Returns
    -------
    `pandas.DataFrame`

    """
    df = mapping[['slot', 'row', 'col', 'xpix', 'ypix']]
    f_rowcol = lambda v: v.values[0] // 8
    f = dict(row=f_rowcol, col=f_rowcol, xpix='mean', ypix='mean')
    df = df.groupby('slot').agg(f, as_index=False).reset_index()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        df.metadata = mapping.metadata.copy()
    df.metadata['n_rows'] = df['row'].max() + 1
    df.metadata['n_columns'] = df['col'].max() + 1
    df.metadata['size'] *= 8
    return df


def get_ctapipe_camera_geometry(mapping, plate_scale=None):
    """
    Obtain a ctapipe CameraGeometry object from the CHECLabPy Mapping object.

    Pixel coordinates are converted into degrees using the plate scale.

    Parameters
    ----------
    mapping : `pandas.DataFrame`
        The mapping for the pixels stored in a pandas DataFrame. Can be
        obtained from either of these options:

        CHECLabPy.io.TIOReader.mapping
        CHECLabPy.io.ReaderR0.mapping
        CHECLabPy.io.ReaderR1.mapping
        CHECLabPy.io.DL1Reader.mapping
        CHECLabPy.utils.mapping.get_clp_mapping_from_tc_mapping

    Returns
    -------
    geom : `ctapipe.instrument.camera.CameraGeometry`
    """
    from ctapipe.instrument import CameraGeometry
    from astropy import units as u

    if plate_scale:
        mapping['xpix'] /= plate_scale
        mapping['ypix'] /= plate_scale
        mapping.metadata['size'] /= plate_scale

    camera = CameraGeometry(
        "CHEC",
        pix_id=np.arange(mapping.metadata['n_pixels']),
        pix_x=mapping['xpix'].values * u.m,
        pix_y=mapping['ypix'].values * u.m,
        pix_area=None,
        pix_type='rectangular',
    )
    return camera
