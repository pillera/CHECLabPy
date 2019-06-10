import numpy as np
from astropy import units as u
from astropy.coordinates import Angle 
from CHECLabPy.utils.mapping import get_ctapipe_camera_geometry
from ctapipe.coordinates import CameraFrame, NominalFrame, HorizonFrame
from ctapipe.image.cleaning import tailcuts_clean
from ctapipe.image.muon.muon_ring_finder import ChaudhuriKunduRingFitter
from ctapipe.image.muon.features import npix_above_threshold
from ctapipe.image.muon.features import npix_composing_ring
from ctapipe.image.muon.features import ring_containment
from ctapipe.image.muon.features import ring_completeness
from ctapipe.image.muon.muon_integrator import MuonLineIntegrate

def analyze_muon_event(image,mapping,alt,azi):
    """
    Adapted for CHECLabPy from 
    ctapipe.image.muon.muon_reco_functions.analyse_muon_event
    
    Parameters:
    -----------
    image: array type event image
    mapping = reader.mapping
    alt 
    azi
    Returns:
    --------
    list with muonringparameters and muonintensityparameter 
    """
    #PARAMETERS FROM SIMTEL EVENT of ASTRI-CHEC
    
    equivalent_focal_length = 2.1500001 * u.m #teldes.optics.equivalent_focal_length
    mirror_area = 14.12623501 * u.m * u.m #teldes.optics.mirror_area    
    tailcuts = (5, 7)
    impact = (0.1, 0.95) * u.m
    ringwidth = (0.02, 0.2) * u.deg
    total_pix = 2048 
    min_pix = 164 # 8% (or 6%) as limit
    # Need to either convert from the pixel area in m^2 or check the camera specs
    ang_pixel_width = 0.163 * u.deg
    # Found from TDRs (or the pixel area)
    hole_rad = 0.171 * u.m  # Assuming approximately spherical hole
    cam_rad =  2.86 * u.deg
    # Above found from the field of view calculation
    sec_rad =  1.8 * u.m   
    # Added cleaning here. All these options should go to an input card
    cleaning = True

    muonringout = None
    muonintensityout = None
    
    geom = get_ctapipe_camera_geometry(mapping)
    x, y = geom.pix_x, geom.pix_y

    camera_coord = CameraFrame(
            x=x, y=y,
            focal_length=equivalent_focal_length,
            rotation=geom.pix_rotation
        )

    # TODO: correct this hack for values over 90
    altval = alt * u.deg #event.mcheader.run_array_direction[1]
    #azi = event.mcheader.run_array_direction[0]
    azi = azi * u.deg
    if Angle(altval) > Angle(90*u.deg):
        altval = Angle(90*u.deg)

    altaz = HorizonFrame(alt=altval,
                         az=azi)
    nom_coord = camera_coord.transform_to(
        NominalFrame(array_direction=altaz, pointing_direction=altaz)
    )
    x = nom_coord.x.to(u.deg)
    y = nom_coord.y.to(u.deg)


    clean_mask = tailcuts_clean(geom, image, picture_thresh=tailcuts[0],
                                    boundary_thresh=tailcuts[1])

    if cleaning:
        image = image * clean_mask

    if not sum(image):  # Nothing left after tail cuts
        return None

    muonring = ChaudhuriKunduRingFitter(None)
    muonringparam = muonring.fit(x, y, image)
    dist = np.sqrt(np.power(x - muonringparam. ring_center_x, 2)
                       + np.power(y - muonringparam.ring_center_y, 2))
    ring_dist = np.abs(dist - muonringparam.ring_radius)

    muonringparam = muonring.fit(
            x, y, img * (ring_dist < muonringparam.ring_radius * 0.4)
        )

    muonringparam.tel_id = {1}
    muonringparam.obs_id = None # where to get it in DL1.h5?
    muonringparam.event_id = None #give the wole reader? or design an event container

    dist_mask = np.abs(dist - muonringparam.
                           ring_radius) < muonringparam.ring_radius * 0.4
    pix_im = image * dist_mask
    nom_dist = np.sqrt(np.power(muonringparam.ring_center_x,
                                    2) + np.power(muonringparam.ring_center_y, 2))

    mir_rad = np.sqrt(mirror_area / np.pi)

    if(npix_above_threshold(pix_im, tailcuts[0]) > 0.1 * minpix
       and npix_composing_ring(pix_im) > minpix
       and nom_dist < cam_rad
       and muonringparam.ring_radius < 1.5 * u.deg
       and muonringparam.ring_radius > 1. * u.deg):
        muonringparam.ring_containment = ring_containment(
            muonringparam.ring_radius,
            cam_rad,
            muonringparam.ring_center_x,
            muonringparam.ring_center_y)


        muonringout = muonringparam
        

        ctel = MuonLineIntegrate(
                mir_rad, hole_radius=hole_rad,
                pixel_width=ang_pixel_width,
                sct_flag=True,
                secondary_radius=sec_rad
            )
        if image.shape[0] == total_pix:
            muonintensityoutput = ctel.fit_muon(muonringparam.ring_center_x,
                                                        muonringparam.ring_center_y,
                                                        muonringparam.ring_radius,
                                                        x[dist_mask], y[dist_mask],
                                                        image[dist_mask])

            muonintensityoutput.tel_id = {1}
            muonintensityoutput.obs_id = None
            muonintensityoutput.event_id = None
            muonintensityoutput.mask = dist_mask

            idx_ring = np.nonzero(pix_im)
            muonintensityoutput.ring_completeness = ring_completeness(
                        x[idx_ring], y[idx_ring], pix_im[idx_ring],
                        muonringparam.ring_radius,
                        muonringparam.ring_center_x,
                        muonringparam.ring_center_y,
                        threshold=30,
                        bins=30)
            muonintensityoutput.ring_size = np.sum(pix_im)
            dist_ringwidth_mask = np.abs(dist - muonringparam.ring_radius
                                                 ) < (muonintensityoutput.ring_width)
            pix_ringwidth_im = image * dist_ringwidth_mask
            idx_ringwidth = np.nonzero(pix_ringwidth_im)

            muonintensityoutput.ring_pix_completeness = npix_above_threshold(
                        pix_ringwidth_im[idx_ringwidth], tailcuts[0]) / len(
                        pix_im[idx_ringwidth])

            conditions = [
                muonintensityoutput.impact_parameter * u.m <
                impact[1] * mir_rad,

                muonintensityoutput.impact_parameter
                > impact[0],

                muonintensityoutput.ring_width
                < ringwidth[1],

                muonintensityoutput.ring_width
                > ringwidth[0]
            ]
            if all(conditions):
                muonintensityout = muonintensityoutput
            
    return [muonringout, muonintensityout]



