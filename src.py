import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches

# from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

from astropy.coordinates import SkyCoord
import astropy.units as u 
import astropy.io 





# Load a FITS file (include the data as a 2d array, the header, and the full path to the file that was loaded)
def load_fits(fullpath):
    print(f"Loading {fullpath}") 
    data= astropy.io.fits.getdata(fullpath, ext=0)  
    header = astropy.io.fits.open(fullpath)[0].header 
    file = {"data": data, "header": header, "fullpath": fullpath}
    return file





def plot_fits( 
        ax, 
        file, 
        vrange_std=None,
        vrange_absval=None,
        cmap='gray', 
    ):

    image = file["data"] 

    # Choose vmin and vmax 
    if vrange_std is None and vrange_absval is None:
        vrange_std = (-1, 4)
    if vrange_std is not None and vrange_absval is not None:
        raise ValueError("Provide either vrange_std or vrange_absval, not both")
    mean = np.mean(image)
    std = np.std(image)
    if vrange_std is not None:
        vmin = mean + min(vrange_std) * std
        vmax = mean + max(vrange_std) * std
    else:
        vmin, vmax = min(vrange_absval), max(vrange_absval)

    # Display image 
    ax.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax)




# Retrive RA, dec, and mag of all stars within a certain radius of a certain RA, dec coordinate 
def retrieve_stars(
        center_ra, 
        center_dec, 
        radius_deg, 
        catalog="I/239/hip_main", 
        mag_limit=18, 
    ): 
    """
    Query a star catalog for objects within a circular field of view and
    return their sky coordinates and magnitudes.

    This function uses VizieR to search a specified astronomical catalog
    (e.g., Gaia DR3 or Hipparcos) for all stars within a given angular
    radius of a central sky position. A magnitude cut is applied at the
    catalog-query level to improve performance by limiting faint sources.

    Parameters
    ----------
    center_ra : float
        Right ascension of the field center in degrees.

    center_dec : float
        Declination of the field center in degrees.

    radius_deg : float
        Search radius around the field center in degrees.

    catalog : str, optional
        VizieR catalog identifier to query. Supported examples:
        - "I/355/gaiadr3" (Gaia DR3)
        - "I/239/hip_main" (Hipparcos main catalog)

    mag_limit : float, optional
        Upper magnitude limit used in the query (fainter stars are excluded
        at the server level). The relevant magnitude band depends on the
        catalog (e.g., Gmag for Gaia, Vmag for Hipparcos).

    Returns
    -------
    dict
        Dictionary containing:
        - "ra" : numpy.ndarray
            Right ascension values in degrees.
        - "dec" : numpy.ndarray
            Declination values in degrees.
        - "mag" : numpy.ndarray
            Apparent magnitudes in the catalog's native band.

        If no stars are found, returns empty arrays.
    """
    
    print(f"Searching '{catalog}'...")

    # Choose keywords to access the variables depending on which catalog is chosen  
    if catalog == "I/355/gaiadr3": 
        filter_str = "Gmag" 
        ra_str = "RA_ICRS" 
        dec_str = "DE_ICRS" 
    elif catalog == "I/239/hip_main": 
        filter_str = "Vmag" 
        ra_str = "RAICRS" 
        dec_str = "DEICRS"

    # Apply magnitude limit to speed up search 
    vizier = Vizier(
        columns=[ra_str, dec_str, filter_str],
        column_filters={filter_str: f"<{mag_limit}"}
        )
    vizier.ROW_LIMIT = -1 

    # Search 
    coord = SkyCoord(
        ra=center_ra*u.deg,
        dec=center_dec*u.deg, 
    )
    tables = vizier.query_region(
        coord,
        radius=radius_deg*u.deg,
        catalog=catalog
    )

    # Handle case where no stars are found 
    try: 
        stars = tables[0]
    except: 
        result = {"ra": [], "dec": [], "mag": []}
        print(f"Found {0} stars with mag < {mag_limit}") 
        return result

    # Store RA, Dec, and Magnitude in a single dictionary with the same keywords 
    ra = np.array(stars[ra_str]) 
    dec = np.array(stars[dec_str])
    mag = np.array(stars[filter_str])
    result = {
        "ra": ra, 
        "dec": dec, 
        "mag": mag, 
    }
    print(f"Found {len(stars)} stars with mag < {mag_limit}") 
    return result






def convert_skycoords_to_pixelcoords(stars, center_ra, center_dec, theta_deg, camera_width_pix, camera_height_pix, plate_scale, invert=False): 
    """
    Convert celestial coordinates (RA, Dec) into pixel coordinates for a
    simulated telescope camera field.

    This function performs a small-angle tangent-plane projection around a
    specified field center, applies an optional rotation (camera orientation),
    and converts angular offsets into pixel units using the detector plate scale.

    The resulting pixel coordinates are added in-place to the input dictionary
    under the keys "x" and "y".

    Parameters
    ----------
    stars : dict
        Dictionary containing at minimum:
        - "ra" : array-like
        - "dec" : array-like

        These are assumed to be in degrees.

    center_ra : float
        Right ascension of the field center in degrees.

    center_dec : float
        Declination of the field center in degrees.

    theta_deg : float
        Rotation angle of the camera in degrees. This defines how the field is
        rotated in the image plane.

    camera_width_pix : float
        Width of the detector/image in pixels.

    camera_height_pix : float
        Height of the detector/image in pixels.

    plate_scale : float
        Plate scale of the instrument in arcseconds per pixel.

    invert : bool, optional
        If True, flips the x-axis direction (used to match different telescope
        or image parity conventions).

    Returns
    -------
    dict
        The same input dictionary with two additional keys:
        - "x" : numpy.ndarray
            Pixel x-coordinates of each star.
        - "y" : numpy.ndarray
            Pixel y-coordinates of each star.

    Notes
    -----
    - Uses the small-angle approximation:
      RA offsets are scaled by cos(dec_center) to account for convergence of meridians.
    - Rotation is applied in the image plane after projection.
    - Plate scale is assumed to be in arcseconds/pixel and is converted internally.
    - Output coordinates are centered on the middle of the detector.
    - This approximation is valid for small fields (typically < a few degrees).
    """

    theta = np.deg2rad(theta_deg) 

    dra = (stars["ra"] - center_ra + 180) % 360 - 180
    dx = (-1)**(invert+1) * dra * np.cos(np.deg2rad(center_dec))
    dy = stars["dec"] - center_dec

    x_rot =  dx*np.cos(theta) + dy*np.sin(theta)
    y_rot = -dx*np.sin(theta) + dy*np.cos(theta)

    x_center = camera_width_pix / 2
    y_center = camera_height_pix / 2

    pixels_per_degree = 1/plate_scale * 3600 # pix / arcsec * 3600 arcsec/deg 
    x = x_center + x_rot * pixels_per_degree
    y = y_center + y_rot * pixels_per_degree
    stars["x"] = x 
    stars["y"] = y 
    return stars 












def plot_stars(
        ax, 
        stars,
        camera_width_pix,
        camera_height_pix,
        plot_lim_mag=12,
        plot_margin=300
    ):
    """
    Generate a simulated star field plot in pixel coordinates for a telescope
    detector, including magnitude-based filtering, brightness scaling, and
    camera field-of-view overlays.

    This function visualizes stars projected onto a detector plane. Stars are
    filtered by an apparent magnitude limit, scaled in size according to
    relative brightness, and plotted in pixel coordinates. The camera field of
    view is shown as a red rectangle, with an optional margin around the edges.

    Parameters
    ----------
    stars : dict
        Dictionary containing at minimum:
        - "x" : array-like
            Pixel x-coordinates of stars.
        - "y" : array-like
            Pixel y-coordinates of stars.
        - "mag" : array-like
            Apparent magnitudes of stars.

    camera_width_pix : int or float
        Width of the camera sensor in pixels.

    camera_height_pix : int or float
        Height of the camera sensor in pixels.

    plot_lim_mag : float, optional
        Limiting magnitude for display. Only stars brighter than this value
        (i.e., mag < plot_lim_mag) are plotted.

    plot_margin : float, optional
        Extra padding (in pixels) added around the camera field when setting
        plot limits.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created Matplotlib figure object.

    ax : matplotlib.axes.Axes
        The axes object containing the plot.

    Notes
    -----
    - Star marker size scales with a square-root function of relative
      brightness derived from magnitude.
    - The camera field of view is drawn as a red rectangle from (0,0) to
      (camera_width_pix, camera_height_pix).
    - A red cross marks the geometric center of the detector.
    - Intended for small-field astronomical simulations where distortion is
      negligible.
    - Assumes pixel coordinates have already been computed externally.
    """
    
    # Apply  mask based on slider value
    mask = stars["mag"] < plot_lim_mag
    stars_new = {key: stars[key][mask] for key in stars}

    # Make brighter stars larger 
    try: 
        brightness = (plot_lim_mag - stars_new["mag"]) / (plot_lim_mag - np.min(stars_new["mag"]))
    except: 
        brightness = [] 
    brightness = np.clip(brightness, 0, 1) 
    max_dot_size = 10 
    size = max_dot_size * brightness**0.5

    # Create plot 
    ax.scatter(stars_new["x"], stars_new["y"], s=size, color="red", zorder=5)

    # Label axes and title 
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    ax.set_title(f"Stars with mag < {plot_lim_mag:.2f} (N={mask.sum()})") 

    # Limits: add margin beyond FOV of camera 
    ax.set_xlim(0-plot_margin, camera_width_pix+plot_margin)
    ax.set_ylim(0-plot_margin, camera_height_pix+plot_margin)

    # Tick labels: evenly spaced
    xticks = np.linspace(0, camera_width_pix, 7)
    ax.set_xticks(xticks.astype(int))
    yticks = np.linspace(0, camera_height_pix, 5)
    ax.set_yticks(yticks.astype(int))

    # Add red rectangle to show the camera FOV 
    rect = patches.Rectangle(
        (0, 0),  # bottom-left corner
        camera_width_pix,
        camera_height_pix,
        linewidth=1, 
        zorder=4, 
        edgecolor="red",
        facecolor="none"
    )
    ax.add_patch(rect)

    # # Add red X crosshair at center of image 
    # plt.scatter(
    #     camera_width_pix/2,
    #     camera_height_pix/2,
    #     marker="x",
    #     color="red",
    #     s=50, 
    #     alpha=0.5, 
    # )



