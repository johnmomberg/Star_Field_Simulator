import numpy as np 
import matplotlib.pyplot as plt 

from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

from astropy.coordinates import SkyCoord
import astropy.units as u 





# Retrive RA, dec, and mag of all stars within a certain radius of a certain RA, dec coordinate 
def retrieve_stars(
        center_ra, 
        center_dec, 
        radius_deg, 
        catalog="I/239/hip_main", 
        mag_limit=17, 
    ): 

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
    
    theta = np.deg2rad(theta_deg) 

    dx = (-1)**(invert+1) * (stars["ra"] - center_ra) * np.cos(np.deg2rad(center_dec))
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












def create_plot(stars, camera_width_pix, camera_height_pix, plot_lim_mag=12, plot_margin=300): 

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
    fig = plt.figure(figsize=(9, 6))
    ax = plt.gca() 
    plt.scatter(stars_new["x"], stars_new["y"], s=size, color="black")

    # Label axes and title 
    plt.xlabel("X (pixels)")
    plt.ylabel("Y (pixels)")
    plt.title(f"Stars with mag < {plot_lim_mag:.2f} (N={mask.sum()})") 

    # Limits: add margin beyond FOV of camera 
    plt.xlim((0-plot_margin, camera_width_pix+plot_margin))
    plt.ylim((0-plot_margin, camera_height_pix+plot_margin)) 

    # Tick labels: evenly spaced
    xticks = np.linspace(0, camera_width_pix, 5)
    ax.set_xticks(xticks.astype(int))
    yticks = np.linspace(0, camera_height_pix, 5)
    ax.set_yticks(yticks.astype(int))

    # Add red lines to show the camera FOV 
    plt.axvline(0, color="red") 
    plt.axvline(camera_width_pix, color="red")
    plt.axhline(0, color="red") 
    plt.axhline(camera_height_pix, color="red") 

    return fig, ax 





# # theta = np.deg2rad(theta_slider.value) 
# plot_lim_mag = mag_lim_slider.value
# Wrap current axes in reactive matplotlib widget
# ax = mo.ui.matplotlib(plt.gca())