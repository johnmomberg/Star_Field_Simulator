# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "astroquery==0.4.11",
#     "numpy==2.4.6",
# ]
# ///

import marimo

__generated_with = "0.23.9"
app = marimo.App(
    width="medium",
    css_file="/usr/local/_marimo/custom.css",
    auto_download=["html"],
)


@app.cell
def _():
    import src 
    import marimo as mo 
    import numpy as np 


    return mo, np, src


@app.cell
def _(
    camera_height_pix_textbox,
    camera_width_pix_textbox,
    center_dec_textbox,
    center_ra_textbox,
    plate_scale_textbox,
):
    # Input 

    plate_scale = plate_scale_textbox.value  
    camera_width_pix = camera_width_pix_textbox.value 
    camera_height_pix = camera_height_pix_textbox.value 

    center_ra = center_ra_textbox.value 
    center_dec = center_dec_textbox.value 

    return (
        camera_height_pix,
        camera_width_pix,
        center_dec,
        center_ra,
        plate_scale,
    )


@app.cell
def _(
    camera_height_pix,
    camera_width_pix,
    center_dec,
    center_ra,
    np,
    plate_scale,
    src,
):
    # Calculate FOV of camera 
    term1 = plate_scale/3600*camera_width_pix 
    term2 = plate_scale/3600*camera_height_pix 
    radius_deg = np.sqrt((term1/2)**2 + (term2/2)**2) 

    # Retrieve stars from Hipparcos (bright)
    stars_hipparcos = src.retrieve_stars(
        center_ra=center_ra, 
        center_dec=center_dec, 
        radius_deg=radius_deg, 
        catalog="I/239/hip_main", 
    )

    # Retrive stars from Gaia (dim) 
    stars_gaia = src.retrieve_stars(
        center_ra=center_ra, 
        center_dec=center_dec, 
        radius_deg=radius_deg, 
        catalog="I/355/gaiadr3", 
    )

    # Combine all stars into single dict 
    stars_combined = {
        "ra":  np.concatenate([stars_gaia["ra"], stars_hipparcos["ra"]]),
        "dec": np.concatenate([stars_gaia["dec"], stars_hipparcos["dec"]]),
        "mag": np.concatenate([stars_gaia["mag"], stars_hipparcos["mag"]]),
    }
    return (stars_combined,)


@app.cell
def _(
    camera_height_pix,
    camera_width_pix,
    center_dec,
    center_ra,
    invert_checkbox,
    plate_scale,
    src,
    stars_combined,
    theta_slider,
):
    # Calculate x, y coordinates of all stars 
    stars = src.convert_skycoords_to_pixelcoords(
        stars_combined, 
        center_ra, 
        center_dec, 
        theta_slider.value, 
        camera_width_pix, 
        camera_height_pix, 
        plate_scale, 
        invert=invert_checkbox.value, 
    ) 

    return (stars,)


@app.cell
def _(camera_height_pix, camera_width_pix, mag_lim_slider, src, stars):
    # Create plot 
    fig, ax = src.create_plot(
        stars, 
        camera_width_pix, 
        camera_height_pix, 
        plot_lim_mag=mag_lim_slider.value, 
    )

    return (ax,)


@app.cell
def _(mo):
    # Input value controls 


    center_ra_textbox = mo.ui.number(
        value=279.23,
        label="Center RA (deg)",
    )

    center_dec_textbox = mo.ui.number(
        value=38.78,
        label="Center Dec (deg)",
    )

    camera_width_pix_textbox = mo.ui.number(
        value=3072, 
        label="Camera width (pixels)", 
    )

    camera_height_pix_textbox = mo.ui.number(
        value=2048, 
        label="Camera height (pixels)", 
    )

    plate_scale_textbox = mo.ui.number(
        value=0.63, 
        label="Plate scale (arcsec/pixel)", 
    )
    return (
        camera_height_pix_textbox,
        camera_width_pix_textbox,
        center_dec_textbox,
        center_ra_textbox,
        plate_scale_textbox,
    )


@app.cell
def _(mo):
    # Plot controls 


    # Magnitude limit slider 
    mag_lim_slider = mo.ui.slider(
        start=4,
        stop=18,
        value=12, 
        step=0.1, 
        label="Mag",
        full_width=True, 
        show_value=True, 
        orientation="vertical", 
    )

    # Theta slider 
    theta_slider = mo.ui.slider(
        start=0, 
        stop=360, 
        value=0, 
        step=0.1, 
        label="Theta", 
        full_width=True, 
        show_value=True, 
        orientation="vertical", 
    )

    # Invert image 
    invert_checkbox = mo.ui.checkbox(
        value=False, 
        label="Invert",
    )
    return invert_checkbox, mag_lim_slider, theta_slider


@app.cell
def _(
    camera_height_pix_textbox,
    camera_width_pix_textbox,
    center_dec_textbox,
    center_ra_textbox,
    mo,
    plate_scale_textbox,
):
    # Input value controls 
    mo.vstack([center_ra_textbox, center_dec_textbox, camera_width_pix_textbox, camera_height_pix_textbox, plate_scale_textbox]) 
    return


@app.cell
def _(ax, invert_checkbox, mag_lim_slider, mo, theta_slider):
    # Plot 
    ax_mo = mo.ui.matplotlib(ax) 
    controls = mo.vstack([mag_lim_slider, theta_slider, invert_checkbox], align="start", justify="center") 
    mo.hstack([ax_mo, controls], widths=[7, 1]) 
    return


if __name__ == "__main__":
    app.run()
