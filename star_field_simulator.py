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
    import matplotlib.pyplot as plt 

    from astropy.coordinates import SkyCoord
    import astropy.units as u




    return SkyCoord, mo, np, plt, src, u


@app.cell
def _(mo):
    # 1: Load fits file displayed in background of plot  
    file_selector_textbox = mo.ui.text(
        value="D:/Observations_(Fits_Files)/SBIG_Speed_Test/NoROI_raw/NoROI_1x1_30sec-0001.fit", 
        label="Full path to FITS file", 
        full_width=True, 
    )

    file_selector_textbox

    return (file_selector_textbox,)


@app.cell
def _(file_selector_textbox, src):
    file = src.load_fits(file_selector_textbox.value)
    return (file,)


@app.cell
def _(mo):
    # 2: Create and display input windows for loading stars 

    plate_scale_textbox = mo.ui.number(
        value=0.63, 
        label="Plate scale (arcsec/pixel)", 
    )

    camera_width_pix_textbox = mo.ui.number(
        value=3072, 
        label="Camera width (pixels)", 
    )

    camera_height_pix_textbox = mo.ui.number(
        value=2048, 
        label="Camera height (pixels)", 
    )

    target_ra_textbox = mo.ui.text(
        value="15:13:16.82",
        label="Target RA (HH:MM:SS)",
    )

    target_dec_textbox = mo.ui.text(
        value="+25:41:10.7",
        label="Target Dec (±DD:MM:SS)",
    )

    mo.vstack([plate_scale_textbox, camera_width_pix_textbox, camera_height_pix_textbox, target_ra_textbox, target_dec_textbox]) 
    return (
        camera_height_pix_textbox,
        camera_width_pix_textbox,
        plate_scale_textbox,
        target_dec_textbox,
        target_ra_textbox,
    )


@app.cell
def _(
    SkyCoord,
    camera_height_pix_textbox,
    camera_width_pix_textbox,
    plate_scale_textbox,
    target_dec_textbox,
    target_ra_textbox,
    u,
):
    # 3: Access values from input windows and convert to usable forms 

    plate_scale = plate_scale_textbox.value 
    camera_width_pix = camera_width_pix_textbox.value 
    camera_height_pix = camera_height_pix_textbox.value 

    target_skycoord = SkyCoord(
        ra=target_ra_textbox.value,
        dec=target_dec_textbox.value,
        unit=(u.hourangle, u.deg)
    )

    target_ra = target_skycoord.ra.deg 
    target_dec = target_skycoord.dec.deg 
    return (
        camera_height_pix,
        camera_width_pix,
        plate_scale,
        target_dec,
        target_ra,
    )


@app.cell
def _(
    camera_height_pix,
    camera_width_pix,
    np,
    plate_scale,
    src,
    target_dec,
    target_ra,
):
    # 4: Load all stars within 2xFOV of target position 



    # Calculate radius to load stars based on FOV of camera 
    term1 = plate_scale/3600*camera_width_pix 
    term2 = plate_scale/3600*camera_height_pix 
    fov_deg = np.sqrt((term1/2)**2 + (term2/2)**2) 
    radius_deg = 2*fov_deg 

    # Retrieve stars from Hipparcos (bright)
    stars_hipparcos = src.retrieve_stars(
        center_ra=target_ra, 
        center_dec=target_dec, 
        radius_deg=radius_deg, 
        catalog="I/239/hip_main", 
    )

    # Retrive stars from Gaia (dim) 
    stars_gaia = src.retrieve_stars(
        center_ra=target_ra, 
        center_dec=target_dec, 
        radius_deg=radius_deg, 
        catalog="I/355/gaiadr3", 
        mag_limit=18 
    )

    # Combine all stars into single dict 
    stars_no_xy = {
        "ra":  np.concatenate([stars_gaia["ra"], stars_hipparcos["ra"]]),
        "dec": np.concatenate([stars_gaia["dec"], stars_hipparcos["dec"]]),
        "mag": np.concatenate([stars_gaia["mag"], stars_hipparcos["mag"]]),
    }
    return fov_deg, stars_no_xy


@app.cell
def _(fov_deg, mo, target_dec, target_ra):
    # Coarse sliders 


    # Magnitude limit slider 
    mag_lim_slider = mo.ui.slider(
        start=0,
        stop=18,
        value=12, 
        step=1, 
        label="Limiting magnitude",
        full_width=True, 
        show_value=True, 
    ) 

    plot_margin_slider = mo.ui.slider(
        start=0, 
        stop=5000, 
        value=400, 
        step=100, 
        label="Plot margin", 
        full_width=True, 
        show_value=True, 
    ) 

    # Theta slider 
    theta_slider_coarse = mo.ui.slider(
        start=-180, 
        stop=360, 
        value=0, 
        step=5, 
        label="Theta", 
        full_width=True, 
        show_value=True, 
    )

    # Invert image 
    invert_checkbox = mo.ui.checkbox(
        value=False, 
        label="Invert",
    )

    show_image = mo.ui.checkbox(
        value=False, 
        label="Show image",
    )


    plot_center_ra_slider_coarse = mo.ui.slider(
        start=round(target_ra - fov_deg, 4),
        stop=round(target_ra + fov_deg, 4),
        value=round(target_ra, 4),
        step=round(fov_deg / 10, 4),
        label="Plot center RA offset (deg)",
        full_width=True,
        show_value=True,
    )

    plot_center_dec_slider_coarse = mo.ui.slider(
        start=round(target_dec - fov_deg, 4),
        stop=round(target_dec + fov_deg, 4),
        value=round(target_dec, 4),
        step=round(fov_deg / 10, 4),
        label="Plot center dec offset (deg)",
        full_width=True,
        show_value=True,
    )


    return (
        invert_checkbox,
        mag_lim_slider,
        plot_center_dec_slider_coarse,
        plot_center_ra_slider_coarse,
        plot_margin_slider,
        show_image,
        theta_slider_coarse,
    )


@app.cell
def _(
    mo,
    plot_center_dec_slider_coarse,
    plot_center_ra_slider_coarse,
    theta_slider_coarse,
):
    # Fine sliders 

    theta_slider_fine = mo.ui.slider(
        start=theta_slider_coarse.value-theta_slider_coarse.step*2, 
        stop=theta_slider_coarse.value+theta_slider_coarse.step*2, 
        value=theta_slider_coarse.value, 
        step=0.1, 
        label="Theta (fine)", 
        full_width=True, 
        show_value=True, 
    )

    plot_center_ra_slider_fine = mo.ui.slider(
        start=round(plot_center_ra_slider_coarse.value-plot_center_ra_slider_coarse.step*2, 4),
        stop=round(plot_center_ra_slider_coarse.value+plot_center_ra_slider_coarse.step*2, 4),
        value=round(plot_center_ra_slider_coarse.value, 4),
        step=round(plot_center_ra_slider_coarse.step/100, 5),
        label="Plot center RA offset (deg) (fine)",
        full_width=True,
        show_value=True,
    )

    plot_center_dec_slider_fine = mo.ui.slider(
        start=round(plot_center_dec_slider_coarse.value-plot_center_dec_slider_coarse.step*2, 4),
        stop=round(plot_center_dec_slider_coarse.value+plot_center_dec_slider_coarse.step*2, 4),
        value=round(plot_center_dec_slider_coarse.value, 4),
        step=round(plot_center_dec_slider_coarse.step/100, 5),
        label="Plot center dec offset (deg) (fine)",
        full_width=True,
        show_value=True,
    )

    return (
        plot_center_dec_slider_fine,
        plot_center_ra_slider_fine,
        theta_slider_fine,
    )


@app.cell
def _(
    invert_checkbox,
    mag_lim_slider,
    mo,
    plot_center_dec_slider_coarse,
    plot_center_dec_slider_fine,
    plot_center_ra_slider_coarse,
    plot_center_ra_slider_fine,
    plot_margin_slider,
    show_image,
    theta_slider_coarse,
    theta_slider_fine,
):
    mo.vstack([
        mag_lim_slider, 
        mo.hstack([plot_margin_slider, show_image]), 
        mo.hstack([theta_slider_coarse, theta_slider_fine, invert_checkbox]), 
        mo.hstack([plot_center_ra_slider_coarse, plot_center_ra_slider_fine]), 
        mo.hstack([plot_center_dec_slider_coarse, plot_center_dec_slider_fine])
    ]) 
    return


@app.cell
def _(
    invert_checkbox,
    mag_lim_slider,
    plot_center_dec_slider_fine,
    plot_center_ra_slider_fine,
    plot_margin_slider,
    theta_slider_fine,
):
    # 6: Access values from input windows 
    plot_lim_mag = mag_lim_slider.value 
    plot_margin = plot_margin_slider.value 
    theta_deg = theta_slider_fine.value 
    invert = invert_checkbox.value 

    plot_center_ra = plot_center_ra_slider_fine.value 
    plot_center_dec = plot_center_dec_slider_fine.value 
    return (
        invert,
        plot_center_dec,
        plot_center_ra,
        plot_lim_mag,
        plot_margin,
        theta_deg,
    )


@app.cell
def _(
    camera_height_pix,
    camera_width_pix,
    invert,
    plate_scale,
    plot_center_dec,
    plot_center_ra,
    src,
    stars_no_xy,
    target_dec,
    target_ra,
    theta_deg,
):
    # 7: Calculate xy coordinates on detector for all stars and target position 

    stars = src.convert_skycoords_to_pixelcoords(
        stars_no_xy, 
        plot_center_ra, 
        plot_center_dec, 
        theta_deg, 
        camera_width_pix, 
        camera_height_pix, 
        plate_scale, 
        invert=invert, 
    ) 

    target = src.convert_skycoords_to_pixelcoords(
        {"ra": target_ra, "dec": target_dec}, 
        plot_center_ra, 
        plot_center_dec, 
        theta_deg, 
        camera_width_pix, 
        camera_height_pix, 
        plate_scale, 
        invert=invert, 
    ) 
    return stars, target


@app.cell
def _(
    camera_height_pix,
    camera_width_pix,
    file,
    plot_lim_mag,
    plot_margin,
    plt,
    show_image,
    src,
    stars,
    target,
):
    # 8: Create plot 

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot() 

    # Plot fits file in background 
    if show_image.value==True: 
        src.plot_fits(ax, file)

    # Plot stars on top of image 
    src.plot_stars(
        ax, 
        stars, 
        camera_width_pix, 
        camera_height_pix, 
        plot_lim_mag=plot_lim_mag, 
        plot_margin=plot_margin 
    )

    # Add crosshair at target position 
    plt.scatter(
        target["x"],
        target["y"],
        facecolors='none',
        edgecolors="tab:green", 
        s=50, 
        label=f"Target: {int(target['x'])}, {int(target['y'])}"
    )
    plt.legend()
    return


@app.cell
def _(mo, target):
    mo.md(f"""
    Target position: {target['x']:.1f}, {target['y']:.1f}
    """)
    return


if __name__ == "__main__":
    app.run()
