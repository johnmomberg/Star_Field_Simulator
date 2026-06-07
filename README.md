# Star Field Simulator 

Molab link currently unavailable 

The goal of this project is to simulate the stars visible in a given telescope observation. It is designed to help verify that you are pointing at the correct target.

First, choose a pointing direction (RA/Dec), which defines the center of the field, along with your telescope’s detector size and plate scale.

The tool uses VizieR to retrieve the RA, Dec, and magnitude of all stars within a circular field of view. It queries both Gaia and Hipparcos catalogs to ensure coverage across a wide range of magnitudes.

It then converts the sky coordinates (RA/Dec) into image coordinates (x/y) and generates a simulated view of what the detector would observe.

You can use sliders to dynamically adjust:

- the rotation angle,
- the limiting magnitude, 
- and invert the image

The interactive sliders make it easy to adjust orientation and brightness limits until the simulated field aligns with what you see in your telescope.

The idea is to start by displaying only the brightest stars in the field and adjust the rotation and orientation until the simulated star pattern matches your telescope image. You can then gradually increase the limiting magnitude to include fainter stars, helping confirm the match and revealing the location of fainter stars if your telescope is sensitive enough to detect them. 

This approach allows you to dynamically control the level of detail in the star field, rather than relying on VizieR, which can contain an overwhelming number of faint sources that your telescope can't detect anyway, which makes it difficult to use it for alignment. 
