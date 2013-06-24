IsoRender
=========
An isosurface visualization program which dynamically computes isosurfaces, allowing the user to observe changes in the isosurface over incremental changes in the isovalue. 
Version 1.1

Copyright 2010 Marc Khoury
Please send any bug reports to marc.khry@gmail.com

IsoRender reads in a scalar grid from a nhdr file and computes an isosurface for a given isovalue using the IJK library. The resulting mesh is then displayed and can be manipulated or written to a file in either the OFF or OBJ mesh formats. 

![ScreenShot](https://raw.github.com/marckhoury/IsoRender/master/screenshot.png)

Build
=====
IsoRender uses the IJK (Isosurface Jeneration Kode) library to construct isosurfaces. Specifically ijkmcube.v0.3.0.tar, which can be downloaded at http://www.cse.ohio-state.edu/research/graphics/isotable/. Extract all the files from this archive into your IsoRender directory. Set the environment variable "IJK_ISOTABLE_DIR" to the location of the isosurface lookup tables, located in the isotable folder. The command looks similar to the following.

export IJK_ISOTABLE_DIR="~/Desktop/IsoRender/isotable"

IsoRender uses OpenGL and GLUI to display the images and create the GUI, as well as zlib and expact for decompressing input and parsing XML files. Be sure to check that all these libraries are installed before running the make command.

Dependencies
============
OpenGL, GLUT, GLUI, ijkmcube, zlib, expat, libITKNrrdIO (library file included).

If you find IsoRender useful in any way I'd love to hear about it. 

Marc Khoury
July 4, 2010


