# 3D Profile Fitting for Neutron TOF Bragg Peaks

This repository contains files to integrate time-of-flight (TOF) neutron single crystal peak intensities using an Ikeda-Carpenter function and a bivariate gaussian.  It is built on the [Mantid framework](http://www.mantidproject.org/Main_Page) and standard Python libraries.  

## How can I use this program?
There is a getting started guide in the docs folder which explains what most of the parameters used for fitting and, with some tuning, should be usable on any TOF neutron crystallography which Mantid can open convert to an MDWorkspace.  Anyone interested in using this repository without Mantid should contact the author.

The scripts ```doICCFit.py``` and ```doBVGFit.py``` are meant to be scripts that can be quickly changed to allow different integration parameters (sample runs, output names, integration parameters).  Ideally, everything that needs to be changed will be at the top of the file, though this may not be the case for features still under development.

## Why use profile fitting?
Neutron crystallography is a powerful technique but remains hindered by high backgrounds and weak scattering cross sections. Bragg peaks on TOF diffractometers are complex, three dimensional intensity profiles in reciprocal space.  While these peaks have traditionally been integrated using peak-minus-background approaches, which rely on counting neutrons in a defined volume and subtracting a geometrically-scaled average background, this approach suffers from several critical disadvantages.  First, the TOF profile is a complex, asymmetric shape which depends on moderator characteristics.  Second, it assumes perfect peak prediction and detector calibration; for small protein peaks being only several pixels off can affect peak intensity upto 50%.  Finally, it fails to reliably integrate peaks that fall near detector edges which must either discarded or incorrectly integrated.  This can affect upto one third of all peaks at some beamlines.  By generating a 3D model of each peak and integrating that peak, these issues can 
be avoided.
