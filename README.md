# ESE156 Remote Sensing Class

**Instructor: Christian Frankenberg**; cfranken@caltech.edu (203 Linde+Robinson), Office hours: Just drop by (usually around after 9am).

**TA: Liyin He:** lhe@caltech.edu, office hours TBD

**Course Material:** Most Python Notebooks on github, additional material on ftp://fluo.gps.caltech.edu/XYZT_ESE156/

To cite Numerical Recipes:
“To be genuinely useful, a fitting procedure should provide (i) parameters,
(ii) error estimates on the parameters, and (iii) a statistical
measure of goodness-of-fit. When the third item suggests that the
model is an unlikely match to the data, then items (i) and (ii) are
probably worthless. Unfortunately, many practitioners of parameter
estimation never proceed beyond item (i). They deem a fit acceptable if a graph of data and model ‘looks good’. This approach is known as
chi-by-eye. Luckily, its practitioners get what they deserve.”

We will try to cover all aspects mentioned above in our Remote Sensing Class here, which is essentially a form of parameter estimation and bears similarity to fitting a line through a set of measured (x,y) points.

## Course Overview

ESE156 will cover the basics of remote sensing with a focus on atmospheric and surface remote sensing using spectroscopy. The class will guide you through typical steps of relevance to all remote sensing problem, i.e. how to set up and formulate your "Forward model", applying constrained and un-constrained, linear and non-linear inverse modeling and perform error characterization of measurements and inverted properties. The focus on this class is more applied, using a scripting language like python to solve numerical problems you might have encountered in more theoretical classes. Even though Python is the preferred language for this class, you can use Julia or Matlab or "real" programming languages such as C/C++, Fortran, etc). Typically, each 1.5hrs class consists of a mix of slides and jupyter notebooks, with interactive class involvement in the last 30min. The class syllabus is still not written in stone but follows these general steps:

**Week 1:** Basic introduction and simple inverse technique by linearizing Lambert-Beer’s law

**Week 2:** Computing layer optical properties using atmospheric model output as well as cross section databases (line-by-line modeling)

**Week 3:** Simplified synthetic retrievals of total columns of strong absorbers (using TCCON as analogy), generating a non-linear forward model and its derivatives

**Week 4:** Based on simple forward model developed in Week 3: Use of prior constraints, information content analysis, error analysis, impact of spectroscopic uncertainties…

**Week 5:** Impact of scattering (Rayleigh, aerosols, clouds) on observations from space

**Week 6:** Nadir sounders in the thermal infrared: Radiative transfer and retrieval principles.

**Week 7:** Solar Induced Chlorophyll Fluorescence and Vegetation Indices.

**Weeks 8-10:** Hands-on examples or paper.

## Homework/Projects

We will use weekly hands-on problem sets, which require python coding in Jupyter. Problem sets can be submitted as Jupyter notebook html export.  

## Workflow

If you're using GitHub Desktop, these general instructions will help:

* <https://guides.github.com/activities/forking/>
* <https://help.github.com/desktop/guides/contributing/>

## Getting started with python: 

The ESE156 Remote Sensing Class will make extensive use of the Jupyter notebook environment as well as python as a coding tool. If you aren't familiar with Python yet, there are good resources online, e.g.

Install for beginners: https://www.python.org/about/gettingstarted/

Data Science: https://jakevdp.github.io/PythonDataScienceHandbook/

Jupyter: https://www.dataquest.io/blog/jupyter-notebook-tutorial/

Extensions for Jupyter notebook: https://jupyter-contrib-nbextensions.readthedocs.io/en/latest/install.html

My goal is to have everyone set up with python and jupyter notebook in week 1. If needed, we will go through some basic python steps in the beginning as well (and/or communicate with the TA). 

## Grading

* Class Participation – 20%
* Homework – 50%
* Project - 30%

