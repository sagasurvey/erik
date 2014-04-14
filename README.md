saga-code
=========

Code for the SAGA survey.

How this code is organized
==========================

The basic structure is based on ``.py`` files and directories, named after
either observatories or tasks.  For example, the `hosts.py` file defines data
structures and functions for the SAGA hosts, including downloading SDSS
catalogs of their environs.  The `targeting.py` file contains the code for
selecting targets from those catalogs.  The `magellan.py`, `mmthecto.py`, and
`wiyn.py` files all contain code to assist in generating configurations/masks
for the respective instruments.

These codes are each designed to place their outputs in appropriate
directories.  For example, `wiyn.py` always places configuration files for
HYDRA in the `wiyn_targets` directory, and similar for `magellan.py` and the
`imacs_targets` directory. `hosts.py` and `targeting.py` use the `catalogs`
directory, which will automatically be created if not present when an SDSS
catalog download is requested.

How to use the code
===================

For detailed instructions on the various pieces, see the ``*.py`` file
comments and docstrings. In general, though, they are intended to be used with
ipython notebooks, so that we can keep a record of how targets are selected,
objects are analyzed, etc.

The prerequesites are:

* Python (2.7+), numpy, scipy, and matplotlib - these can usually be installed
  with a package manager if you have one, but if you're not sure, google (or
  asking someone like Erik) will generally get you the necessary installation
  functions.
* [Astropy](http://www.astropy.org/) - the front page of the web site has
  [detailed install instructions.
* [IPython](http://ipython.org/) - you'll need a version that supports html
  [notebooks to follow the instructions below.  I think that means v0.13 or
  [later, but your best bet is just the latest (v2.0).

If you have already downloaded the NSA catalog, you might want to symlink it
to the saga directory (with the default name: ``nsa_v0_1_2.fits``).  The code
will automatically download it, but it's a pretty big file, so you probably
don't want more than one.

Once you've got everything, you can cd into the ``ipython_notebooks``
directory, and start ipython as ``ipython notebook --matplotlib=inline``. You
can then either run pre-existing notebooks, or create new ones if you're
preparing for a new observing run or doing some new analysis task.
