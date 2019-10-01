.. _software-reqs:

Software Requirements
=====================

======
Python
======

Our scripts are tested on python 3.6. Most will likely work with lower versions
of python but some external modules may present issues that you'll need to
solve independently. We recommend downloading and installing `anaconda
<https://www.anaconda.com/>`_, which ships with most of the heavily used
modules in LLC_Membranes.

Some modules that you may need to install manually include:

* `mdtraj <http://mdtraj.org>`_
* `tqdm <https://tqdm.github.io/>`_
* `ruptures <http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html>`_
* `fbm <https://github.com/crflynn/fbm>`_

All are available on `PyPI <https://pypi.org>`_ and can therefore be installed using pip.

=======
GROMACS
=======

Most of the python packages in LLC_Membranes have been tested on trajectories
and coordinate files generated using GROMACS version 2018.3.  The scripts
should work with any version of GROMACS greater than or equal to 5.x (when the
all tools became modules of a binary names gmx).

See their documentation for `Installation instructions <http://http://manual.gromacs.org/documentation/2018.3/install-guide/index.html>`_.

==========
AmberTools
==========

We parameterize monomers and solutes using the generalized amber force field
(GAFF). One can parameterize organic molecules with GAFF using the `antechamber
<http://ambermd.org/antechamber/ac.html>`_ package which comes with `AmberTools
<http://ambermd.org/AmberTools.php>`_.

See their `download page <http://ambermd.org/GetAmber.php#ambertools>`_ to obtain 
the software for you machine.

