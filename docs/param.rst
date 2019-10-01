.. _parameterize:

Parameterize Residues
======================

With the exception of water, we parameterize all monomers and solutes with
the general Amber force field (GAFF). To add GAFF parameters, we use the
`antechamber <http://ambermd.org/antechamber/ac.html>`_ package which comes with
`AmberTools <http://ambermd.org/AmberTools.php>`_. See their `download page
<http://ambermd.org/GetAmber.php#ambertools>`_ to obtain the software for your
machine. `antechamber`'s primary role is to identify GAFF atom types and
assign the corresponding bond, pair, angle, dihedral and improper dihedral
parameters. It also uses the AM1-BCC method in order to assign charges to the
atoms.

The AM1-BCC method is flawed in that it can easily assign charges asymetrically
to an otherwise symmetric molecule. We avoid this issue by re-assigning charges
with am1bccsym method of `molcharge
<https://docs.eyesopen.com/applications/quacpac/molcharge/molcharge.html>`_, a
python script shipped with `QUACPAC <https://www.eyesopen.com/quacpac>`_ from
OpenEye Scientific. 

Our general workflow for parameterizing a molecule with GAFF is:

#. Generate an initial structure. You will have the easiest time with a reasonable approximation of the 3D structure. I recommend using a tool like `MarvinSketch <https://chemaxon.com/products/marvin>`_ which is capable performing a crude energy minimization in order to arrange a molecule in 3D space.

#. Run antechamber to assign parameters and initial charges.

#. Energy minimize the structure using GAFF params.

#. Re-assign charges with molcharge.

#. Energy minimize.

#. (optional) Run an MD simulation in order to get new conformations. Sometimes this is necessary in order to get a monomer into a configuration that can easily be built into an H\ :sub:`II` unit cell.

The above steps are written in parameterize.py

.. code-block:: python

    from LLC_Membranes.setup import parameterize

======================
Command Line Interface
======================

.. argparse:: 
   :filename: ../LLC_Membranes/setup/parameterize.py
   :func: initialize
   :prog: parameterize.py

=======
Classes
=======

.. autoclass:: parameterize.Parameterize
   :members: __init__, assign_gaff_parameters, assign_molcharge_charges, convert_parameter_format, energy_minimize, insert_mol2_charges, make_box
