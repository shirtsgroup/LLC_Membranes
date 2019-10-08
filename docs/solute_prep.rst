.. _solute-prep:

Preparation of Solutes
======================

Similar to :ref:`monomer_prep`, we annotate solute topologies
to aid specific calculations. At present these annotations are
only needed to label solute atoms involved in hydrogen bonding.

.. _annotation-table-solutes:

Summary of Annotations
-----------------------

+-------+--------------------------------------------+-----------------------------+
|Symbol | Description                                |  More info                  |
+=======+============================================+=============================+
|H      |Specifies a hydrogen that can participate in|  :ref:`hbond-annotations`   |        
|       |a hydrogen bonding interaction              |                             |
+-------+--------------------------------------------+-----------------------------+
|D      |Specifies a hydrogen bond donor atom        |  :ref:`hbond-annotations`   |
+-------+--------------------------------------------+-----------------------------+
|A      |Specifies a hydrogen bond acceptor atom     |  :ref:`hbond-annotations`   |
+-------+--------------------------------------------+-----------------------------+

.. _hbond-annotations:

Hydrogen Bonding
----------------

In order to identify hydrogen bonds in a trajectory, one must specify which
hydrogen atoms are able to participate in such an interaction as well as which
atoms are hydrogen bond donors and which are hydrogen bond acceptors.
