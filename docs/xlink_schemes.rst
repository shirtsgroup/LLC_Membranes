.. _xlink-schemes:

Designing Cross-Linking Reactions
=================================

Programming cross-linking reactions can be very tedious because it requires
consistent numbering, accurate chemistry and adding and subtracting bonds,
dihedrals, angles, virtual sites etc. The classes in xlink_schemes.py attempt
to simplify the process so that minimal work is required to cross-link a new
system. As of now, only the :ref:`diene scheme<diene-scheme>` used for the QI
phase is implemented. In order to implement new reaction schemes, one should
follow this reaction as an example.

In addition to designing reactions, different monomers, even those with the
same cross-linkable groups, may be number differently. This is hard to avoid
without manually renumbering everything. Instead of manually renumbering, one
can define a class for a monomer which defines all of the necessary monomer
attributes. Follow the recommendations in the docstring for the QI monomer
:ref:`Dibrpyr14<monomer-example>` as an example.

=======
Classes
=======

.. autoclass:: xlink_schemes.XlinkReaction
   :members: __init__

.. _diene-scheme:

.. autoclass:: xlink_schemes.DieneScheme
   :members: __init__, determine_reaction_type, head2tail, radical_c2, terminate

.. _monomer-example:

.. autoclass:: xlink_schemes.Dibrpyr14
   :members: __init__

