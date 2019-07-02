Cross-linking
=============

Our cross-linking algorithm is designed to simulate the result of a 
free radical polymerization (FRP) between vinyl group containing monomers.
There are three main components to a free radical polymerization:

1. Initiation

        In an FRP, an initiator decomposes into free radicals which readily
        react with one electron from the pi bond of a C=C group, creating a new bond
        with one of the C atoms. The extra pi-bond electron moves to the other C atom
        which becomes the new active radical site.

        In an an effort to simplify our proces, we represent the initiator as hydrogen
        so that we only have to add single hydrogen atoms at initiation sites.

        Atoms are added to the system during cross-linking. In order to make this work
        efficiently with molecular simulations, we added dummy hydrogen atoms to all
        sites where an initiator could potentially be added. The dummy atoms are only
        made real when necessary and are eliminated after cross-linking has terminated.

2. Progagation

        During propagation, monomers sequentially add to active centers.
        Monomers with radical components react with other monomers in the same way as
        the intiator. 

3. Termination

        For the reaction to be terminated, all of the radicals must react.  The
        two most common termination mechanisms are combination and disproportionation.
        Combination involves two monomers with active centers coupling together.
        Disproportionation involves abstraction of a hydrogen atom by one radical
        species from another species.

Our cross-linking algorithm happens iteratively. Each iteration, a group of
potentially bonding atoms are selected based on their proximity to each other.
We update the topology file with the new atom types, bonds, angles, etc, then
perform an energy minimization to zip together the new bonds and finally a
short NVT simulation. This process is repeated until the convergence criteria,
the number of cross-links, is met.

This algorithm is meant to be used with any vinyl containing system, however
the user will have to define any undefined reactions. 

See :ref:`xlink-schemes`

.. argparse:: 
   :filename: ../LLC_Membranes/setup/xlink.py
   :func: initialize
   :prog: xlink.py

=======
Classes
=======

.. autoclass:: xlink.System
   :show-inheritance:
   :members: __init__

.. autoclass:: xlink.Topology
   :members: __init__

