Inside FM_potentials you can find potential tables for all FM potentials (non-bonded, pairwise) involved in the hybrid MARTINI/FM force field developed for ZIF-8/PVDF. These are labelled as CGX-CGY.dat, where X is the bead type for ZIF-8 and Y is the bead type for PVDF. In the case of ZIF-8, X follows the definition given in figure 1(a). In the case of PVDF, Y follows the definition given in figure 1(b) summed by 2.

In Nanoparticle_configurations you can find LAMMPS data files (atom_style format = FULL) containing a configuration for the ZIF-8 nanoparticle and its connectivity. Please note: (1) there is no polylmer around it and (2) the 212 angles are subcategorized in 2 different types, following the approach underlying the force field used, so if you want to use it for another force field, you need to pay attention to that.
