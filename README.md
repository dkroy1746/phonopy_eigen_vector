This program is useful for visualising phonon eigen modes in crystal structures:
- In the pre-processing stage you generate the supercell file (s) using [Phonopy](https://phonopy.github.io/phonopy/)
- You then run [VASP](https://www.vasp.at/wiki/index.php/The_VASP_Manual) for calculating forces using either perturbation theory (IBRION 7 or 8 for undistorted single supercell) or using the finite difference method (IBRION -1 for multiple distorted structures)
- After a successful calculation you generate force constants [VASP & Phonopy](https://phonopy.github.io/phonopy/vasp.html)
- You calculate the phonon bandstructure
- Now you use the code "eigen_vect_phono.sh". The code requires "band.yaml" and POSCAR-unitcell in direct coordinates. 
  You input the frequency without rounding off. 
  Do note that only first q-point would be considered for a given frequency. You either give your phonon bandstructure path in such a way that the required q-point is towards the beginning or you delete other parts from band.yaml file.

- I have included an example band.yaml and POSCAR-unitcell.

Caution: Non-orthogonal lattice vectors would result in phonon modes whose directions may not be accurate.
