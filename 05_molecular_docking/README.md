# Molecular Docking & MD Simulation Workflow

This folder documents the computational workflow used to generate Figures 7 and 8 of the manuscript.

## Inputs

- **PFOA 3D structure**: `Conformer3D_COMPOUND_CID_9554.sdf` (downloaded from PubChem, CID 9554)
- **Protein structures** (RCSB PDB):
  - SHBG: PDB 7REK (or the structure used in the study)
  - FABP4, BCL6, CASP1, IL10, NR4A2: PDB structures as listed in the manuscript Methods

## Steps

1. **Protein preparation**: remove water molecules and ligands, add hydrogens (AutoDock Tools).
2. **Ligand preparation**: convert PFOA SDF to PDBQT (Open Babel / AutoDock Tools).
3. **Docking**: AutoDock Vina with a grid box covering the reported binding pocket; exhaustiveness = 20 (default). The top-ranked pose by binding affinity was selected for each target.
4. **Complex generation**: protein–ligand complex PDB files (e.g., `FABP4_Conformer3D_COMPOUND_CID_9554_out_1.-8.2.complex.pdb`) were built from the best docking pose.
5. **Molecular dynamics (GROMACS)**: 100 ns MD simulations of the FABP4–PFOA and SHBG–PFOA complexes (force field: as stated in the manuscript Methods); trajectory analyses (RMSD, RMSF, radius of gyration, H-bonds, SASA) were performed with GROMACS analysis tools (`.xvg` files in the authors' local archive).

## Notes

- Docking scores (binding affinities) for the six core targets are reported in the manuscript Results and Supporting Information.
- Docking/MD are hypothesis-generating; direct binding was experimentally confirmed by surface plasmon resonance (SPR) for FABP4 (Kd = 11 µM) and SHBG (Kd = 7.82 µM).
- Detailed parameter tables are available from the corresponding author; the primary outputs (complex structures, trajectories) are large and therefore not included in this repository.
