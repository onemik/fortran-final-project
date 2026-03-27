program  main
  use NumberKinds
  use Types
  use InputOutput
  use MoleculeBuilder
  use ComputeEnergy
  use Metropolis
  implicit none 

  type(Molecule) :: mol

  !read molecule from file
  call ReadFile(mol)

  !build molecular structure
  call BuildMolecule(mol)

  !set initial energy (before optimisation)
  call ComputeAllEnergies(mol)
  
end program 
