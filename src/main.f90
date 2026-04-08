program  main
  use NumberKinds
  use Types
  use InputOutput
  use MoleculeBuilder
  use ComputeEnergy
  use Metropolis
  implicit none 

  type(Molecule) :: mol
  type(Energies) :: Einitial, Efinal

  integer :: NumberOfSteps
  real(KREAL) :: stepSize
  real(KREAL) :: temperature
  real(KREAL) :: bestEnergy
  real(KREAL) :: acceptanceRate

  !read molecule from file
  call ReadFile(mol)

  !build molecular structure
  call BuildMolecule(mol)

  !set initial energy (before optimisation)
  Einitial = ComputeAllEnergies(mol)
  
  !print starting information of the atom setup and energies
  call printStartingInfo(mol,Einitial)

  !ask user input for metropolis parameters
  call readMetropolisParams(NumberOfSteps, stepSize, temperature)

  !run metropolis
  call RunMetropolis(mol, NumberOfSteps, stepSize, temperature, bestEnergy, acceptanceRate)

  !calculate final energies
  Efinal = ComputeAllEnergies(mol)

  !print final information and results
  call printFinalInfo(Efinal, bestEnergy, acceptanceRate)
  call printResults(mol)

end program 
