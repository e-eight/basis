/****************************************************************
  jjjpn_scheme_general_test.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>

#include "jjjpn_scheme_general.h"

////////////////////////////////////////////////////////////////
// test code
////////////////////////////////////////////////////////////////

void test_orbitals_Nmax(const std::string& filename)
{

  ////////////////////////////////////////////////////////////////
  // orbital tests (should come after suborbital tests are passed)
  ////////////////////////////////////////////////////////////////

  std::cout << "Orbitals -- Nmax scheme" << std::endl;

  // set up space
  std::cout << "Space" << std::endl;
  int Nmax = 4;
  basis::OrbitalSpacePN space(Nmax);
  std::cout << space.DebugStr();

  // check subspaces
  std::cout << "Subspaces" << std::endl;
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const basis::OrbitalSubspacePN& subspace = space.GetSubspace(subspace_index);

      std::cout << " index " << subspace_index
                << " species " << int(subspace.orbital_species())
                << std::endl;

      std::cout << subspace.DebugStr();
    }

  // check file output
  std::ofstream os(filename.c_str());
  os << space.OrbitalDefinitionStr();

}

void test_two_body_Nmax()
{

  ////////////////////////////////////////////////////////////////
  // two-body basis tests
  ////////////////////////////////////////////////////////////////

  std::cout << "Two-body basis -- Nmax scheme" << std::endl;

  // set up orbitals
  int orbital_Nmax = 4;
  basis::OrbitalSpacePN orbital_space(orbital_Nmax);

  // example subspace
  std::cout << "Example subspace" << std::endl;
  std::cout << "  basis::TwoBodySpeciesPN::kPN,2,0,basis::WeightMax(2,4)" << std::endl;
  basis::TwoBodySubspaceJJJPN subspace(
      orbital_space,
      basis::TwoBodySpeciesPN::kPN,2,0,
      basis::WeightMax(2,4)
    );
  std::cout << subspace.DebugStr();
  std::cout << "Orbital subspace sizes"
            << " " << subspace.orbital_subspace1().size()
            << " " << subspace.orbital_subspace2().size()
            << std::endl;

  // set up space

  std::cout << "Two-body space" << std::endl;
  std::cout << "      basis::WeightMax(2,2)" << std::endl;
  basis::TwoBodySpaceJJJPN space(
      orbital_space,
      basis::WeightMax(2,2)
    );
  std::cout << space.DebugStr();

  // then set up allowed sectors
  std::cout << "Two-body operator sectors" << std::endl;
  int J0 = 0;  // try: J0=0 for interaction, J0=2 for quadrupole operator
  int g0 = 0;
  basis::TwoBodySectorsJJJPN sectors(space,J0,g0);

  std::cout << " J0 " << J0 << " g0 " << g0 << std::endl;
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
    {
      int bra_subspace_index = sectors.GetSector(sector_index).bra_subspace_index();
      const basis::TwoBodySubspaceJJJPN& bra_subspace = sectors.GetSector(sector_index).bra_subspace();
      int ket_subspace_index = sectors.GetSector(sector_index).ket_subspace_index();
      const basis::TwoBodySubspaceJJJPN& ket_subspace = sectors.GetSector(sector_index).ket_subspace();

      std::cout 
        << " sector " 
        << std::setw(3) << sector_index 
        << "     "
        << " index "
        << std::setw(3) << bra_subspace_index
        << " sJg "
        << std::setw(3) << int(bra_subspace.two_body_species())
        << std::setw(3) << bra_subspace.J() 
        << std::setw(3) << bra_subspace.g()
        << " dim "
        << std::setw(3) << bra_subspace.size()
        << "     "
        << " index "
        << std::setw(3) << ket_subspace_index
        << " sJg "
        << std::setw(3) << int(ket_subspace.two_body_species())
        << std::setw(3) << ket_subspace.J() 
        << std::setw(3) << ket_subspace.g()
        << " dim "
        << std::setw(3) << ket_subspace.size()
        << std::endl;
    }
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string filename("test/jjjpn_scheme_general_test_orbitals_Nmax04.dat");
  test_orbitals_Nmax(filename);
  test_two_body_Nmax();

  // termination
  return 0;
}
