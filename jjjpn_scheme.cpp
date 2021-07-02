/****************************************************************
  jjjpn_scheme.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <cstddef>
#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "jjjpn_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJPN::TwoBodySubspaceJJJPN(
      const OrbitalSpacePN& orbital_space,
      TwoBodySpeciesPN two_body_species, int J, int g,
      const WeightMax& weight_max
    )
  {

    // std::cout << "Subspace construction"
    //           << " " << int(two_body_species)
    //           << " " << J
    //           << " " << g
    //           << std::endl;

    // set values
    labels_ = SubspaceLabelsType(two_body_species,J,g);
    weight_max_ = weight_max;

    // identify orbital subspaces
    OrbitalSpeciesPN orbital_species1 =
      (two_body_species==TwoBodySpeciesPN::kPP) || (two_body_species==TwoBodySpeciesPN::kPN)
      ? OrbitalSpeciesPN::kP : OrbitalSpeciesPN::kN;
    OrbitalSpeciesPN orbital_species2 =
      (two_body_species==TwoBodySpeciesPN::kPP)
      ? OrbitalSpeciesPN::kP : OrbitalSpeciesPN::kN;
    orbital_subspace1_ptr_ = &(orbital_space.GetSubspace(int(orbital_species1)));
    orbital_subspace2_ptr_ = &(orbital_space.GetSubspace(int(orbital_species2)));

    // set up indexing
    // iterate over orbitals for particle 1
    for (std::size_t index1 = 0; index1 < orbital_subspace1().size(); ++index1)
      // iterate over orbitals for particle 2
      for (std::size_t index2 = 0; index2 < orbital_subspace2().size(); ++index2)
        {

          // retrieve orbitals
          OrbitalStatePN orbital1(orbital_subspace1(),index1);
          OrbitalStatePN orbital2(orbital_subspace2(),index2);

          // for like-particle subspaces
          if ((two_body_species==TwoBodySpeciesPN::kPP) || (two_body_species==TwoBodySpeciesPN::kNN))
            {
              // impose canonical ordering on single-particle states
              if (!(index1<=index2))
                continue;

              // impose antisymmetry
              if ((index1==index2)&&(!(J%2==0)))
                continue;
            }

          // impose triangularity
          if (!(am::AllowedTriangle(orbital1.j(),orbital2.j(),J)))
            continue;

          // impose parity
          if (!((orbital1.g()+orbital2.g()+g)%2==0))
            continue;

          // impose weight cutoffs
          double w1 = orbital1.weight();
          double w2 = orbital2.weight();
          double w = w1 + w2;
          double w1_max = weight_max.one_body[int(orbital_species1)];
          double w2_max = weight_max.one_body[int(orbital_species2)];
          double w_max = weight_max.two_body[int(two_body_species)];
          bool weight_allowed = ((w1<=w1_max)&&(w2<=w2_max)&&(w<=w_max));
          // std::cout << " " << w1
          //           << " " << w2
          //           << " " << w
          //           << " " << w1_max
          //           << " " << w2_max
          //           << " " << w_max
          //           << " " << weight_allowed
          //           << std::endl;
          if (!weight_allowed)
            continue;

          // keep surviving states
          PushStateLabels(StateLabelsType(index1,index2));
        }
  }

  std::string TwoBodySubspaceJJJPN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(two_body_species())]
       << " " << std::setw(width) << J()
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }

  std::string TwoBodySubspaceJJJPN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {

        // retrieve two-body state
        TwoBodyStateJJJPN state(*this,state_index);

        // retrieve orbitals
        OrbitalStatePN orbital1(state.orbital_subspace1(),state.index1());
        OrbitalStatePN orbital2(state.orbital_subspace2(),state.index2());

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "index1 index2"
          << " " << std::setw(width) << state.index1()
          << " " << std::setw(width) << state.index2()
          << " " << "nlj1 nlj2"
          << " " << std::setw(width) << orbital1.n()
          << " " << std::setw(width) << orbital1.l()
          << " " << std::setw(width+2) << orbital1.j().Str()
          << " " << std::setw(width) << orbital2.n()
          << " " << std::setw(width) << orbital2.l()
          << " " << std::setw(width+2) << orbital2.j().Str()
          << std::endl;
      }

    return os.str();

  }

  std::string TwoBodyStateJJJPN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << kTwoBodySpeciesPNCodeTz[int(two_body_species())]
       << " " << std::setw(width) << J()
       << " " << std::setw(width) << g()
       << " " << ";"
       << " " << std::setw(width) << index1()
       << " " << std::setw(width) << index2()
       << " " << "("
       << " " << std::setw(width) << GetOrbital1().N()
       << " " << std::setw(width) << GetOrbital1().l()
       << " " << std::setw(width) << GetOrbital1().j()
       << " " << std::setw(width) << GetOrbital2().N()
       << " " << std::setw(width) << GetOrbital2().l()
       << " " << std::setw(width) << GetOrbital2().j()
       << " " << ")"
       << " " << "]";

    return os.str();
  }


  TwoBodySpaceJJJPN::TwoBodySpaceJJJPN(
      const OrbitalSpacePN& orbital_space,
      const WeightMax& weight_max,
      basis::TwoBodySpaceJJJPNOrdering ordering
    )
    : weight_max_(weight_max), space_ordering_(ordering) // orbital_space_(orbital_space)
  {

    // find putative Jmax from maximal j among orbitals
    HalfInt jmax=HalfInt(1,2);
    for (std::size_t subspace_index=0; subspace_index<orbital_space.size(); ++subspace_index)
      {
        // retrieve subspace
        const OrbitalSubspacePN& subspace = orbital_space.GetSubspace(subspace_index);

        // iterate over states
        for (std::size_t state_index=0; state_index<subspace.size(); ++state_index)
          {
            OrbitalStatePN state(subspace,state_index);
            jmax = std::max(jmax,state.j());
          }
      }
    int Jmax = TwiceValue(jmax);

    // iterate over two-body species
    std::array<basis::TwoBodySpeciesPN,3> species_list;
    if (space_ordering_==basis::TwoBodySpaceJJJPNOrdering::kPN)
      species_list = {TwoBodySpeciesPN::kPP,TwoBodySpeciesPN::kNN,TwoBodySpeciesPN::kPN};
    else  // if (space_ordering_==basis::TwoBodySpaceJJJPNOrdering::kTz)
      species_list = {TwoBodySpeciesPN::kPP,TwoBodySpeciesPN::kPN,TwoBodySpeciesPN::kNN};
    for (TwoBodySpeciesPN two_body_species : species_list)
      // iterate over J
      for (int J=0; J<=Jmax; ++J)
      {
        // iterate over g
        for (int g=0; g<=1; ++g)
          {
            TwoBodySubspaceJJJPN subspace(
                orbital_space,
                two_body_species,J,g,
                weight_max
              );

            if (subspace.size()!=0)
              PushSubspace(subspace);
          }
      }
  }


  std::string TwoBodySpaceJJJPN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "sJg"
          << " " << std::setw(width) << int(subspace.two_body_species())
          << " " << std::setw(width) << subspace.J()
          << " " << std::setw(width) << subspace.g()
          << " " << "dim"
          << " " << std::setw(width) << subspace.size()
          << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJPN::TwoBodySectorsJJJPN(
      const TwoBodySpaceJJJPN& space,
      int J0, int g0, int Tz0,
      basis::SectorDirection sector_direction
    )
    : BaseSectors(space), J0_(J0), g0_(g0), Tz0_(Tz0)
  {

    // enforce canonical
    assert(sector_direction == basis::SectorDirection::kCanonical);
    // enforce nonzero-Tz0 only if {pp,pn,nn} ordering
    assert((Tz0_ == 0) || (space.space_ordering() == basis::TwoBodySpaceJJJPNOrdering::kTz));

    for (std::size_t bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (std::size_t ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {

          // enforce canonical ordering
          if (
              (sector_direction == basis::SectorDirection::kCanonical)
              && !(bra_subspace_index<=ket_subspace_index)
            )
            continue;

          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // enforce particle conservation
          int bra_Tz = bra_subspace.Tz();
          int ket_Tz = ket_subspace.Tz();
          if (!(bra_Tz == Tz0_ + ket_Tz))
              continue;

          // verify angular momentum, isosopin, and parity selection rules
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
          if (allowed)
            PushSector(bra_subspace_index,ket_subspace_index);
        }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
