/************************************************************//**
  @file jjjpn_scheme.h

  Defines two-body state indexing in jjJpn coupling scheme, based on
  general single-particle orbital sets.

  Nomenclature: The "orb" in the filename flags that these definitions
  are for general orbitals, rather than oscillator orbitals (compare,
  e.g., lsjt_scheme.h or jjjt_scheme.h).  However, we do not propagate
  this "orb" nomenclature to the actual typedefs and function names
  defined in this file.  Naming conflicts would therefore arise if we
  were to wish to define a traditional "jjjpn_scheme" with hard-coded
  (N,j) orbital labels, in the spirit of jjjt_scheme.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 7/7/16 (mac): Created, building on code from jjjt_scheme (jjjpnorb_scheme).
  + 7/19/16 (mac):
    - Add default constructors.
    - Use enum Rank for truncation rank.
    - Add two-body species code definitions.
    - Add GetOrbital accessors.
  + 7/22/16 (mac):
    - Fix reference error in TwoBodySpaceJJJPN.
    - Add debugging strings.
  + 9/28/16 (mac,pjf): Extract orbital definitions into nlj_orbital.
  + 10/14/16 (mac): Add constructors for WeightMax.
  + 10/14/16 (mac): Store operator properties with TwoBodySectorsJJJPN.
  + 10/25/16 (mac):
    - Rename to jjjpn_scheme.
    - Add Tz0 argument to sectors constructor.
    - Add reference to orbital_space from space (currently disabled due to
      initializer issues).
  + 7/1/17 (mac): Extract generic proton-neutron definitions to
    proton_neutron.
  + 1/22/18 (mac): Enable nonzero Tz0 in sector enumeration.
  + 02/12/19 (pjf): Allow space ordering {pp,pn,nn} for h2v15200.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 09/06/19 (pjf): Move WeightMax to many_body.h.
  + 12/12/19 (pjf): Fix UB in TwoBodySpaceJJJPN construction.
****************************************************************/

#ifndef BASIS_JJJPN_SCHEME_H_
#define BASIS_JJJPN_SCHEME_H_

#include <cstddef>
#include <array>
#include <string>

#include "am/halfint.h"

#include "basis/basis.h"
#include "basis/many_body.h"
#include "basis/proton_neutron.h"
#include "basis/nlj_orbital.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (species,J,g)    P=(-)^g
  //
  //   species (enum): two-body species (equivalent to Tz)
  //   J (int): total angular momentum
  //   g (int): grade (=0,1) for the parity P
  //
  // state labels within subspace: (index1,index2)
  //
  //   index1 (std::size_t): index of particle 1 within appropriate
  //     (proton or neutron) orbital set
  //   index2 (std::size_t): index of particle 2 within appropriate
  //     (proton or neutron) orbital set
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing species (species=pp,nn,pn)
  //    -- increasing J
  //    -- increasing g (g=0,1)
  //
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (species,J,g).
  //
  // Truncation of the space is by one-body and two-body weights.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing index1
  //   -- increasing index2
  // and subject to:
  //   -- triangularity constraint on (j1,j2,J)
  //   -- parity constraint g1+g2~g
  //   -- in pp/nn subspaces, antisymmetry constraint J~0
  //      if index1==index2
  //
  // This basis is for *identical* particle states:
  //   -- In the pp/nn subspaces, the labels are subject to
  //      the antisymmetry constraint (J~0) if the orbitals are identical.
  //   -- In the pp/nn subspaces, a canonical (lexicographic) ordering
  //      constraint is applied to the single-particle quantum numbers.
  //      That is, when enumerating the basis, the states
  //
  //        |(index1,index2)...>  and  |(index2,index1)...>
  //
  //      would be redundant, and only the first (for index1<=index2) is
  //      retained.
  //
  ////////////////////////////////////////////////////////////////

  // labels

  typedef std::tuple<TwoBodySpeciesPN,int,int> TwoBodySubspaceJJJPNLabels;
  typedef std::tuple<std::size_t,std::size_t> TwoBodyStateJJJPNLabels;

  // subspace

  class TwoBodySubspaceJJJPN
    : public BaseSubspace<TwoBodySubspaceJJJPNLabels,TwoBodyStateJJJPNLabels>
    {

      public:

      // constructor

      TwoBodySubspaceJJJPN() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceJJJPN(
          const OrbitalSpacePN& orbital_space,
          TwoBodySpeciesPN two_body_species, int J, int g,
          const WeightMax& weight_max
        );
      // Set up indexing.

      // accessors
      TwoBodySpeciesPN two_body_species() const {return std::get<0>(labels_);}
      int J() const {return std::get<1>(labels_);}
      int g() const {return std::get<2>(labels_);}
      int Tz() const {return basis::kTwoBodySpeciesPNCodeTz[int(two_body_species())];}
      const WeightMax& weight_max() const {return weight_max_;}
      const OrbitalSubspacePN& orbital_subspace1() const {return *orbital_subspace1_ptr_;}
      const OrbitalSubspacePN& orbital_subspace2() const {return *orbital_subspace2_ptr_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      // truncation
      WeightMax weight_max_;

      // direct access to orbital subspaces
      const OrbitalSubspacePN* orbital_subspace1_ptr_;
      const OrbitalSubspacePN* orbital_subspace2_ptr_;
    };

  // state

  class TwoBodyStateJJJPN
    : public BaseState<TwoBodySubspaceJJJPN>
  {

    public:

    // pass-through constructors

    TwoBodyStateJJJPN(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateJJJPN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    TwoBodySpeciesPN two_body_species() const {return subspace().two_body_species();}
    int J() const {return subspace().J();}
    int g() const {return subspace().g();}
    const OrbitalSubspacePN& orbital_subspace1() const {return subspace().orbital_subspace1();}
    const OrbitalSubspacePN& orbital_subspace2() const {return subspace().orbital_subspace2();}

    // state label accessors
    std::size_t index1() const {return std::get<0>(labels());}
    std::size_t index2() const {return std::get<1>(labels());}

    // state retrieval
    const OrbitalStatePN GetOrbital1() const {return OrbitalStatePN(orbital_subspace1(),index1());}
    const OrbitalStatePN GetOrbital2() const {return OrbitalStatePN(orbital_subspace2(),index2());}

    // diagnostic string
    std::string LabelStr() const;

  };

  // subspace ordering in space
  //
  // kPN -> {pp,nn,pn}
  // kTz -> {pp,pn,nn}
  enum class TwoBodySpaceJJJPNOrdering : int { kPN=0, kTz=1 };

  // space

  class TwoBodySpaceJJJPN
    : public BaseSpace<TwoBodySubspaceJJJPN>
  {

    public:

    // constructor

    TwoBodySpaceJJJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceJJJPN(
          const OrbitalSpacePN& orbital_space,
          const WeightMax& weight_max,
          basis::TwoBodySpaceJJJPNOrdering ordering = TwoBodySpaceJJJPNOrdering::kPN
      );
    // Enumerate subspaces.

    // accessors
    // const OrbitalSpacePN& orbital_space() {return orbital_space_;}
    const WeightMax& weight_max() const {return weight_max_;}
    const TwoBodySpaceJJJPNOrdering& space_ordering() const {return space_ordering_;}

    // diagnostic string
    std::string DebugStr() const;

    private:

    // convenience reference to underlying orbitals
    //
    // Caveat: Any reference member interferes with defining a default
    // constructor, since references must be explicitly initialized.

    // const OrbitalSpacePN& orbital_space_;

    // truncation
    WeightMax weight_max_;
    TwoBodySpaceJJJPNOrdering space_ordering_;

  };

  // sectors

  class TwoBodySectorsJJJPN
    : public BaseSectors<TwoBodySpaceJJJPN>
  {

    public:

    // constructor

    TwoBodySectorsJJJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsJJJPN(
        const TwoBodySpaceJJJPN& space,
        int J0, int g0, int Tz0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

    // accessors
    int J0() const {return J0_;};
    int g0() const {return g0_;};
    int Tz0() const {return Tz0_;};

    private:

    // operator properties
    int J0_, g0_, Tz0_;
  };



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
