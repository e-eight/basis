/****************************************************************
  @file nlj_orbital.h

  Defines general single-particle orbital sets.

  Language: C++11

  Mark A. Caprio, Patrick J. Fasano
  University of Notre Dame

  + 7/7/16 (mac): Created (jjjpnorb_scheme), building on code from jjjt_scheme.
  + 7/19/16 (mac):
    - Add default constructors.
    - Use enum Rank for truncation rank.
    - Add two-body species code definitions.
    - Add GetOrbital accessors.
  + 7/22/16 (mac):
    - Fix reference error in TwoBodySpaceJJJPN.
    - Add debugging strings.
  + 9/28/16 (mac,pjf): Break out into nlj_orbital.
  + 10/6/16 (pjf): Add LJPN classes.
  + 10/7/16 (pjf): Add LJPN sectors and general constructors.
  + 10/13/16 (mac): Add default constructors.
  + 10/18/16 (pjf): Add PN general constructors and serializers (OrbitalInfo).
  + 10/18/16 (pjf):
    - Add general constructors and serializers for Orbital*PN (OrbitalPNInfo).
    - Refactor OrbitalDefinitionStr out of spaces.
    - Add oscillator-likeness tests to OrbitalSpacePN and OrbitalSubspacePN.
    - Add == operator for OrbitalPNInfo.
  + 10/19/16 (pjf): Update documentation to Doxygen format.
  + 10/19/16 (mac):
    - Add default constructors yet again.
    - Make orbital I/O switchable between standalone and embedded modes.
  + 10/22/16 (mac): Augment debug strings.
  + 10/25/16 (pjf):
    - Fix constrained sector enumeration to use l0max and Tz0.
    - Add sector enumeration between different bra and ket spaces.
    - Define array mappings for orbital species.
    - Define Tz accessor for states and subspaces.
  + 10/26/16 (mac): Add equality test for orbitals.
  + 10/26/16 (pjf): Add stream operators to OrbitalPNInfo.
  + 10/26/16 (mac): Add OrbitalStatePN::LabelStr().
  + 11/5/16 (mac): Update MFDn orbital file version code (and define enum).
  + 11/13/16 (mac): Fix sorting order of orbitals.
  + 11/21/16 (mac): Add FullLabels accessor for orbitals.
  + 11/24/16 (pjf):
    - Remove sorting from OrbitalDefinitionStr().
    - Remove operator< from OrbitalPNInfo.
  + 09/20/17 (pjf): Fix Tz0 constraint in OrbitalSectorsLJPN
  + 10/13/17 (pjf): Fix j0 type in OrbitalSectorsLJPN.
  + 10/15/17 (pjf): Copy Nmax sanity checks to OrbitalSubspaceLJPN and OrbitalSpaceLJPN.
  + 02/11/18 (pjf): Remove radial sector constraint mode.
  + 09/07/18 (pjf):
    - Insert space before header lines for MFDn compatibility.
    - Add OrbitalPNList typedef.
  + 02/01/19 (pjf):
    - Mark OrbitalSpacePN and OrbitalSpaceLJPN constructors explicit as appropriate.
    - Add constructor for converting from OrbitalSpacePN to OrbitalSpaceLJPN.
  + 02/19/19 (pjf): Add v15200 orbital format.
  + 04/03/19 (pjf): Define type conversion between OrbitalPNInfo and
      FullOrbitalLabels.
  + 04/07/19 (pjf): Add LabelStr() for OrbitalStateLJPN.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 05/27/19 (pjf): Update to initialize BaseSectors with spaces.
  + 06/03/19 (pjf): Add OrbitalPNInfo constructor which takes tz.
  + 09/06/19 (pjf): Add TruncateOrbitalList.
  + 10/10/19 (pjf): Fix standalone output of 15200 orbital files.

****************************************************************/

#ifndef BASIS_NLJ_ORBITAL_H_
#define BASIS_NLJ_ORBITAL_H_

#include <cstddef>
#include <array>
#include <string>

#include "am/halfint.h"
#include "mcutils/parsing.h"

#include "basis/basis.h"
#include "basis/many_body.h"
#include "basis/proton_neutron.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // orbital label containers
  ////////////////////////////////////////////////////////////////

  typedef std::tuple<OrbitalSpeciesPN,int,int,HalfInt> FullOrbitalLabels;
  ///< Full (species,n,l,j) labels for an orbital.

  struct OrbitalPNInfo
  ///
  /// A flattened orbital container.
  ///
  /// All quantum numbers used by MFDn are contained in this simple struct.
  ///
  {
    OrbitalSpeciesPN orbital_species;
    int n;
    int l;
    HalfInt j;
    double weight;

    OrbitalPNInfo() = default;
    // default constructor

    OrbitalPNInfo(OrbitalSpeciesPN os, int n, int l, HalfInt j, double weight)
      : orbital_species(os), n(n), l(l), j(j), weight(weight) {};

    OrbitalPNInfo(int n, int l, HalfInt j, HalfInt tz, double weight)
      : orbital_species(kTzCodeOrbitalSpeciesPN.at(tz)), n(n), l(l), j(j), weight(weight)
    {}

    inline operator FullOrbitalLabels() const
    {
      return FullOrbitalLabels(orbital_species, n, l, j);
    }
    inline FullOrbitalLabels Key() const {
      return static_cast<FullOrbitalLabels>(*this);
    }
    static constexpr double kWeightTolerance = 1e-8;
    inline bool operator==(const OrbitalPNInfo& rhs) const
    {
      bool equiv = true;
      equiv &= (this->Key()==rhs.Key());
      equiv &= (fabs(this->weight-rhs.weight)<kWeightTolerance);
      return equiv;
    }

    // friend std::ostream& operator<<(std::ostream &out, const OrbitalPNInfo& orbital_info);
    // friend std::istream& operator>>(std::istream &in, OrbitalPNInfo& orbital_info);
  };

  // a simple list of orbitals
  typedef std::vector<OrbitalPNInfo> OrbitalPNList;

  // canonical orbital sorting
  inline bool OrbitalSortCmpWeight(const OrbitalPNInfo& lhs, const OrbitalPNInfo& rhs) {
    std::tuple<OrbitalSpeciesPN,double,int,HalfInt,int>
      lhs_labels(lhs.orbital_species, lhs.weight, lhs.l, lhs.j, lhs.n);
    std::tuple<OrbitalSpeciesPN,double,int,HalfInt,int>
      rhs_labels(rhs.orbital_species, rhs.weight, rhs.l, rhs.j, rhs.n);
    return (lhs_labels < rhs_labels);
  }

  OrbitalPNList TruncateOrbitalList(
      const WeightMax& weight_max, const OrbitalPNList& orbital_list
    );
  // Apply weight max truncation to a list of orbitals.
  //
  // Arguments:
  //   weight_max (WeightMax): truncation parameters
  //   orbital_list (OrbitalPNList): input list of orbitals

  // orbital I/O

  // orbital file version codes
  enum class MFDnOrbitalFormat : int {kVersion15099=15099, kVersion15200=15200};

  OrbitalPNList ParseOrbitalPNStream(
      std::istream& is,
      bool standalone,
      MFDnOrbitalFormat format = MFDnOrbitalFormat::kVersion15099
    );
  /// Read orbital definitions from a stream.
  ///
  /// Arguments:
  ///   is (std::istream, input) :
  ///     input stream containing MFDn-formatted orbital definitions
  ///   standalone (bool): whether or not to expect initial comments and version number
  ///     as for standalone orbital file
  ///   format (MFDnOrbitalFormat, optional): orbital list format
  ///     - (ignored/deduced for standalone)
  ///
  /// Returns:
  ///   (OrbitalPNList) : list of flattened orbital parameters

  std::string OrbitalDefinitionStr(
      const OrbitalPNList& orbitals,
      bool standalone = false,
      MFDnOrbitalFormat format = MFDnOrbitalFormat::kVersion15099
    );
  /// Output orbital info as a string suitable for MFDn version 15.
  ///
  /// Arguments:
  ///   orbitals (const OrbitalPNList&, input) :
  ///     list of flattened orbital parameters
  ///   standalone (bool, optional): whether or not to include initial comments and version number
  ///     as for standalone orbital file
  ///   format (MFDnOrbitalFormat, optional): orbital list format
  ///
  /// Returns:
  ///   (std::string) output stream containing MFDn-formatted orbital definitions

  ////////////////////////////////////////////////////////////////
  /// @defgroup pn-subspaces PN Subspaces
  /// Single particle orbitals divided into subspaces with definite
  /// orbital species.
  ///
  /// ## Labeling ##
  ///
  /// subspace labels: (species)
  ///
  ///  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  ///
  /// state labels within subspace: (n,l,j)
  ///
  ///  * n (int): radial quantum number (0,1,...)
  ///  * l (int): orbital angular momentum
  ///  * j (HalfInt): total angular momentum
  ///
  /// The parity for each orbital is deduced from the l quantum number:
  ///
  ///  * g (int): grade (=0,1) for the parity P, given by g~l
  ///
  /// Each orbital (state) also has a floating point "weight"
  /// associated with it:
  ///
  ///  * weight (double): orbital weight
  ///
  /// This quantity is not considered a "label", since it is not used
  /// for lookup purposes.  In a traditional oscillator scheme, the
  /// weight is set to the oscillator quantum number w -> N=2n+l.
  ///
  /// The hard-coded oscillator quantum number is also deduced from the
  /// n and l quantum numbers, as an integer, to be used only when the
  /// orbitals are known to be oscillator orbitals:
  ///
  ///  * N (int): oscillator quanta (N=2n+l)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## Subspaces ##
  ///
  /// Within the full space, subspaces are ordered by:
  ///    * enumerated species
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## States ##
  ///
  /// Within a subspace, the states are ordered by:
  ///   * lexicographically increasing (n,l,j)
  ///
  ////////////////////////////////////////////////////////////////
  /// @{

  // labels

  typedef std::tuple<OrbitalSpeciesPN> OrbitalSubspacePNLabels;
  typedef std::tuple<int,int,HalfInt> OrbitalStatePNLabels;


  // subspace

  class OrbitalSubspacePN
    : public BaseSubspace<OrbitalSubspacePNLabels,OrbitalStatePNLabels>
    {

      public:

      // constructor

      OrbitalSubspacePN() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax);
      // Set up indexing and weights in traditional oscillator Nmax
      // truncation.

      OrbitalSubspacePN(OrbitalSpeciesPN orbital_species,
        const OrbitalPNList& states);
      // Set up indexing for a list of states.

      // produce flattened orbital information
      OrbitalPNList OrbitalInfo() const;

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
      HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
      double weight_max() const {return weight_max_;}
      bool is_oscillator_like() const {return is_oscillator_like_;}
      int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
      // only meaningful if oscillator scheme constructor used
      const std::vector<double>& weights() const {return weights_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      // truncation
      double weight_max_;
      bool is_oscillator_like_;
      int Nmax_;  // only meaningful if oscillator scheme constructor used

      // test for truncation scheme
      bool IsOscillatorLike_() const;
      // Test if labeling and weights match oscillator truncation.

      // weights
      std::vector<double> weights_;
    };

  // state

  class OrbitalStatePN
    : public BaseState<OrbitalSubspacePN>
  {

    public:

    // pass-through constructors

    OrbitalStatePN(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStatePN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // produce flattened orbital information
    OrbitalPNInfo OrbitalInfo() const;

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
    HalfInt Tz() const {return subspace().Tz();}

    // state label accessors
    int n() const {return std::get<0>(labels());}
    int l() const {return std::get<1>(labels());}
    HalfInt j() const {return std::get<2>(labels());}
    int g() const {return l()%2;}
    FullOrbitalLabels full_labels() const
    {
      return std::make_tuple(orbital_species(),n(),l(),j());
    }

    // state weight accessors
    double weight() const {return subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.

    // comparison
    friend bool operator == (const basis::OrbitalStatePN& a1, const basis::OrbitalStatePN& a2)
    // Equality test based on labels (so permits comparison across different subspace indexings).
    {
      return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
    }

  };

  // space

  class OrbitalSpacePN
    : public BaseSpace<OrbitalSubspacePN>
  {

    public:

    // constructor

    OrbitalSpacePN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OrbitalSpacePN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    explicit OrbitalSpacePN(const OrbitalPNList& states);
    // Set up indexing for a list of states.

    // produce flattened orbital information
    OrbitalPNList OrbitalInfo() const;

    // accessors
    double weight_max() const {return weight_max_;}
    bool is_oscillator_like() const {return is_oscillator_like_;}
    int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
    // only meaningful if oscillator scheme constructor used

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    double weight_max_;
    bool is_oscillator_like_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

    // test for truncation scheme
    bool IsOscillatorLike_() const;
    // Test if labeling and weights match oscillator truncation for
    // all subspaces.

  };

  // sectors -- not applicable

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////
  /// @defgroup ljpn-subspaces LJPN Subspaces
  /// Single particle orbitals divided into subspaces with definite l, j and
  /// orbital species.
  ///
  /// ## Labeling ##
  ///
  /// subspace labels: (species,l,j)
  ///
  ///  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  ///  * l (int): orbital angular momentum
  ///  * j (HalfInt): total angular momentum
  ///
  /// state labels within subspace: (n,l,j)
  ///
  ///  * n (int): radial quantum number (0,1,...)
  ///
  /// The parity for each orbital is deduced from the l quantum number:
  ///
  ///  * g (int): grade (=0,1) for the parity P, given by g~l
  ///
  /// Each orbital (state) also has a floating point "weight"
  /// associated with it:
  ///
  ///  * weight (double): orbital weight
  ///
  /// This quantity is not considered a "label", since it is not used
  /// for lookup purposes.  In a traditional oscillator scheme, the
  /// weight is set to the oscillator quantum number w -> N=2n+l.
  ///
  /// The hard-coded oscillator quantum number is also deduced from the
  /// n and l quantum numbers, as an integer, to be used only when the
  /// orbitals are known to be oscillator orbitals:
  ///
  ///   *N (int): oscillator quanta (N=2n+l)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## Subspaces ##
  ///
  /// Within the full space, subspaces are ordered by:
  ///    * lexicographically increasing (species,l,j)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## States ##
  ///
  /// Within a subspace, the states are ordered by:
  ///   * increasing n
  ///
  ////////////////////////////////////////////////////////////////
  /// @{

  // labels

  typedef std::tuple<OrbitalSpeciesPN,int,HalfInt> OrbitalSubspaceLJPNLabels;
  typedef std::tuple<int> OrbitalStateLJPNLabels;

  // subspace

  class OrbitalSubspaceLJPN
    : public BaseSubspace<OrbitalSubspaceLJPNLabels,OrbitalStateLJPNLabels>
    {

      public:

      // constructor

      OrbitalSubspaceLJPN() = default;
      /// default constructor -- provided since required for certain
      /// purposes by STL container classes (e.g., std::vector::resize)

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j, int Nmax);
      /// Set up indexing and weights in traditional oscillator Nmax
      /// truncation.

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j,
        const OrbitalPNList& states);
      /// Set up indexing for a list of states.

      OrbitalPNList OrbitalInfo() const;
      /// produce flattened orbital information

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
      HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
      int l() const {return std::get<1>(labels_);}
      HalfInt j() const {return std::get<2>(labels_);}
      int g() const {return l()%2;}
      double weight_max() const {return weight_max_;}
      bool is_oscillator_like() const {return is_oscillator_like_;}
      int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
        // only meaningful if oscillator scheme constructor used
      const std::vector<double>& weights() const {return weights_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      // truncation
      double weight_max_;
      bool is_oscillator_like_;
      int Nmax_;  // only meaningful if oscillator scheme constructor used

      // test for truncation scheme
      bool IsOscillatorLike_() const;
      // Test if labeling and weights match oscillator truncation.

      // weights
      std::vector<double> weights_;
    };

  // state

  class OrbitalStateLJPN
    : public BaseState<OrbitalSubspaceLJPN>
  {

    public:

    // pass-through constructors

    OrbitalStateLJPN(const SubspaceType& subspace, std::size_t index)
      /// Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStateLJPN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      /// Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    OrbitalPNInfo OrbitalInfo() const;
    /// produce flattened orbital information

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
    HalfInt Tz() const {return subspace().Tz();}
    int l() const {return subspace().l();}
    HalfInt j() const {return subspace().j();}
    int g() const {return l()%2;}
    FullOrbitalLabels full_labels() const
    {
      return std::make_tuple(orbital_species(),n(),l(),j());
    }

    // state label accessors
    int n() const {return std::get<0>(labels());}

    // state weight accessors
    double weight() const {return subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.
  };

  // space

  class OrbitalSpaceLJPN
    : public BaseSpace<OrbitalSubspaceLJPN>
  {

    public:

    // constructor

    OrbitalSpaceLJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OrbitalSpaceLJPN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    explicit OrbitalSpaceLJPN(const OrbitalPNList& states);
    // Set up indexing for a list of states.

    explicit OrbitalSpaceLJPN(const OrbitalSpacePN& space)
      : OrbitalSpaceLJPN(space.OrbitalInfo())
    {}
    // Set up indexing from space divided into pn-sectors.

    // accessors
    double weight_max() const {return weight_max_;}
    bool is_oscillator_like() const {return is_oscillator_like_;}
    int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
      // only meaningful if oscillator scheme constructor used

    // produce flattened orbital information
    OrbitalPNList OrbitalInfo() const;

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    double weight_max_;
    bool is_oscillator_like_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

    // test for truncation scheme
    bool IsOscillatorLike_() const;
    // Test if labeling and weights match oscillator truncation for
    // all subspaces.
  };

  // sectors

  class OrbitalSectorsLJPN
    : public BaseSectors<OrbitalSpaceLJPN>
  {

    public:

    // constructor

    OrbitalSectorsLJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    DEPRECATED("all-to-all enumeration is now deprecated")
    OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space
    );
    DEPRECATED("all-to-all enumeration is now deprecated")
    explicit OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space
    ) : OrbitalSectorsLJPN(space, space) {}
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).

    OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space,
      int j0, int g0, int Tz0
    );
    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space,
        int j0, int g0, int Tz0
    ) : OrbitalSectorsLJPN(space, space, j0, g0, Tz0) {}
    // Enumerate sector pairs between two spaces connected by an operator of
    // given tensorial and parity character ("constrained" sector enumeration).

    // diagnostic string
    std::string DebugStr() const;

    // accessors
    int j0()     const {return j0_;}
    int g0()     const {return g0_;}
    int Tz0()    const {return Tz0_;}

   private:
    // operator properties
    int j0_, g0_, Tz0_;
  };

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // NLJ_ORBITAL_H_
