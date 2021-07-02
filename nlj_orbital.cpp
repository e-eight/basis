/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio, Patrick J. Fasano
  University of Notre Dame

****************************************************************/


#include <cstddef>
#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>
#include <set>
#include <algorithm>

#include "am/am.h"
#include "mcutils/parsing.h"

#include "nlj_orbital.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // orbital truncation
  ////////////////////////////////////////////////////////////////

  OrbitalPNList TruncateOrbitalList(
      const WeightMax& weight_max, const OrbitalPNList& orbital_list
    )
  {
    OrbitalPNList output_orbitals;

    for (const auto& orbital : orbital_list)
    {
      if (orbital.weight <= weight_max.one_body[int(orbital.orbital_species)])
        output_orbitals.push_back(orbital);
    }

    return output_orbitals;
  }

  ////////////////////////////////////////////////////////////////
  // single-particle definition file parsing
  ////////////////////////////////////////////////////////////////

  /**
   * Read orbital definitions from a stream.
   *
   * @param[in] is input stream containing MFDn-formatted orbital definitions
   * @return list of flattened orbital parameters
   */
  OrbitalPNList ParseOrbitalPNStream(
      std::istream& is,
      bool standalone,
      MFDnOrbitalFormat format
    )
  {

    // set up line counter for use in error messages
    std::string line;
    int line_count = 0;


    int num_orbitals_p, num_orbitals_n;
    if (standalone)
      {
        // line 1: version
        {
          mcutils::GetLine(is,line,line_count);
          std::istringstream line_stream(line);
          int version;
          line_stream >> version;
          mcutils::ParsingCheck(line_stream,line_count,line);
          format = static_cast<MFDnOrbitalFormat>(version);
          assert((format==MFDnOrbitalFormat::kVersion15099)
                 ||(format==MFDnOrbitalFormat::kVersion15200));
        }

        // line 2: number of p,n orbitals
        {
          mcutils::GetLine(is,line,line_count);
          std::istringstream line_stream(line);
          line_stream >> num_orbitals_p >> num_orbitals_n;
          mcutils::ParsingCheck(line_stream,line_count,line);
        }
      }

    // lines 3+: orbital definitions
    OrbitalPNList states;
    int num_orbitals_p_extracted=0, num_orbitals_n_extracted=0;
    while (mcutils::GetLine(is,line,line_count))
      {
        // set up for parsing
        std::istringstream line_stream(line);

        int index;
        int n, l, twice_j;
        basis::OrbitalSpeciesPN species;
        double weight;
        if (format==MFDnOrbitalFormat::kVersion15099)
          {
            int species_code;
            line_stream >> index >> n >> l >> twice_j >> species_code >> weight;
            species = kDecimalCodeOrbitalSpeciesPN.at(species_code);
          }
        else if (format==MFDnOrbitalFormat::kVersion15200)
          {
            int twice_tz;
            line_stream >> index >> n >> l >> twice_j >> twice_tz >> weight;
            species = kTzCodeOrbitalSpeciesPN.at(HalfInt(twice_tz,2));
          }
        mcutils::ParsingCheck(line_stream,line_count,line);

        OrbitalPNInfo state(species, n, l, HalfInt(twice_j, 2), weight);
        // count orbitals by type
        num_orbitals_p_extracted +=
          static_cast<int>(state.orbital_species == OrbitalSpeciesPN::kP);
        num_orbitals_n_extracted +=
          static_cast<int>(state.orbital_species == OrbitalSpeciesPN::kN);

        states.push_back(state);
      }

    if (standalone)
      {
        assert(num_orbitals_p==num_orbitals_p_extracted);
        assert(num_orbitals_n==num_orbitals_n_extracted);
      }

    return states;
  }

  /**
   * Output orbital info as a string suitable for MFDn version 15.
   *
   * @param[in] orbitals list of flattened orbital parameters
   * @return output stream containing MFDn-formatted orbital definitions
   */
  std::string OrbitalDefinitionStr(
      const OrbitalPNList& orbitals,
      bool standalone,
      MFDnOrbitalFormat format
    )
  {
    // sanity check
    assert((format==MFDnOrbitalFormat::kVersion15099)
           ||(format==MFDnOrbitalFormat::kVersion15200));

    std::ostringstream header;
    std::ostringstream body;
    std::ostringstream os;

    // construct body
    const int width = 3;
    const int precision = 8;
    // orbital indices
    int p_index = 0;
    int n_index = 0;
    int output_index = 0;
    for (const OrbitalPNInfo& orbital_info: orbitals)
      // iterate over states
      {
        if (orbital_info.orbital_species == OrbitalSpeciesPN::kP) {
          output_index = ++p_index;
        } else if (orbital_info.orbital_species == OrbitalSpeciesPN::kN) {
          output_index = ++n_index;
        }
        if (format==MFDnOrbitalFormat::kVersion15099)
          {
            body << " " << std::setw(width) << output_index
                 << " " << std::setw(width) << orbital_info.n
                 << " " << std::setw(width) << orbital_info.l
                 << " " << std::setw(width) << TwiceValue(orbital_info.j)
                 << " " << std::setw(width) << kOrbitalSpeciesPNCodeDecimal[static_cast<int>(orbital_info.orbital_species)]  // 1-based
                 << " " << std::fixed << std::setw(width+1+precision)
                 << std::setprecision(precision) << orbital_info.weight
                 << std::endl;
          }
        else if (format==MFDnOrbitalFormat::kVersion15200)
          {
            output_index = p_index + n_index;
            body << " " << std::setw(width) << output_index
                 << " " << std::setw(width) << orbital_info.n
                 << " " << std::setw(width) << orbital_info.l
                 << " " << std::setw(width) << TwiceValue(orbital_info.j)
                 << " " << std::setw(width) << kOrbitalSpeciesPNCodeTz[static_cast<int>(orbital_info.orbital_species)].TwiceValue()
                 << " " << std::fixed << std::setw(width+1+precision)
                 << std::setprecision(precision) << orbital_info.weight
                 << std::endl;
          }
      }

    // construct header
    //
    // We defer constructing the header until after constructing the
    // body, since we need statistics obtained while constructing the
    // body.
    if (standalone)
      {
        // header comments
        header << "# MFDn SPorbital file" << std::endl;
        header << "#   version" << std::endl;
        header << "#   norb_p norb_n" << std::endl;
        header << "#   index n l 2*j species weight" << std::endl;

        // header line 1: version
        int version = static_cast<int>(format);
        header << " " << version << std::endl;

        // header line 2: dimensions
        header << " " << p_index << " " << n_index << std::endl;
      }

    // assemble file
    os << header.str() << body.str();

    return os.str();

  }

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  /**
   * Construct an Nmax-truncated single-particle subspace with a particular
   * species.
   *
   * @param[in] orbital_species species type for subspace
   * @param[in] Nmax number of oscillator quanta
   */
  OrbitalSubspacePN::OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(orbital_species);
    weight_max_ = double(Nmax);
    is_oscillator_like_ = true;
    Nmax_ = Nmax;

    // iterate over total oscillator quanta
    for (int N = 0; N <= Nmax; ++N)
      // iterate over j within shell
      for (HalfInt j = HalfInt(1,2); j <= N+HalfInt(1,2); ++j)
        {
          // recover derived quantum numbers (n,l) from (N,j)
          int l = (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
          int n = (N-l)/2;

          // save state
          PushStateLabels(StateLabelsType(n,l,j));

          // save oscillator quantum number as weight
          weights_.push_back(double(N));
        }
  }

  /**
   * Construct a subspace with a particular species from a list of orbitals.
   *
   * @param[in] orbital_species species type for subspace
   * @param[in] states vector of orbitals
   */
  OrbitalSubspacePN::OrbitalSubspacePN(OrbitalSpeciesPN orbital_species,
                                       const OrbitalPNList& states)
  {

    // set values
    labels_ = SubspaceLabelsType(orbital_species);
    weight_max_ = 0.0;

    // iterate over all states, picking out those which belong to this subspace
    for (const OrbitalPNInfo& state : states) {
      if (state.orbital_species == orbital_species) {
        PushStateLabels(StateLabelsType(state.n,state.l,state.j));
        weights_.push_back(state.weight);
        weight_max_ = std::max(weight_max_,state.weight);
      }
    }

    // check if oscillator-like and set Nmax
    if (IsOscillatorLike_()) {
      is_oscillator_like_ = true;
      Nmax_ = int(weight_max_);
    } else {
      is_oscillator_like_ = false;
      Nmax_ = -1;
    }

  }

  /**
   * Do a deep comparison to oscillator-truncated basis.
   * @return true if truncation is like Nmax-truncated oscillator
   */
  bool OrbitalSubspacePN::IsOscillatorLike_() const
  {
    // oscillators have states
    if (size()==0)
      return false;

    // maximum weight should be an integer
    int Nmax = int(weight_max());
    if (Nmax != weight_max())
      return false;

    // only positive Nmax is allowed
    if (Nmax < 0)
      return false;

    // compare to equivalent Nmax-truncated subspace
    OrbitalSubspacePN reference_subspace(orbital_species(),Nmax);
    if (reference_subspace.OrbitalInfo() != OrbitalInfo())
      return false;

    // if it looks like an oscillator, swims like an oscillator, and quacks
    // like an oscillator, it's probably like an oscillator
    return true;
  }

  /**
   * Generate a string representation of the subspace labels.
   * @return subspace labels as a string
   */
  std::string OrbitalSubspacePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << "]";

    return os.str();
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string OrbitalSubspacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    std::string oscillator_like_indicator = (is_oscillator_like() ? "true" : "false");
    os << " weight_max " << weight_max();
    if (is_oscillator_like())
      os << " Nmax " << Nmax();
    os << " (oscillator-like: " << oscillator_like_indicator << ")"
       << std::endl;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStatePN state(*this,state_index);

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "nlj"
          << " " << std::setw(width) << state.n()
          << " " << std::setw(width) << state.l()
          << " " << std::setw(width+2) << state.j().Str()
          << " " << "weight"
          << " " << state.weight()
          << std::endl;
      }

    return os.str();

  }

  /**
   * Flatten subspace into a vector of OrbitalPNInfo objects.
   *
   * @return vector representation of subspace
   */
  OrbitalPNList OrbitalSubspacePN::OrbitalInfo() const
  {
    OrbitalPNList orbitals;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStatePN state(*this,state_index);
        orbitals.push_back(state.OrbitalInfo());
      }

    return orbitals;

  }

  /**
   * Flatten state into an OrbitalPNInfo object.
   *
   * @return OrbitalPNInfo representation of state
   */
  OrbitalPNInfo OrbitalStatePN::OrbitalInfo() const
  {
    OrbitalPNInfo orbital;

    orbital.orbital_species = orbital_species();
    orbital.n = n();
    orbital.l = l();
    orbital.j = j();
    orbital.weight = weight();

    return orbital;
  }

  /**
   * Generate a string representation of the orbital labels.
   * @return orbital labels as a string
   */
  std::string OrbitalStatePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << std::setw(width) << index()
       << " :"
       << " " << std::setw(width) << n()
       << " " << std::setw(width) << l()
       << " " << std::setw(width) << j()
       << " :"
       << " " << std::setw(width) << weight()
       << " " << "]";

    return os.str();
  }

  /**
   * Construct an Nmax-truncated single-particle space with species subspaces.
   *
   * @param[in] Nmax number of oscillator quanta
   */
  OrbitalSpacePN::OrbitalSpacePN(int Nmax)
  {

    // save truncation
    weight_max_ = double(Nmax);
    is_oscillator_like_ = true;
    Nmax_ = Nmax;

    // iterate over species
    for (OrbitalSpeciesPN orbital_species : {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN})
      {
        OrbitalSubspacePN subspace(orbital_species,Nmax);
        PushSubspace(subspace);
      }
  }

  /**
   * Construct a space with species subspaces from a list of orbitals.
   *
   * @param[in] states vector of orbitals
   */
  OrbitalSpacePN::OrbitalSpacePN(const OrbitalPNList& states)
  {
    weight_max_ = 0.0;
    // collect orbital_species subspace labels sorted in canonical order
    std::set<OrbitalSubspacePNLabels> subspace_labels_set;
    for (std::size_t state_index=0; state_index<states.size(); ++state_index)
      {
        OrbitalPNInfo state = states[state_index];
        OrbitalSubspacePNLabels labels(state.orbital_species);
        subspace_labels_set.insert(labels);
      }

    // construct subspaces
    for (const OrbitalSubspacePNLabels& labels : subspace_labels_set)
      {
        OrbitalSpeciesPN orbital_species;
        std::tie(orbital_species) = labels;
        OrbitalSubspacePN subspace(orbital_species,states);
        PushSubspace(subspace);
        weight_max_ = std::max(weight_max_,subspace.weight_max());
      }

    // check if oscillator-like and set Nmax
    if (IsOscillatorLike_()) {
      is_oscillator_like_ = true;
      Nmax_ = int(weight_max_);
    } else {
      is_oscillator_like_ = false;
      Nmax_ = -1;
    }

  }

  /**
   * Check if space is truncated like an Nmax oscillator truncation.
   *
   * @return true if all subspaces are Nmax truncated with the same Nmax.
   */
  bool OrbitalSpacePN::IsOscillatorLike_() const
  {

    // first check that the space contains subspaces
    if (size() < 1)
      return false;

    // check if subspaces are individually Nmax truncated
    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        if (!subspace.is_oscillator_like())
          return false;
      }

    // see if can extract viable Nmax
    int Nmax = GetSubspace(0).Nmax();
    for (std::size_t subspace_index=1; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        if (subspace.Nmax() != Nmax)
          return false;
      }

    // orbitals are oscillator like!
    return true;
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string OrbitalSpacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    std::string oscillator_like_indicator = (is_oscillator_like() ? "true" : "false");
    os << " weight_max " << weight_max();
    if (is_oscillator_like())
      os << " Nmax " << Nmax();
    os << " (oscillator-like: " << oscillator_like_indicator << ")"
       << std::endl;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);

        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "species"
          << " " << std::setw(width) << int(subspace.orbital_species())
          << " " << "dim"
          << " " << std::setw(width) << subspace.size()
          << " " << std::endl;
      }

    return os.str();

  }

  /**
   * Flatten space into a vector of OrbitalPNInfo objects.
   *
   * @return vector representation of space
   */
  OrbitalPNList OrbitalSpacePN::OrbitalInfo() const
  {
    OrbitalPNList orbitals;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        OrbitalPNList subspace_orbitals;

        // get orbitals for subspace and append to vector
        subspace_orbitals = subspace.OrbitalInfo();
        orbitals.insert(orbitals.end(),
                        subspace_orbitals.begin(),
                        subspace_orbitals.end()
                       );

      }

    return orbitals;

  }


  ////////////////////////////////////////////////////////////////
  // single-particle orbitals - lj subspaces
  ////////////////////////////////////////////////////////////////

  /**
   * Construct an Nmax-truncated single-particle subspace with a particular
   * species.
   *
   * @param[in] orbital_species species type for subspace
   * @param[in] l orbital angular momentum quantum number
   * @param[in] j total angular momentum quantum number
   * @param[in] Nmax number of oscillator quanta
   */
  OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species,
                                           int l, HalfInt j, int Nmax) {
    // set values
    labels_ = SubspaceLabelsType(orbital_species,l,j);
    weight_max_ = double(Nmax);
    is_oscillator_like_ = true;
    Nmax_ = Nmax;

    // iterate over radial quantum number
    for (int n = 0; (2*n+l) <= Nmax; ++n) {
      // save state
      PushStateLabels(StateLabelsType(n));

      // save oscillator quantum number as weight
      weights_.push_back(double(2*n+l));
    }
  }

  /**
   * Construct a subspace with a particular l, j, and species from a list
   * of orbitals.
   *
   * @param[in] orbital_species species type for subspace
   * @param[in] l orbital angular momentum quantum number
   * @param[in] j total angular momentum quantum number
   * @param[in] states vector of orbitals
   */
  OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species,
                                           int l, HalfInt j,
                                           const OrbitalPNList& states) {
    labels_ = SubspaceLabelsType(orbital_species,l,j);
    weight_max_ = 0.0;
    for (auto&& state : states) {
      if (state.orbital_species == orbital_species
          && state.l == l && state.j == j) {
        PushStateLabels(StateLabelsType(state.n));
        weights_.push_back(state.weight);
        weight_max_ = std::max(weight_max_,state.weight);
      }
    }

    // check if oscillator-like and set Nmax
    if (IsOscillatorLike_()) {
      is_oscillator_like_ = true;
      Nmax_ = static_cast<int>(weight_max_);
    } else {
      is_oscillator_like_ = false;
      Nmax_ = -1;
    }
  }

  /**
   * Do a deep comparison to oscillator-truncated basis.
   * @return true if truncation is like Nmax-truncated oscillator
   */
  bool OrbitalSubspaceLJPN::IsOscillatorLike_() const
  {
    // oscillators have states
    if (size()==0)
      return false;

    // maximum weight should be an integer
    int Nmax = int(weight_max());
    if (Nmax != weight_max())
      return false;

    // only positive Nmax is allowed
    if (Nmax < 0)
      return false;

    // compare to equivalent Nmax-truncated subspace
    OrbitalSubspacePN reference_subspace(orbital_species(),Nmax);
    if (reference_subspace.OrbitalInfo() != OrbitalInfo())
      return false;

    // if it looks like an oscillator, swims like an oscillator, and quacks
    // like an oscillator, it's probably like an oscillator
    return true;
  }

  /**
   * Generate a string representation of the subspace labels.
   * @return subspace labels as a string
   */
  std::string OrbitalSubspaceLJPN::LabelStr() const {
    std::ostringstream os;

    const int width = 3;

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << std::setw(width) << l()
       << " " << std::setw(width+2) << j().Str()
       << " " << "]";

    return os.str();
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string OrbitalSubspaceLJPN::DebugStr() const {
    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index) {
      OrbitalStateLJPN state(*this,state_index);

      os
        << " " << "index"
        << " " << std::setw(width) << state_index
        << " " << "nlj"
        << " " << std::setw(width) << state.n()
        << " " << std::setw(width) << state.l()
        << " " << std::setw(width+2) << state.j().Str()
        << " " << "weight"
        << " " << state.weight()
        << std::endl;
    }

    return os.str();

  }

  /**
   * Flatten subspace into a vector of OrbitalPNInfo objects.
   *
   * @return vector representation of subspace
   */
  OrbitalPNList OrbitalSubspaceLJPN::OrbitalInfo() const
  {
    OrbitalPNList orbitals;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStateLJPN state(*this,state_index);
        orbitals.push_back(state.OrbitalInfo());
      }

    return orbitals;

  }

  /**
   * Generate a string representation of the orbital labels.
   * @return orbital labels as a string
   */
  std::string OrbitalStateLJPN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << std::setw(width) << index()
       << " :"
       << " " << std::setw(width) << n()
       << " " << std::setw(width) << l()
       << " " << std::setw(width) << j()
       << " :"
       << " " << std::setw(width) << weight()
       << " " << "]";

    return os.str();
  }

  /**
   * Flatten state into an OrbitalPNInfo object.
   *
   * @return OrbitalPNInfo representation of state
   */
  OrbitalPNInfo OrbitalStateLJPN::OrbitalInfo() const
  {
    OrbitalPNInfo orbital;

    orbital.orbital_species = orbital_species();
    orbital.n = n();
    orbital.l = l();
    orbital.j = j();
    orbital.weight = weight();

    return orbital;
  }

  /**
   * Construct an Nmax-truncated single-particle space divided into LJPN
   * subspaces.
   *
   * @param[in] Nmax number of oscillator quanta
   */
  OrbitalSpaceLJPN::OrbitalSpaceLJPN(int Nmax) {
    // save truncation
    weight_max_ = double(Nmax);
    is_oscillator_like_ = true;
    Nmax_ = Nmax;

    // iterate over species
    for (OrbitalSpeciesPN orbital_species :
          {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN}) {
      for (int l=0; l<=Nmax; ++l) {
        for (HalfInt j = l-HalfInt(1,2); j<=(l+HalfInt(1,2)); ++j) {
          if (j<0) continue;
          OrbitalSubspaceLJPN subspace(orbital_species,l,j,Nmax);
          PushSubspace(subspace);
        }
      }
    }
  }

  /**
   * Construct a space with LJPN subspaces from a list of orbitals.
   *
   * @param[in] states vector of orbitals
   */
  OrbitalSpaceLJPN::OrbitalSpaceLJPN(const OrbitalPNList& states)
  {
    weight_max_ = 0.0;

    // collect (l,j) subspace labels sorted in canonical order
    std::set<OrbitalSubspaceLJPNLabels> subspace_labels_set;
    for (std::size_t state_index=0; state_index<states.size(); ++state_index)
      {
        OrbitalPNInfo state = states[state_index];
        OrbitalSubspaceLJPNLabels labels(state.orbital_species,state.l,state.j);
        subspace_labels_set.insert(labels);
      }

    // construct subspaces
    for (const OrbitalSubspaceLJPNLabels& labels : subspace_labels_set)
      {
        OrbitalSpeciesPN orbital_species;
        int l;
        HalfInt j;
        std::tie(orbital_species,l,j) = labels;
        OrbitalSubspaceLJPN subspace(orbital_species,l,j,states);
        PushSubspace(subspace);
        weight_max_ = std::max(weight_max_,subspace.weight_max());
      }

    // check if oscillator-like and set Nmax
    if (IsOscillatorLike_()) {
      is_oscillator_like_ = true;
      Nmax_ = int(weight_max_);
    } else {
      is_oscillator_like_ = false;
      Nmax_ = -1;
    }
  }

  /**
   * Check if space is truncated like an Nmax oscillator truncation.
   *
   * @return true if all subspaces are Nmax truncated with the same Nmax.
   */
  bool OrbitalSpaceLJPN::IsOscillatorLike_() const
  {

    // first check that the space contains subspaces
    if (size() < 1)
      return false;

    // check if subspaces are individually Nmax truncated
    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        if (!subspace.is_oscillator_like())
          return false;
      }

    // see if can extract viable Nmax
    int Nmax = GetSubspace(0).Nmax();
    for (std::size_t subspace_index=1; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        if (subspace.Nmax() != Nmax)
          return false;
      }

    // orbitals are oscillator like!
    return true;
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string OrbitalSpaceLJPN::DebugStr() const {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index) {
      const SubspaceType& subspace = GetSubspace(subspace_index);

      os
        << " " << "index"
        << " " << std::setw(width) << subspace_index
        << " " << "species"
        << " " << std::setw(width) << int(subspace.orbital_species())
        << " " << "dim"
        << " " << std::setw(width) << subspace.size()
        << " " << std::endl;
    }

    return os.str();

  }

  /**
   * Flatten space into a vector of OrbitalPNInfo objects.
   *
   * @return vector representation of space
   */
  OrbitalPNList OrbitalSpaceLJPN::OrbitalInfo() const
  {
    OrbitalPNList orbitals;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        OrbitalPNList subspace_orbitals;

        // get orbitals for subspace and append to vector
        subspace_orbitals = subspace.OrbitalInfo();
        orbitals.insert(orbitals.end(),
                        subspace_orbitals.begin(),
                        subspace_orbitals.end()
                       );

      }

    return orbitals;

  }

  OrbitalSectorsLJPN::OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space
    )
    : BaseSectors(bra_space, ket_space), j0_(-1), g0_(-1), Tz0_(1)
  {
    for (std::size_t bra_subspace_index=0; bra_subspace_index<bra_space.size(); ++bra_subspace_index) {
      for (std::size_t ket_subspace_index=0; ket_subspace_index<ket_space.size(); ++ket_subspace_index) {
        // retrieve subspaces
        const SubspaceType& bra_subspace = bra_space.GetSubspace(bra_subspace_index);
        const SubspaceType& ket_subspace = ket_space.GetSubspace(ket_subspace_index);

        // push sector
        PushSector(bra_subspace_index,ket_subspace_index);
      }
    }
  }

  OrbitalSectorsLJPN::OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space,
      int j0, int g0, int Tz0)
    : BaseSectors(bra_space, ket_space), j0_(j0), g0_(g0), Tz0_(Tz0)
  {
    for (std::size_t bra_subspace_index=0; bra_subspace_index<bra_space.size(); ++bra_subspace_index) {
      for (std::size_t ket_subspace_index=0; ket_subspace_index<ket_space.size(); ++ket_subspace_index) {

        // retrieve subspaces
        const SubspaceType& bra_subspace = bra_space.GetSubspace(bra_subspace_index);
        const SubspaceType& ket_subspace = ket_space.GetSubspace(ket_subspace_index);

        bool allowed = true;
        allowed &= am::AllowedTriangle(ket_subspace.j(), j0, bra_subspace.j());
        allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2 == 0);
        allowed &= ((bra_subspace.Tz() - ket_subspace.Tz()) == Tz0);

        // push sector
        if (allowed) {
          PushSector(bra_subspace_index,ket_subspace_index);
        }
      }
    }
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  std::string OrbitalSectorsLJPN::DebugStr() const
  {
    std::ostringstream os;
    int width = 3;

    for (std::size_t sector_index=0; sector_index<size(); ++sector_index) {
      // const basis::OrbitalSectorLJPN& sector = GetSector(sector_index);
      // os << sector_index+1 << sector.DebugStr();
      const SectorType& sector = GetSector(sector_index);
      os << std::setw(width) << sector_index
         << " bra " << std::setw(width) << sector.bra_subspace_index()
         << " (" << int(sector.bra_subspace().orbital_species())
         << ", " << sector.bra_subspace().l()
         << ", " << sector.bra_subspace().j().Str()  << ")"
         << " ket " << std::setw(width) << sector.ket_subspace_index()
         << " (" << int(sector.ket_subspace().orbital_species())
         << ", " << sector.ket_subspace().l()
         << ", " << sector.ket_subspace().j().Str() << ")"
         << std::endl;
    }
    return os.str();
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
