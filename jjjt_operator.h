/****************************************************************
  jjjt_operator.h

  Defines functions for I/O and manipulation of two-body operator
  matrices in jjJT coupling scheme.  Written for use in Moshinsky
  transformation.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/13/16 (mac): Created, adapting code from lsjt_operator.

****************************************************************/

#ifndef JJJT_OPERATOR_H_
#define JJJT_OPERATOR_H_

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "eigen3/Eigen/Core"

#include "basis/lsjt_operator.h"  // for OperatorLabelsJT
#include "basis/jjjt_scheme.h"
#include "basis/operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body JJJT operator -- gather N blocks
  ////////////////////////////////////////////////////////////////

  void GatherOperatorTwoBodyNJJJTToTwoBodyJJJT(
      const basis::OperatorLabelsJT& operator_labels,
      const basis::TwoBodySpaceNJJJT& two_body_njjjt_space,
      const std::array<basis::TwoBodySectorsNJJJT,3>& two_body_njjjt_component_sectors,
      const std::array<basis::MatrixVector,3>& two_body_njjjt_component_matrices,
      const basis::TwoBodySpaceJJJT& two_body_jjjt_space,
      std::array<basis::TwoBodySectorsJJJT,3>& two_body_jjjt_component_sectors,
      std::array<basis::MatrixVector,3>& two_body_jjjt_component_matrices
    );
  // Assemble two-body representation of operator in JJJT basis, from
  // two-body representation in NJJJT basis, i.e., gathering the
  // matrix elements from different N blocks.
  //
  // See notes on "internal representation of an operator in JT
  // scheme" in jjjt_operator.h for the general principles of how the
  // operators are represented.
  //
  // Symmetry: The lower triangle of diagonal sectors is
  // zero-initialized, but not populated.
  //
  // Arguments:
  //   operator_labels (basis::OperatorLabejjjt) : tensorial properties of operator
  //   two_body_njjjt_space (...) : source space
  //   two_body_njjjt_component_sectors (...) : source sectors
  //   two_body_njjjt_component_matrices (...) : source matrices
  //   two_body_jjjt_space (...) : target space
  //   two_body_jjjt_component_sectors (..., output) : target sectors
  //   two_body_jjjt_component_matrices (..., output) : target matrices

  ////////////////////////////////////////////////////////////////
  // two-body LSJT operator output
  ////////////////////////////////////////////////////////////////

  // Note that the primary intention of the output for two-body
  // operators in JJJT scheme is for diagnostic purposes.  If we were
  // to store operators more permanently in this format, we would also
  // want to define an appropriate file header format.  We would also
  // probably replace the human-readable "n/2" fractional j values
  // with integer twice-values.
  //
  // Data lines are of the form:
  //
  //   T0  N1' l1' j1' N2' l2' J' T' g'  N1 l1 j1 N2 l2 j2 J T g  JT-RME
  //
  // Although the g label is redundant (it can be deduced from l1 and
  // l2), it is included to make the sector structure more easily
  // apparent to a human reader.
  //
  // Iteration follows the usual scheme within the basis module:
  // sectors are lexicographic by (bra,ket) subspace indices, then
  // matrix elements within a sector are lexicographic by (bra,ket)
  // state indices.
  //
  // Reminder: One should be sure to document whether one is writing
  // AS or NAS matrix elements!

  void WriteTwoBodyOperatorComponentJJJT(
      std::ostream& os,
      int T0,
      const basis::TwoBodySectorsJJJT& sectors,
      const basis::MatrixVector& matrices,
      basis::NormalizationConversion conversion_mode
    );
  // Write single isospin component of a two-body operator in JJJT
  // scheme.
  //
  // Side effect: The floating point precision attribute of the output
  // stream is modified.
  //
  // Arguments:
  //   os (std::ostream) : text-mode output stream
  //   T0 (int) : isospin for this isospin component
  //   sector (basis::TwoBodySectorsJJJT) : sectors defining operator
  //   matrices (basis::MatrixVector) : matrices defining operator
  //   conversion (basis::NormalizationConversion) : specifies any 
  //     conversion between AS and NAS for output




  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif