/************************************************************//**
  @file many_body.h

  Provides convenience definitions for many-body indexing schemes.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 07/19/16 (mac): Created, incorporating code from operator.h.
  + 09/06/19 (pjf): Moved WeightMax from jjjpn_scheme.h.

****************************************************************/

#ifndef BASIS_MANY_BODY_H_
#define BASIS_MANY_BODY_H_

#include <cassert>
#include <tuple>

namespace basis {

  ////////////////////////////////////////////////////////////////
  // many-body truncation rank
  ////////////////////////////////////////////////////////////////

  enum class Rank {kOneBody=1, kTwoBody=2, kThreeBody=3, kFourBody=4};
  // Truncation rank for many-body bases.
  //
  // To be used as a mode parameter for basis constructors, to
  // determine the many-body rank of the truncation scheme, e.g.,
  // one-body ("square") truncation, two-body ("triangular")
  // truncation.

  std::tuple<int,int> TwoBodyCutoffs(basis::Rank truncation_rank, int truncation_cutoff);
  // Calculate one-body N1max and two-body N2max cutoffs for truncation scheme.
  //
  // Example:
  //    int N1max, N2max;
  //    std::tie(N1max,N2max) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);
  //
  // Arguments:
  //   truncation_rank (basis::Rank) : many-body rank (kOneBody or kTwoBody)
  //   truncation_cutoff (int) : oscillator cutoff for that rank
  //
  // Returns:
  //   (std::tuple<int,int>) : one-body and two-body cutoffs (N1max,N2max)

  // maximum weight collection
  struct WeightMax
  {

    // constructor

    WeightMax() {};
    // default constructor

    WeightMax(double wp, double wn, double wpp, double wnn, double wpn)
    // trivial constructor
    {
      one_body[0] = wp;
      one_body[1] = wn;
      two_body[0] = wpp;
      two_body[1] = wnn;
      two_body[2] = wpn;
    }

    WeightMax(int N1max, int N2max)
    // Set conventional oscillator one-body/two-body truncation from
    // separate N1max and N2max.
    {
      one_body[0] = N1max;
      one_body[1] = N1max;
      two_body[0] = N2max;
      two_body[1] = N2max;
      two_body[2] = N2max;
    }

    WeightMax(basis::Rank truncation_rank, int truncation_cutoff)
    // Set conventional oscillator one-body/two-body truncation from
    // either a given one-body truncation or a given two-body
    // truncation.
    {
      // extract one-body and two-body cutoffs
      int N1max, N2max;
      std::tie(N1max,N2max) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

      // save cutoffs
      one_body[0] = N1max;
      one_body[1] = N1max;
      two_body[0] = N2max;
      two_body[1] = N2max;
      two_body[2] = N2max;
    }

    // maximum weights
    std::array<double,2> one_body;
    std::array<double,3> two_body;

    // truncation information -- TODO (mac)? but convert to class
    //
    // bool is_oscillator_like() const {return is_oscillator_like_;}
    // int N1max() const {assert(is_oscillator_like()); return N1max_;}
    // int N2max() const {assert(is_oscillator_like()); return N2max_;}
    // // only meaningful if oscillator scheme constructor used
    // bool is_oscillator_like_;
    // int N1max_, N2max_;  // only meaningful if oscillator scheme constructor used

  };

  ////////////////////////////////////////////////////////////////
  // many-body normalization conversion flag
  ////////////////////////////////////////////////////////////////

  enum class NormalizationConversion {kNone, kASToNAS, kNASToAS};
  // Normalization scheme converson flag.
  //
  // To be used as value for a mode parameter to requests on-the-fly
  // conversion between antisymmetrized (AS) and normalized
  // antisymmetrized (NAS) matrix elements on input/output.



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
