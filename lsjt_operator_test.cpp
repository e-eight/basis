/****************************************************************
  lsjt_operator_test.cpp

  Mark A. Caprio, University of Notre Dame.
  See interaction_lsjt.h for documentation.

****************************************************************/

#include <fstream>
#include <iomanip>

#include "lsjt_operator.h"

////////////////////////////////////////////////////////////////
// test parts
////////////////////////////////////////////////////////////////

void WriteTestRelative(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  std::cout << "Setup" << std::endl;

  // set operator and file header parameters
  basis::RelativeOperatorParametersLSJT parameters;
  // file format
  parameters.version=1;
  // operator tensorial parameters
  parameters.J0=0;
  parameters.g0=0;
  parameters.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
  parameters.T0_min=0;
  parameters.T0_max=2;
  // relative basis parameters
  parameters.Nmax=2;
  parameters.Jmax=parameters.Nmax+1;

  // set up relative space
  basis::RelativeSpaceLSJT space(parameters.Nmax,parameters.Jmax);

  // set up operator containers
  //
  // These are arrays to store information for T0=0/1/2 components.
  std::array<basis::RelativeSectorsLSJT,3> component_sectors;
  std::array<basis::MatrixVector,3> component_matrices;

  // populate operator containers
  //
  // Note: This can now instead be done automatically using
  // ConstructIdentityOperatorRelativeLSJT.

  for (int T0=parameters.T0_min; T0<=parameters.T0_max; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::RelativeSectorsLSJT(space,parameters.J0,T0,parameters.g0);
      std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;
          
      // populate matrices
      if (T0==0)
        basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
      else
        basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
    }

  ////////////////////////////////////////////////////////////////
  // write test
  ////////////////////////////////////////////////////////////////

  // set up stream for output
  std::ostringstream os;

  // write header parameters
  basis::WriteRelativeOperatorParametersLSJT(os,parameters);

  // write matrices
  for (int T0=parameters.T0_min; T0<=parameters.T0_max; ++T0)
    {
      basis::WriteRelativeOperatorComponentLSJT(
          os,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
    }

  // dump to terminal for inspection
  std::cout << os.str();

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();
}

void ReadTestRelative(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // read test
  ////////////////////////////////////////////////////////////////

  std::cout << "Readback test" << std::endl;

  // set up stream for readback
  std::ifstream is(filename.c_str());

  // read header parameters
  basis::RelativeOperatorParametersLSJT parameters;
  basis::ReadRelativeOperatorParametersLSJT(is,parameters);
  int J0=parameters.J0;
  int g0=parameters.g0;
  // and inspect...
  basis::WriteRelativeOperatorParametersLSJT(std::cout,parameters);

  // set up space
  basis::RelativeSpaceLSJT space(parameters.Nmax,parameters.Jmax);

  // read matrices
  std::array<basis::RelativeSectorsLSJT,3> component_sectors;
  std::array<basis::MatrixVector,3> component_matrices;

  for (int T0=parameters.T0_min; T0<=parameters.T0_max; ++T0)
    // for each isospin component
    {
      // enumerate sectors
      component_sectors[T0] = basis::RelativeSectorsLSJT(space,J0,T0,g0);

      // read matrices
      basis::ReadRelativeOperatorComponentLSJT(
          is,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
      // and inspect...
      basis::WriteRelativeOperatorComponentLSJT(
          std::cout,
          T0,
          component_sectors[T0],component_matrices[T0]
        );
    }
}

void WriteTestTwoBody(const std::string& filename)
{
  ////////////////////////////////////////////////////////////////
  // defining operator
  ////////////////////////////////////////////////////////////////

  // since we set the operator to a "naive identity" operator, the
  // matrix elements are taken to be NAS

  std::cout << "Setup" << std::endl;

  // set up space
  int Nmax = 2;
  basis::TwoBodySpaceLSJT space(Nmax);

  // set up operator containers
  //
  // These are vectors to store information for T0=0/1/2 components.
  std::array<basis::TwoBodySectorsLSJT,3> component_sectors;
  std::array<basis::MatrixVector,3> component_matrices;

  // populate operator containers
  int J0 = 0;
  int g0 = 0;
  for (int T0=0; T0<=2; ++T0)
    // for each isospin component
    {

      // enumerate sectors
      component_sectors[T0] = basis::TwoBodySectorsLSJT(space,J0,T0,g0);
      std::cout << " T0 " << T0 << " size " << component_sectors[T0].size() << std::endl;
          
      // populate matrices
      if (T0==0)
        basis::SetOperatorToIdentity(component_sectors[T0],component_matrices[T0]);
      else
        basis::SetOperatorToZero(component_sectors[T0],component_matrices[T0]);
    }

  ////////////////////////////////////////////////////////////////
  // write test
  ////////////////////////////////////////////////////////////////

  // set up stream for output
  std::ostringstream os;

  // write matrices
  int T0_min=0;
  int T0_max=2;
  for (int T0=T0_min; T0<=T0_max; ++T0)
    {
      basis::WriteTwoBodyOperatorComponentLSJT(
          os,
          T0,
          component_sectors[T0],component_matrices[T0],
          basis::NormalizationConversion::kNone
        );
    }

  // dump to file
  std::ofstream ofile(filename.c_str());
  ofile << os.str();
}

void IdentityTest()
{
  ////////////////////////////////////////////////////////////////
  // construct relative identity operator
  ////////////////////////////////////////////////////////////////

  // define operator properties
  int Nmax_relative = 4;
  int Jmax_relative = Nmax_relative+1;

  basis::OperatorLabelsJT operator_labels;
  operator_labels.J0 = 0;
  operator_labels.g0 = 0;
  operator_labels.T0_min = 0;
  operator_labels.T0_max = 2;
  operator_labels.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

  // define space and operator containers
  basis::RelativeSpaceLSJT relative_space(Nmax_relative,Jmax_relative);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;

  // do construction
  ConstructIdentityOperatorRelativeLSJT(
      operator_labels,
      relative_space,relative_component_sectors,relative_component_matrices
    );

  // and now it would be nice to inspect the contents, wouldn't it,
  // ...

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  std::string relative_filename("test/lsjt_operator_test_relative_identity_Nmax02.dat");
  WriteTestRelative(relative_filename);
  ReadTestRelative(relative_filename);
 
  std::string two_body_filename("test/lsjt_operator_test_two_body_identity_nas_Nmax02.dat");
  WriteTestTwoBody(two_body_filename);

  IdentityTest();

  // termination
  return 0;
}
