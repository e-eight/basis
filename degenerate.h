/**************************************************************
  @file degenerate.h

  Defines basis containers in which each state has associated with it
  a substate degeneracy.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 06/06/17 (mac): Created as multibasis.h, abstracted from code in spncci
    branching_u3s.
  + 06/11/17 (mac): Rename TotalFullDimension to FullDimension.
  + 06/17/17 (mac): Rename to degenerate.h.  Rename multiplicity to degeneracy.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 06/24/21 (pjf):
    - Add dimension() accessor to BaseDegenerateSubspace.
    - Deprecate full_dimension() in BaseDegenerateSubspace and FullDimension()
      in BaseDegenerateSpace.
  + 06/25/21 (pjf):
    - Complete renaming of multiplicity to degeneracy.
    - Add GetStateOffset() and GetStateDegeneracy() to BaseDegenerateSubspace.
    - Repurpose BaseDegenerateSpace to fit its name:
      * Allow multiple degenerate subspaces within BaseDegenerateSpace.
      * Add subspace degeneracy.
      * Override PushSubspace and EmplaceSubspace to account for degeneracies.
  + 06/28/21 (pjf): Pass tSpaceLabelsType through BaseDegenerateSpace
    to BaseSpace.
  + 06/29/21 (pjf): Fix GetSubspaceOffset accessor.
  + 07/02/21 (pjf):
    - Add private convenience typedefs to access base class.
    - Replace `typedef`s with alias declarations (`using`).
    - Add constructors for BaseSubspace and BaseSpace which accept labels.
****************************************************************/

#ifndef BASIS_DEGENERATE_H_
#define BASIS_DEGENERATE_H_

#include <cstddef>
#include <utility>
#include <vector>

#include "basis/basis.h"

namespace basis {


  ////////////////////////////////////////////////////////////////
  // generic subspace
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateSubspace -- holds indexing of states within a symmetry
  // subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  // Template arguments:
  //   tSubspaceLabelsType : tuple for subspace labels, e.g., std::tuple<int,int,int,int,int>
  //   tStateLabelsType : tuple for state labels, e.g., std::tuple<int>
  //
  // Note: Even if only a single integer label is needed, we must use
  // tuple<int> (as opposed to plain int) to make the two forms of the
  // state constructor syntactically distinct.

  template <typename tSubspaceLabelsType, typename tStateLabelsType>
    class BaseDegenerateSubspace
    : public BaseSubspace<tSubspaceLabelsType,tStateLabelsType>
  {

    private:

    ////////////////////////////////////////////////////////////////
    // private (convenience) typedefs
    ////////////////////////////////////////////////////////////////
    using BaseSubspaceType = BaseSubspace<tSubspaceLabelsType,tStateLabelsType>;

    public:

    ////////////////////////////////////////////////////////////////
    //  common typedefs
    ////////////////////////////////////////////////////////////////

    // these typedefs apparently need to be repeated, not inherited
    // from BaseSubspace (they are not recognized as inherited types
    // below)

    using SubspaceLabelsType = tSubspaceLabelsType;
    using StateLabelsType = tStateLabelsType;

    ////////////////////////////////////////////////////////////////
    // general constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    //   Implicitly invoked by derived class.
    BaseDegenerateSubspace() : full_dimension_(0) {}

    // pass-through constructor accepting labels
    explicit BaseDegenerateSubspace(const SubspaceLabelsType& labels)
      : full_dimension_{0}, BaseSubspaceType{labels}
    {}

    ////////////////////////////////////////////////////////////////
    // accessors for substate information
    ////////////////////////////////////////////////////////////////

    const std::vector<std::size_t>& state_offsets() const
    {
      return state_offsets_;
    }

    const std::vector<int>& state_degeneracies() const
    {
      return state_degeneracies_;
    }

    std::size_t GetStateOffset(std::size_t i, int degeneracy_index=1) const
    {
      assert(degeneracy_index <= state_degeneracies_.at(i));
      return state_offsets()[i]+(degeneracy_index-1);
    }

    std::size_t GetStateDegeneracy(std::size_t i) const
    {
      return state_degeneracies_.at(i);
    }

    std::size_t dimension() const
    {
      return full_dimension_;
    }

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use dimension() instead")
    std::size_t full_dimension() const
    {
      return full_dimension_;
    }
#endif

    protected:

    ////////////////////////////////////////////////////////////////
    // state label push (for initial construction)
    ////////////////////////////////////////////////////////////////

    void PushStateLabels(const StateLabelsType& state_labels, int degeneracy)
    // Create indexing information (in both directions, index <->
    // labels) for a state.
    {
      // push state itself
      BaseSubspaceType::PushStateLabels(state_labels);

      // push substate information
      state_offsets_.push_back(full_dimension_);
      state_degeneracies_.push_back(degeneracy);
      full_dimension_ += degeneracy;
    };

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // degeneracy counting information
    std::vector<std::size_t> state_offsets_;  // offset to given state's starting substate
    std::vector<int> state_degeneracies_;  // given state's number of substates
    std::size_t full_dimension_;  // total number of substates

  };

  ////////////////////////////////////////////////////////////////
  // generic state realized within subspace
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateState -- realization of a state withinin a given subspace
  //
  // The derived class is expected to set up a constructor and
  // friendlier accessors for the individual labels.
  //
  // The space (and the indexing it provides) is *not* copied into the
  // state but rather stored by pointer reference.  It should
  // therefore exist for the lifetime of the state object.
  //
  // Template arguments:
  //   tSubspaceType : subspace type in which this state lives

  template <typename tSubspaceType>
    class BaseDegenerateState
    : public BaseState<tSubspaceType>
    {

      private:

      ////////////////////////////////////////////////////////////////
      // private (convenience) typedefs
      ////////////////////////////////////////////////////////////////
      using BaseStateType = BaseState<tSubspaceType>;

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;
      using StateLabelsType = typename SubspaceType::StateLabelsType;

      ////////////////////////////////////////////////////////////////
      // general constructors
      ////////////////////////////////////////////////////////////////

      // default constructor -- disabled
      BaseDegenerateState();

      // copy constructor -- synthesized

      // constructors

      BaseDegenerateState(const SubspaceType& subspace, std::size_t index)
        // Construct state, given index within subspace.
        : BaseStateType(subspace,index)
      {
      }

      BaseDegenerateState(const SubspaceType& subspace, const StateLabelsType& state_labels)
        // Construct state, by reverse lookup on labels within subspace.
        : BaseStateType(subspace,state_labels)
        {
        }

      ////////////////////////////////////////////////////////////////
      // retrieval of substate information
      ////////////////////////////////////////////////////////////////

      std::size_t offset() const
      {
        return BaseStateType::subspace().state_offsets()[BaseStateType::index()];
      }

      int degeneracy() const
      {
        return BaseStateType::subspace().state_degeneracies()[BaseStateType::index()];
      }


    };



  ////////////////////////////////////////////////////////////////
  // generic space
  ////////////////////////////////////////////////////////////////

  // BaseDegenerateSpace -- container to hold subspaces, with reverse lookup by
  // subspace labels
  //
  // Template arguments:
  //   tSubspaceType (typename) : type for subspace

  template <typename tSubspaceType, typename tSpaceLabelsType = void>
    class BaseDegenerateSpace
    : public BaseSpace<tSubspaceType,tSpaceLabelsType>
    {

      private:

      ////////////////////////////////////////////////////////////////
      // private (convenience) typedefs
      ////////////////////////////////////////////////////////////////
      using BaseSpaceType = BaseSpace<tSubspaceType,tSpaceLabelsType>;

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;
      using SpaceLabelsType = tSpaceLabelsType;

      ////////////////////////////////////////////////////////////////
      // general constructors
      ////////////////////////////////////////////////////////////////

      // default constructor
      //   Implicitly invoked by derived class.
      BaseDegenerateSpace() = default;

      // pass-through constructor accepting labels
      //   n.b. this template here is SFINAE black magic
      template<typename T = SpaceLabelsType>
      explicit BaseDegenerateSpace(const T& labels)
        : BaseSpaceType{labels}
      {}

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

#ifdef BASIS_ALLOW_DEPRECATED
      DEPRECATED("use dimension() instead")
      std::size_t FullDimension() const
      // Return the total dimension of all subspaces within the space,
      // taking into account substate degeneracies.
      {
        std::size_t full_dimension = 0;
        for (std::size_t subspace_index=0; subspace_index<BaseSpace<tSubspaceType>::size(); ++subspace_index)
          full_dimension += BaseSpace<tSubspaceType>::GetSubspace(subspace_index).full_dimension();
        return full_dimension;
      }
#endif

      ////////////////////////////////////////////////////////////////
      // subspace lookup and retrieval
      ////////////////////////////////////////////////////////////////

      int GetSubspaceDegeneracy(std::size_t i) const
      /// Given the index for a subspace, return the degeneracy of the
      /// subspace within the space.
      {
        return subspace_degeneracies_.at(i);
      };

      int GetSubspaceOffset(std::size_t i, int degeneracy_index) const
      /// Given the index for a subspace and the degeneracy index, return
      /// the offset of the subspace (with given degeneracy index) within
      /// the space.
      {
        assert(degeneracy_index <= subspace_degeneracies_.at(i));
        return (BaseSpaceType::GetSubspaceOffset(i)
          +(degeneracy_index-1)*this->GetSubspace(i).dimension());
      }

      protected:

      ////////////////////////////////////////////////////////////////
      // subspace push (for initial construction)
      ////////////////////////////////////////////////////////////////

      void PushSubspace(const SubspaceType& subspace) = delete;
      // Prevent use of PushSubspace without degeneracy.

      void PushSubspace(const SubspaceType& subspace, int degeneracy)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        this->subspace_offsets_.push_back(this->dimension_);
        (*(this->lookup_))[subspace.labels()] = this->subspaces_->size();  // index for lookup
        this->subspaces_->push_back(subspace);  // save space
        subspace_degeneracies_.push_back(degeneracy);  // save degeneracy
        this->dimension_ += subspace.dimension()*degeneracy;
      };

      template <class... Args>
      void EmplaceSubspace(Args&&... args) = delete;
      // Prevent use of EmplaceSubspace without degeneracy.

      template <class... Args>
      void EmplaceSubspace(Args&&... args, int degeneracy)
      /// Create indexing information (in both directions, index <->
      /// labels) for a subspace.
      {
        std::size_t index = this->subspaces_->size();  // index for lookup
        this->subspace_offsets_.push_back(this->dimension_);
        this->subspaces_->emplace_back(std::forward<Args>(args)...);
        (*(this->lookup_))[this->subspaces_->back().labels()] = index;
        subspace_degeneracies_.push_back(degeneracy);  // save degeneracy
        this->dimension_ += this->subspaces_->back().dimension()*degeneracy;
      }

    ////////////////////////////////////////////////////////////////
    // private storage
    ////////////////////////////////////////////////////////////////

    // degeneracy counting information
    std::vector<int> subspace_degeneracies_;  // given subspace's number of sub-subspaces

    };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // BASIS_DEGENERATE_H_
