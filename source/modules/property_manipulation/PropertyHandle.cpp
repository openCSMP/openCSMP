#include "PropertyHandle.h"

#include "FEM_Data.h"
#include "TensorVariable.h"
#include "VectorVariable.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "ModelSubDomain.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "TextInterface.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace csmp {

// ============================================================================
//  Shared constructor body
// ============================================================================

template<uint32_t dim>
void PropertyHandle<dim>::Initialise( const char*   subdomain_name,
                                        const char*   var_name,
                                        VARIABLE_TYPE type,
                                        PLACEMENT     place,
                                        uint32_t      vsize )
{
    subdomain_name_    = subdomain_name;
    region_            = nullptr;
    boundary_          = nullptr;
    split_boundary_    = nullptr;
    is_split_boundary_ = false;

    if ( model_.ContainsRegion( subdomain_name ) )
    {
        region_ = &model_.Region( subdomain_name );
    }
    else if ( model_.ContainsBoundary( subdomain_name ) )
    {
        boundary_ = &model_.Boundary( subdomain_name );
    }
    else if ( model_.ContainsSplitBoundary( subdomain_name ) )
    {
        split_boundary_    = &model_.SplitBoundary( subdomain_name );
        is_split_boundary_ = true;
    }

    if ( place == BOUNDARY )
        throw Exception( ERROR,
            "PropertyHandle::Initialise", var_name,
            "Use FACE or INTER_FACE placement rather than BOUNDARY." );

    // ------------------------------------------------------------------
    //  Finite-volume placements require the FV stencils to be
    //  initialised before any FV property can be created. The sector
    //  and facet counts on each element are set during
    //  InstantiateFiniteVolumes() — if this has not been called,
    //  Element::Sectors() and Element::Facets() return 0, which causes
    //  AddProperty to crash with an out-of-bounds access in
    //  LocalVariableStorage.
    //
    //  We call InstantiateFiniteVolumes() automatically here so that
    //  the user does not need to remember to call it before constructing
    //  a PropertyHandle for a FV placement.
    // ------------------------------------------------------------------
    const bool isFVPlacement =
        ( place == SECTOR_INTEGRATION_POINT              ||
          place == FACET_INTEGRATION_POINT               ||
          place == FACE_SECTOR_INTEGRATION_POINT         ||
          place == FACE_FACET_INTEGRATION_POINT          ||
          place == INTER_FACE_SECTOR_INTEGRATION_POINT   ||
          place == INTER_FACE_FACET_INTEGRATION_POINT    );

    if ( isFVPlacement )
    {
        bool fv_initialised = false;
        WithSubdomain( [&]( auto& sd )
        {
            if ( !sd.CellVector().empty() )
                fv_initialised = ( sd.CellVector().front()->Sectors() > 0 );
        });

        if ( !fv_initialised )
        {
            cerr << "\nPropertyHandle::Initialise: '"
                 << var_name
                 << "' requires finite-volume stencils. "
                    "Calling model_.InstantiateFiniteVolumes() automatically.\n";
            model_.InstantiateFiniteVolumes();
        }
    }

    var_name_ = var_name;

    if ( model_.Database().IsDefined( var_name ) )
    {
        const csmp::Index existing = model_.Database().StorageKey( var_name );

        if ( type != SCALAR && existing.type != type )
            cerr << "\nPropertyHandle: type of existing variable '"
                 << var_name << "' differs from requested type.\n";

        if ( place != NODE && existing.place != place )
            cerr << "\nPropertyHandle: placement of existing variable '"
                 << var_name << "' differs from requested placement.\n";

        owns_variable_ = false;
    }
    else
    {
        model_.CreateProperty( var_name, var_name, "SI", type, place, vsize );
        owns_variable_ = true;
    }

    key_ = model_.Database().StorageKey( var_name );
}




// ============================================================================
//  Constructors / destructor
// ============================================================================

template<uint32_t dim>
PropertyHandle<dim>::PropertyHandle( Model<dim>&   model,
                                        const char*   var_name,
                                        VARIABLE_TYPE type,
                                        PLACEMENT     place,
                                        uint32_t      vsize )
    : model_( model ),
      region_( nullptr ),
      boundary_( nullptr ),
      split_boundary_( nullptr ),
      is_split_boundary_( false ),
      flag_output_( ANY ),
      owns_variable_( false )
{
    Initialise( "Model", var_name, type, place, vsize );
}



template<uint32_t dim>
PropertyHandle<dim>::PropertyHandle( Model<dim>&   model,
                                        const char*   subdomain_name,
                                        const char*   var_name,
                                        VARIABLE_TYPE type,
                                        PLACEMENT     place,
                                        uint32_t      vsize )
    : model_( model ),
      region_( nullptr ),
      boundary_( nullptr ),
      split_boundary_( nullptr ),
      is_split_boundary_( false ),
      flag_output_( ANY ),
      owns_variable_( false )
{
    Initialise( subdomain_name, var_name, type, place, vsize );
}




template<uint32_t dim>
PropertyHandle<dim>::PropertyHandle( const PropertyHandle<dim>& other )
    : model_( other.model_ ),
      region_( other.region_ ),
      boundary_( other.boundary_ ),
      split_boundary_( other.split_boundary_ ),
      subdomain_name_( other.subdomain_name_ ),
      is_split_boundary_( other.is_split_boundary_ ),
      var_name_( other.var_name_ ),
      key_( other.key_ ),
      flag_output_( other.flag_output_ ),
      owns_variable_( false )
{
}



template<uint32_t dim>
PropertyHandle<dim>::~PropertyHandle()
{
    if ( owns_variable_ &&
         model_.Database().IsDefined( var_name_.c_str() ) )
        model_.DeleteProperty( var_name_.c_str() );
}

// ============================================================================
//  Queries
// ============================================================================

template<uint32_t dim>
const char* PropertyHandle<dim>::VariableName() const
{
    return var_name_.c_str();
}

template<uint32_t dim>
const csmp::Index& PropertyHandle<dim>::Key() const
{
    return key_;
}

template<uint32_t dim>
VARIABLE_FLAG PropertyHandle<dim>::OutputCondition() const
{
    return flag_output_;
}

template<uint32_t dim>
void PropertyHandle<dim>::OutputCondition( VARIABLE_FLAG c )
{
    flag_output_ = c;
}



template<uint32_t dim>
void PropertyHandle<dim>::Range( double& omin, double& omax ) const
{
    model_.MinMaxOf( var_name_.c_str(), omin, omax );
}




template<uint32_t dim>
bool PropertyHandle<dim>::IsWithinRange() const
{
    double pmin, pmax;
    model_.Database().RangeOf( var_name_.c_str(), pmin, pmax );

    // ------------------------------------------------------------------
    //  NODE placement on SplitBoundary — handled separately for all
    //  types because SplitBoundary has no NodeVector().
    // ------------------------------------------------------------------
    if ( key_.place == NODE && is_split_boundary_ )
    {
        bool in_range = true;
        bool warned   = false;

        std::unordered_set<Node<dim>*> visited;
        for ( auto& ifit : split_boundary_->CellVector() )
            for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
            {
                auto checkNode = [&]( Node<dim>* nd ) -> bool
                {
                    if ( !visited.insert( nd ).second ) return true;

                    switch ( key_.type )
                    {
                        case SCALAR:
                        {
                            double sc = nd->Read( key_ );
                            return ( sc >= pmin && sc <= pmax );
                        }
                        case VECTOR:
                        {
                            VectorVariable<dim> vc;
                            nd->Read( key_, vc );
                            const double len = vc.Length();
                            return ( len >= pmin && len <= pmax );
                        }
                        case TENSOR:
                        {
                            TensorVariable<dim> ts;
                            nd->Read( key_, ts );
                            if constexpr ( dim == 1 )
                                return ( ts(0,0) >= pmin && ts(0,0) <= pmax );
                            else if constexpr ( dim == 2 )
                            {
                                for ( uint32_t ii = 0; ii < 2; ++ii )
                                    for ( uint32_t j = 0; j < 2; ++j )
                                        if ( ts(ii,j) < pmin || ts(ii,j) > pmax )
                                            return false;
                                return true;
                            }
                            else
                            {
                                double e0, e1, e2;
                                if ( !ts.EigenValuesPositiveDefiniteSymmetricMatrix(
                                         e0, e1, e2 ) )
                                {
                                    if ( !warned )
                                    {
                                        ErrorHandler::Instance().Note( WARNING,
                                            "PropertyHandle3::IsWithinRange",
                                            ( string("Tensor variable '")
                                              + var_name_
                                              + "' is not symmetric positive "
                                              + "definite at one or more "
                                              + "SplitBoundary nodes." ).c_str() );
                                        warned = true;
                                    }
                                    return true;
                                }
                                if ( e0 < pmin || e0 > pmax ) return false;
                                if ( e2 < pmin || e2 > pmax ) return false;
                                return true;
                            }
                        }
                        case ARRAY:
                        {
                            ArrayVariable av;
                            nd->Read( key_, av );
                            return av.IsWithinRange( pmin, pmax );
                        }
                        case FLAGGEDARRAY:
                        {
                            FlaggedArrayVariable fav;
                            nd->Read( key_, fav );
                            return fav.IsWithinRange( pmin, pmax );
                        }
                        default:
                            return true;
                    }
                };

                for ( auto side : { INSIDE, OUTSIDE } )
                {
                    if ( !checkNode( ifit->N( i, side ) ) )
                    {
                        in_range = false;
                        goto done_splitboundary_node;
                    }
                }
                if ( ifit->HasInterveningElement() )
                {
                    if ( !checkNode( ifit->N( i, MIDDLE ) ) )
                    {
                        in_range = false;
                        goto done_splitboundary_node;
                    }
                }
            }
        done_splitboundary_node:
        return in_range;
    }

    // ------------------------------------------------------------------
    //  All other placements and subdomains.
    // ------------------------------------------------------------------
    switch ( key_.type )
    {
        case SCALAR:
        {
            double omin, omax;
            WithSubdomain( [&]( auto& sd )
            {
                sd.MinMaxOf( var_name_.c_str(), omin, omax );
            });
            return ( omin >= pmin && omax <= pmax );
        }

        case VECTOR:
        {
            double omin, omax;
            WithSubdomain( [&]( auto& sd )
            {
                sd.MinMaxOf( var_name_.c_str(), omin, omax );
            });
            return ( omin >= pmin && omax <= pmax );
        }

        case TENSOR:
        {
            if constexpr ( dim == 1 )
            {
                double omin, omax;
                WithSubdomain( [&]( auto& sd )
                {
                    sd.MinMaxOf( var_name_.c_str(), omin, omax );
                });
                return ( omin >= pmin && omax <= pmax );
            }
            else if constexpr ( dim == 2 )
            {
                double omin, omax;
                WithSubdomain( [&]( auto& sd )
                {
                    sd.MinMaxOf( var_name_.c_str(), omin, omax );
                });
                return ( omin >= pmin && omax <= pmax );
            }
            else
            {
                TensorVariable<dim> ts;
                bool warned   = false;
                bool in_range = true;

                auto checkTensor = [&]( TensorVariable<dim>& t ) -> bool
                {
                    double e0, e1, e2;
                    if ( !t.EigenValuesPositiveDefiniteSymmetricMatrix(
                             e0, e1, e2 ) )
                    {
                        if ( !warned )
                        {
                            ErrorHandler::Instance().Note( WARNING,
                                "PropertyHandle3::IsWithinRange",
                                ( string("Tensor variable '")
                                  + var_name_
                                  + "' is not symmetric positive definite "
                                  + "at one or more points; eigenvalue "
                                  + "range check skipped for those points."
                                ).c_str() );
                            warned = true;
                        }
                        return true;
                    }
                    if ( e0 < pmin || e0 > pmax ) return false;
                    if ( e2 < pmin || e2 > pmax ) return false;
                    return true;
                };

                switch ( key_.place )
                {
                    case NODE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& nit : sd.NodeVector() )
                            {
                                nit->Read( key_, ts );
                                if ( !checkTensor( ts ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                        });
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                            {
                                cit->Read( key_, ts );
                                if ( !checkTensor( ts ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0;
                                      ip < cit->IntegrationPoints(); ++ip )
                                {
                                    cit->Read( ip, key_, ts );
                                    if ( !checkTensor( ts ) )
                                    {
                                        in_range = false;
                                        return;
                                    }
                                }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0;
                                      s < cit->Sectors(); ++s )
                                {
                                    cit->Read( s, 0U, key_, ts );
                                    if ( !checkTensor( ts ) )
                                    {
                                        in_range = false;
                                        return;
                                    }
                                }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0;
                                      fac < cit->Facets(); ++fac )
                                {
                                    cit->Read( fac, 0U, key_, ts );
                                    if ( !checkTensor( ts ) )
                                    {
                                        in_range = false;
                                        return;
                                    }
                                }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            sd.Read( key_, ts );
                            if ( !checkTensor( ts ) )
                                in_range = false;
                        });
                        break;

                    default:
                        return true;
                }
                return in_range;
            }
        }

        // ------------------------------------------------------------------
        //  ARRAY: use ArrayVariable::IsWithinRange(pmin,pmax) which checks
        //  the min and max element values against the registered range.
        // ------------------------------------------------------------------
        case ARRAY:
        {
            bool in_range = true;
            ArrayVariable av;

            switch ( key_.place )
            {
                case NODE:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& nit : sd.NodeVector() )
                        {
                            nit->Read( key_, av );
                            if ( !av.IsWithinRange( pmin, pmax ) )
                            {
                                in_range = false;
                                return;
                            }
                        }
                    });
                    break;

                case ELEMENT:
                case FACE:
                case INTER_FACE:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                        {
                            cit->Read( key_, av );
                            if ( !av.IsWithinRange( pmin, pmax ) )
                            {
                                in_range = false;
                                return;
                            }
                        }
                    });
                    break;

                case ELEMENT_INTEGRATION_POINT:
                case FACE_INTEGRATION_POINT:
                case INTER_FACE_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t ip = 0;
                                  ip < cit->IntegrationPoints(); ++ip )
                            {
                                cit->Read( ip, key_, av );
                                if ( !av.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case SECTOR_INTEGRATION_POINT:
                case FACE_SECTOR_INTEGRATION_POINT:
                case INTER_FACE_SECTOR_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                            {
                                cit->Read( s, 0U, key_, av );
                                if ( !av.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case FACET_INTEGRATION_POINT:
                case FACE_FACET_INTEGRATION_POINT:
                case INTER_FACE_FACET_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                            {
                                cit->Read( fac, 0U, key_, av );
                                if ( !av.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case REGION:
                    WithSubdomain( [&]( auto& sd )
                    {
                        sd.Read( key_, av );
                        if ( !av.IsWithinRange( pmin, pmax ) )
                            in_range = false;
                    });
                    break;

                default:
                    return true;
            }
            return in_range;
        }

        // ------------------------------------------------------------------
        //  FLAGGEDARRAY: use FlaggedArrayVariable::IsWithinRange(pmin,pmax)
        //  which checks all elements regardless of their flags.
        // ------------------------------------------------------------------
        case FLAGGEDARRAY:
        {
            bool in_range = true;
            FlaggedArrayVariable fav;

            switch ( key_.place )
            {
                case NODE:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& nit : sd.NodeVector() )
                        {
                            nit->Read( key_, fav );
                            if ( !fav.IsWithinRange( pmin, pmax ) )
                            {
                                in_range = false;
                                return;
                            }
                        }
                    });
                    break;

                case ELEMENT:
                case FACE:
                case INTER_FACE:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                        {
                            cit->Read( key_, fav );
                            if ( !fav.IsWithinRange( pmin, pmax ) )
                            {
                                in_range = false;
                                return;
                            }
                        }
                    });
                    break;

                case ELEMENT_INTEGRATION_POINT:
                case FACE_INTEGRATION_POINT:
                case INTER_FACE_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t ip = 0;
                                  ip < cit->IntegrationPoints(); ++ip )
                            {
                                cit->Read( ip, key_, fav );
                                if ( !fav.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case SECTOR_INTEGRATION_POINT:
                case FACE_SECTOR_INTEGRATION_POINT:
                case INTER_FACE_SECTOR_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                            {
                                cit->Read( s, 0U, key_, fav );
                                if ( !fav.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case FACET_INTEGRATION_POINT:
                case FACE_FACET_INTEGRATION_POINT:
                case INTER_FACE_FACET_INTEGRATION_POINT:
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& cit : sd.CellVector() )
                            for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                            {
                                cit->Read( fac, 0U, key_, fav );
                                if ( !fav.IsWithinRange( pmin, pmax ) )
                                {
                                    in_range = false;
                                    return;
                                }
                            }
                    });
                    break;

                case REGION:
                    WithSubdomain( [&]( auto& sd )
                    {
                        sd.Read( key_, fav );
                        if ( !fav.IsWithinRange( pmin, pmax ) )
                            in_range = false;
                    });
                    break;

                default:
                    return true;
            }
            return in_range;
        }

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::IsWithinRange",
                var_name_.c_str(),
                "Unrecognised variable type." );
    }
}





// ============================================================================
//  Interpolate
// ============================================================================

template<uint32_t dim>
void PropertyHandle<dim>::Interpolate( const char* src_name,
                                         const char* dst_name,
                                         PLACEMENT   src,
                                         PLACEMENT   dst )
{
    // ------------------------------------------------------------------
    //  SplitBoundary: only non-node interpolation paths are supported.
    //  InterpolateNodeToCellProperty and
    //  InterpolateNodeToIntegrationPointProperty are deleted on
    //  SplitBoundary, and extrapolation to nodes is not available
    //  because SplitBoundary has no NodeVector().
    // ------------------------------------------------------------------
    if ( is_split_boundary_ )
    {
        if ( src == NODE || dst == NODE )
            throw csmp::Exception( ERROR,
                "PropertyHandle::Interpolate",
                src_name,
                ( std::string(
                    "Interpolation involving NODE placement is not "
                    "supported for SplitBoundary. "
                    "InterpolateNodeToCellProperty and "
                    "InterpolateNodeToIntegrationPointProperty are "
                    "deleted on SplitBoundary, and extrapolation to "
                    "nodes is not available because SplitBoundary has "
                    "no NodeVector(). "
                    "Each InterFace has two distinct node sets: INSIDE "
                    "nodes belonging to Parent(INSIDE) and OUTSIDE nodes "
                    "belonging to Parent(OUTSIDE). "
                    "Use SplitBoundary::InsideNodes() or "
                    "SplitBoundary::OutsideNodes() to obtain the node "
                    "vector for the appropriate side, then use "
                    "InterFace::Parent(INSIDE) or "
                    "InterFace::Parent(OUTSIDE) with the finite element "
                    "interpolation machinery directly." )
                ).c_str() );

        if      ( dst == INTER_FACE &&
                  src == INTER_FACE_INTEGRATION_POINT )
            split_boundary_->InterpolateIntegrationPointToCellProperty(
                src_name, dst_name );
        else if ( dst == INTER_FACE_INTEGRATION_POINT &&
                  src == INTER_FACE )
            split_boundary_->ExtrapolateCellToIntegrationPointProperty(
                src_name, dst_name );
        else if ( dst == INTER_FACE_FACET_INTEGRATION_POINT &&
                  src == INTER_FACE )
            split_boundary_->ExtrapolateCellToFacetIntegrationPointProperty(
                src_name, dst_name );
        else
            throw csmp::Exception( ERROR,
                "PropertyHandle::Interpolate",
                src_name,
                ( std::string(
                    "No interpolation path exists between the two "
                    "placements for SplitBoundary. "
                    "Supported paths are: "
                    "INTER_FACE <-> INTER_FACE_INTEGRATION_POINT and "
                    "INTER_FACE -> INTER_FACE_FACET_INTEGRATION_POINT." )
                ).c_str() );
        return;
    }

    // ------------------------------------------------------------------
    //  Region and Boundary: dispatch directly to region_ or boundary_
    //  to avoid instantiating deleted SplitBoundary methods through
    //  the WithSubdomain generic lambda.
    // ------------------------------------------------------------------
    auto interpolateOn = [&]( auto& sd )
    {
        if      ( dst == ELEMENT && src == NODE )
            sd.InterpolateNodeToCellProperty( src_name, dst_name );
        else if ( dst == NODE && src == ELEMENT )
            sd.ExtrapolateCellToNodeProperty( src_name, dst_name );
        else if ( dst == ELEMENT_INTEGRATION_POINT && src == NODE )
            sd.InterpolateNodeToIntegrationPointProperty( src_name, dst_name );
        else if ( dst == NODE && src == ELEMENT_INTEGRATION_POINT )
            sd.ExtrapolateIntegrationPointToNodeProperty( src_name, dst_name );
        else if ( dst == ELEMENT && src == ELEMENT_INTEGRATION_POINT )
            sd.InterpolateIntegrationPointToCellProperty( src_name, dst_name );
        else if ( dst == ELEMENT_INTEGRATION_POINT && src == ELEMENT )
            sd.ExtrapolateCellToIntegrationPointProperty( src_name, dst_name );
        else if ( dst == FACET_INTEGRATION_POINT && src == ELEMENT )
            sd.ExtrapolateCellToFacetIntegrationPointProperty( src_name, dst_name );
        // FACE placements
        else if ( dst == FACE && src == NODE )
            sd.InterpolateNodeToCellProperty( src_name, dst_name );
        else if ( dst == NODE && src == FACE )
            sd.ExtrapolateCellToNodeProperty( src_name, dst_name );
        else if ( dst == FACE_INTEGRATION_POINT && src == NODE )
            sd.InterpolateNodeToIntegrationPointProperty( src_name, dst_name );
        else if ( dst == NODE && src == FACE_INTEGRATION_POINT )
            sd.ExtrapolateIntegrationPointToNodeProperty( src_name, dst_name );
        else if ( dst == FACE && src == FACE_INTEGRATION_POINT )
            sd.InterpolateIntegrationPointToCellProperty( src_name, dst_name );
        else if ( dst == FACE_INTEGRATION_POINT && src == FACE )
            sd.ExtrapolateCellToIntegrationPointProperty( src_name, dst_name );
        else if ( dst == FACE_FACET_INTEGRATION_POINT && src == FACE )
            sd.ExtrapolateCellToFacetIntegrationPointProperty( src_name, dst_name );
        else
            throw csmp::Exception( ERROR,
                "PropertyHandle::Interpolate",
                src_name,
                "No interpolation path exists between the two placements." );
    };

    if      ( region_   ) interpolateOn( *region_ );
    else if ( boundary_ ) interpolateOn( *boundary_ );
}





// ============================================================================
//  Cross-type assignment helpers
// ============================================================================

template<uint32_t dim>
void PropertyHandle<dim>::AssignVectorLengthToScalar(
    const csmp::Index& rkey )
{
    VectorVariable<dim> vc;

    auto store = [&]( auto* pt )
    {
        pt->Read( rkey, vc );
        pt->Store( key_, makeScalar( flag_output_, vc.Length() ) );
    };

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( key_ ) == flag_output_ )
                            store( nit );
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                                store( nd );
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                                store( nd );
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( key_ ) == flag_output_ )
                        store( cit );
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, key_ ) == flag_output_ )
                        {
                            cit->Read( ip, rkey, vc );
                            cit->Store( ip, key_,
                                makeScalar( flag_output_, vc.Length() ) );
                        }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( s, 0U, rkey, vc );
                            cit->Store( s, 0U, key_,
                                makeScalar( flag_output_, vc.Length() ) );
                        }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( fac, 0U, rkey, vc );
                            cit->Store( fac, 0U, key_,
                                makeScalar( flag_output_, vc.Length() ) );
                        }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                if ( sd.Status( key_ ) == flag_output_ )
                    store( &sd );
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::AssignVectorLengthToScalar",
                var_name_.c_str(), "Unsupported placement." );
    }
}




template<uint32_t dim>
void PropertyHandle<dim>::AssignTensorDetToScalar(
    const csmp::Index& rkey )
{
    TensorVariable<dim> ts;

    auto store = [&]( auto* pt )
    {
        pt->Read( rkey, ts );
        pt->Store( key_, makeScalar( flag_output_, ts.Determinant() ) );
    };

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( key_ ) == flag_output_ )
                            store( nit );
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                                store( nd );
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                                store( nd );
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( key_ ) == flag_output_ )
                        store( cit );
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, key_ ) == flag_output_ )
                        {
                            cit->Read( ip, rkey, ts );
                            cit->Store( ip, key_,
                                makeScalar( flag_output_, ts.Determinant() ) );
                        }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( s, 0U, rkey, ts );
                            cit->Store( s, 0U, key_,
                                makeScalar( flag_output_, ts.Determinant() ) );
                        }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( fac, 0U, rkey, ts );
                            cit->Store( fac, 0U, key_,
                                makeScalar( flag_output_, ts.Determinant() ) );
                        }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                if ( sd.Status( key_ ) == flag_output_ )
                    store( &sd );
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::AssignTensorDetToScalar",
                var_name_.c_str(), "Unsupported placement." );
    }
}



template<uint32_t dim>
void PropertyHandle<dim>::AssignFlaggedArrayToArray( const csmp::Index& rkey )
{
    ArrayVariable        av;
    FlaggedArrayVariable fav;

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( key_ ) == flag_output_ )
                        {
                            nit->Read( rkey, fav );
                            nit->Read( key_, av );
                            av.CopyValuesOnly( fav );
                            nit->Store( key_, av );
                        }
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                            {
                                nd->Read( rkey, fav );
                                nd->Read( key_, av );
                                av.CopyValuesOnly( fav );
                                nd->Store( key_, av );
                            }
                        }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( key_ ) == flag_output_ )
                    {
                        cit->Read( rkey, fav );
                        cit->Read( key_, av );
                        av.CopyValuesOnly( fav );
                        cit->Store( key_, av );
                    }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, key_ ) == flag_output_ )
                        {
                            cit->Read( ip, rkey, fav );
                            cit->Read( ip, key_, av );
                            av.CopyValuesOnly( fav );
                            cit->Store( ip, key_, av );
                        }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( s, 0U, rkey, fav );
                            cit->Read( s, 0U, key_, av );
                            av.CopyValuesOnly( fav );
                            cit->Store( s, 0U, key_, av );
                        }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( fac, 0U, rkey, fav );
                            cit->Read( fac, 0U, key_, av );
                            av.CopyValuesOnly( fav );
                            cit->Store( fac, 0U, key_, av );
                        }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                if ( sd.Status( key_ ) == flag_output_ )
                {
                    sd.Read( rkey, fav );
                    sd.Read( key_, av );
                    av.CopyValuesOnly( fav );
                    sd.Store( key_, av );
                }
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::AssignFlaggedArrayToArray",
                var_name_.c_str(), "Unsupported placement." );
    }
}



// ============================================================================
//  operator= (PropertyHandle)
// ============================================================================

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const PropertyHandle<dim>& other )
{
    if ( this == &other )
        return *this;

    const bool same_subdomain = ( subdomain_name_ == other.subdomain_name_ );
    const bool other_is_model = ( other.subdomain_name_ == "Model" );
    const bool this_is_model  = ( subdomain_name_ == "Model" );

    if ( !same_subdomain && !other_is_model )
    {
        if ( this_is_model )
            throw Exception( ERROR,
                "PropertyHandle::operator=(PropertyHandle)",
                var_name_.c_str(),
                ( string("Cannot assign sub-domain variable '")
                  + other.var_name_ + "' (sub-domain '"
                  + other.subdomain_name_
                  + "') into whole-model variable '"
                  + var_name_
                  + "': only part of the model would be overwritten. "
                  + "If this is intended, construct a whole-model handle "
                  + "for the source variable and assign from that."
                ).c_str() );
        else
            throw Exception( ERROR,
                "PropertyHandle::operator=(PropertyHandle)",
                var_name_.c_str(),
                ( string("Cannot assign between two different sub-domains: '")
                  + other.subdomain_name_ + "' -> '"
                  + subdomain_name_ + "'. "
                  + "Overlap between non-Model sub-domains is undefined."
                ).c_str() );
    }

    if ( !same_subdomain && other_is_model )
    {
        // Cross-subdomain assignment from whole model into a sub-domain
        // is only valid for NODE placement. FACE and INTER_FACE placements
        // exist exclusively on Boundary and SplitBoundary respectively.
        // ELEMENT and integration point placements exist exclusively on
        // Region. There is no meaningful cross-entity assignment between
        // these placement types.
        if ( key_.place != NODE || other.key_.place != NODE )
            throw Exception( ERROR,
                "PropertyHandle::operator=(PropertyHandle)",
                var_name_.c_str(),
                ( string("Cross-subdomain assignment from whole-model variable '")
                  + other.var_name_
                  + "' into sub-domain variable '"
                  + var_name_
                  + "' is only supported for NODE placement. "
                  + "FACE placement exists exclusively on Boundary subdomains. "
                  + "INTER_FACE placement exists exclusively on SplitBoundary "
                  + "subdomains. "
                  + "ELEMENT and integration point placements exist exclusively "
                  + "on Region subdomains. "
                  + "Source placement: "
                  + std::to_string( static_cast<int>( other.key_.place ) )
                  + ", destination placement: "
                  + std::to_string( static_cast<int>( key_.place ) )
                  + "."
                ).c_str() );
    }

    if ( key_.type != other.key_.type )
    {
    const bool this_is_array  = ( key_.type == ARRAY ||
                                   key_.type == FLAGGEDARRAY );
    const bool other_is_array = ( other.key_.type == ARRAY ||
                                   other.key_.type == FLAGGEDARRAY );

    // ------------------------------------------------------------------
    //  Array types are not interoperable with scalar/vector/tensor.
    // ------------------------------------------------------------------
    if ( this_is_array != other_is_array )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::operator=(PropertyHandle3)",
            var_name_.c_str(),
            ( string("Cannot assign between array and non-array "
                      "variable types. "
                      "Source type: ")
              + std::to_string( static_cast<int>( other.key_.type ) )
              + ", destination type: "
              + std::to_string( static_cast<int>( key_.type ) )
            ).c_str() );

    // ------------------------------------------------------------------
    //  ARRAY = FLAGGEDARRAY: copy values, discard per-element flags.
    // ------------------------------------------------------------------
    if ( key_.type == ARRAY && other.key_.type == FLAGGEDARRAY )
    {
        cerr << "\nPropertyHandle3::operator=: "
                "assigning FlaggedArrayVariable values to ArrayVariable "
                "(per-element flags discarded).\n";

        static std::atomic<uint32_t> tmp_counter_af{0};
        const string tmp_name = var_name_
                              + "__ph3_tmp_"
                              + std::to_string( tmp_counter_af.fetch_add(1) );

        const PLACEMENT dst     = key_.place;
        const PLACEMENT src     = other.key_.place;
        const bool needs_interp = ( dst != src );

        if ( needs_interp )
        {
            uint32_t array_size = 1U;
            WithSubdomain( [&]( auto& sd )
            {
                if ( !sd.CellVector().empty() )
                {
                    FlaggedArrayVariable fav;
                    sd.CellVector().front()->Read( other.key_, fav );
                    array_size = fav.Size();
                }
            });
            model_.CreateProperty( tmp_name.c_str(), tmp_name.c_str(),
                                   "SI", FLAGGEDARRAY, dst, array_size );
        }

        const char* rhs_name = needs_interp ? tmp_name.c_str()
                                            : other.VariableName();
        try
        {
            if ( needs_interp )
                Interpolate( other.VariableName(), rhs_name, src, dst );
            const csmp::Index rkey =
                model_.Database().StorageKey( rhs_name );
            AssignFlaggedArrayToArray( rkey );
        }
        catch ( ... )
        {
            if ( needs_interp )
                model_.DeleteProperty( tmp_name.c_str() );
            throw;
        }
        if ( needs_interp )
            model_.DeleteProperty( tmp_name.c_str() );

        IsWithinRange();
        return *this;
    }

    // ------------------------------------------------------------------
    //  FLAGGEDARRAY = ARRAY: copy values, set all flags to flag_output_.
    // ------------------------------------------------------------------
    if ( key_.type == FLAGGEDARRAY && other.key_.type == ARRAY )
    {
        cerr << "\nPropertyHandle3::operator=: "
                "assigning ArrayVariable values to FlaggedArrayVariable "
                "(all element flags set to OutputCondition()).\n";

        static std::atomic<uint32_t> tmp_counter_fa{0};
        const string tmp_name = var_name_
                              + "__ph3_tmp_"
                              + std::to_string( tmp_counter_fa.fetch_add(1) );

        const PLACEMENT dst     = key_.place;
        const PLACEMENT src     = other.key_.place;
        const bool needs_interp = ( dst != src );

        if ( needs_interp )
        {
            uint32_t array_size = 1U;
            WithSubdomain( [&]( auto& sd )
            {
                if ( !sd.CellVector().empty() )
                {
                    ArrayVariable av;
                    sd.CellVector().front()->Read( other.key_, av );
                    array_size = av.Size();
                }
            });
            model_.CreateProperty( tmp_name.c_str(), tmp_name.c_str(),
                                   "SI", ARRAY, dst, array_size );
        }

        const char* rhs_name = needs_interp ? tmp_name.c_str()
                                            : other.VariableName();
        try
        {
            if ( needs_interp )
                Interpolate( other.VariableName(), rhs_name, src, dst );
            const csmp::Index rkey =
                model_.Database().StorageKey( rhs_name );
            AssignFlaggedArrayToArray( rkey );
        }
        catch ( ... )
        {
            if ( needs_interp )
                model_.DeleteProperty( tmp_name.c_str() );
            throw;
        }
        if ( needs_interp )
            model_.DeleteProperty( tmp_name.c_str() );

        IsWithinRange();
        return *this;
    }

    // ------------------------------------------------------------------
    //  SCALAR = VECTOR (length) or SCALAR = TENSOR (determinant).
    // ------------------------------------------------------------------
    static std::atomic<uint32_t> tmp_counter{0};
    const string tmp_name = var_name_
                          + "__ph3_tmp_"
                          + std::to_string( tmp_counter.fetch_add(1) );

    const PLACEMENT dst     = key_.place;
    const PLACEMENT src     = other.key_.place;
    const bool needs_interp = ( dst != src );

    if ( needs_interp )
        model_.CreateProperty( tmp_name.c_str(), tmp_name.c_str(), "SI",
                               other.key_.type, dst, 1U );

    const char* rhs_name = needs_interp ? tmp_name.c_str()
                                        : other.VariableName();
    try
    {
        if ( needs_interp )
            Interpolate( other.VariableName(), rhs_name, src, dst );

        const csmp::Index rkey =
            model_.Database().StorageKey( rhs_name );

        if ( key_.type == SCALAR && other.key_.type == VECTOR )
        {
            cerr << "\nPropertyHandle3::operator=: "
                    "assigning vector length to scalar.\n";
            AssignVectorLengthToScalar( rkey );
        }
        else if ( key_.type == SCALAR && other.key_.type == TENSOR )
        {
            cerr << "\nPropertyHandle3::operator=: "
                    "assigning tensor determinant to scalar.\n";
            AssignTensorDetToScalar( rkey );
        }
        else
        {
            if ( needs_interp )
                model_.DeleteProperty( tmp_name.c_str() );
            throw csmp::Exception( ERROR,
                "PropertyHandle3::operator=(PropertyHandle3)",
                var_name_.c_str(),
                "No rule exists to assign between these variable types." );
        }
    }
    catch ( ... )
    {
        if ( needs_interp )
            model_.DeleteProperty( tmp_name.c_str() );
        throw;
    }
    if ( needs_interp )
        model_.DeleteProperty( tmp_name.c_str() );

    IsWithinRange();
    return *this;
}

    return ApplyBinaryOpAligned( other,
        []( double& a,                       double b )                { a = b; },
        []( VectorVariable<dim>& a,  const VectorVariable<dim>& b )   { a = b; },
        []( TensorVariable<dim>& a,  const TensorVariable<dim>& b )   { a = b; } );
}

// ============================================================================
//  operator= (uniform value / typed variable)
// ============================================================================

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const ScalarVariable& s )
{
    const double val = s();
    ApplyByType(
        [val]( double& sc )              { sc = val; },
        [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) = val; },
        [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                               for ( uint32_t j=0; j<dim; ++j ) ts(i,j) = val; }
    );
    IsWithinRange();
    return *this;
}

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const VectorVariable<dim>& v )
{
    switch ( key_.type )
    {
        case SCALAR:
        {
            const double len = v.Length();
            cerr << "\nPropertyHandle::operator=(VectorVariable): "
                    "assigning vector length to scalar.\n";
            ApplyScalar( [len]( double& sc ) { sc = len; } );
            break;
        }
        case VECTOR:
            ApplyVector( [&v]( VectorVariable<dim>& vc )
            {
                for ( uint32_t j = 0; j < dim; ++j ) vc(j) = v[j];
            });
            break;
        case TENSOR:
            cerr << "\nPropertyHandle::operator=(VectorVariable): "
                    "no rule to assign vector to tensor; nothing done.\n";
            break;
        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::operator=(VectorVariable)",
                var_name_.c_str(), "Unknown variable type." );
    }
    IsWithinRange();
    return *this;
}

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const TensorVariable<dim>& t )
{
    switch ( key_.type )
    {
        case SCALAR:
        {
            const double det = t.Determinant();
            cerr << "\nPropertyHandle::operator=(TensorVariable): "
                    "assigning determinant to scalar.\n";
            ApplyScalar( [det]( double& sc ) { sc = det; } );
            break;
        }
        case VECTOR:
            cerr << "\nPropertyHandle::operator=(TensorVariable): "
                    "no rule to assign tensor to vector; nothing done.\n";
            break;
        case TENSOR:
            ApplyTensor( [&t]( TensorVariable<dim>& ts )
            {
                for ( uint32_t i = 0; i < dim; ++i )
                    for ( uint32_t j = 0; j < dim; ++j )
                        ts(i,j) = t(i,j);
            });
            break;
        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::operator=(TensorVariable)",
                var_name_.c_str(), "Unknown variable type." );
    }
    IsWithinRange();
    return *this;
}



template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const std::vector<VectorVariable<dim>>& vc )
{
    if ( key_.type != VECTOR )
        throw csmp::Exception( ERROR,
            "PropertyHandle::operator=(vector<VectorVariable>)",
            var_name_.c_str(),
            "Type mismatch: handle does not hold a VECTOR variable." );

    FEM_Data<VectorVariable<dim>> var_data( key_.place, vc );
    model_.InputVariableFrom( model_.Database().Name( key_ ), var_data );
    return *this;
}



template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const std::vector<TensorVariable<dim>>& ts )
{
    if ( key_.type != TENSOR )
        throw csmp::Exception( ERROR,
            "PropertyHandle::operator=(vector<TensorVariable>)",
            var_name_.c_str(),
            "Type mismatch: handle does not hold a TENSOR variable." );

    FEM_Data<TensorVariable<dim>> var_data( key_.place, ts );
    model_.InputVariableFrom( model_.Database().Name( key_ ), var_data );
    return *this;
}


template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const ArrayVariable& av )
{
    if ( key_.type != ARRAY )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::operator=(ArrayVariable)",
            var_name_.c_str(),
            "Type mismatch: handle does not hold an ARRAY variable." );

    ApplyArray( [&av]( ArrayVariable& dest )
    {
        dest = av;
    });
    IsWithinRange();
    return *this;
}


template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( const FlaggedArrayVariable& fav )
{
    if ( key_.type != FLAGGEDARRAY )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::operator=(FlaggedArrayVariable)",
            var_name_.c_str(),
            "Type mismatch: handle does not hold a FLAGGEDARRAY variable." );

    ApplyFlaggedArray( [&fav]( FlaggedArrayVariable& dest )
    {
        dest = fav;
    });
    IsWithinRange();
    return *this;
}


// for array types
template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator=( double val )
{
    if ( key_.type == ARRAY )
    {
        ApplyArray( [val]( ArrayVariable& av )
        {
            av = val;  // ArrayVariable::operator=(double) sets all elements
        });
    }
    else if ( key_.type == FLAGGEDARRAY )
    {
        ApplyFlaggedArray( [val, this]( FlaggedArrayVariable& fav )
        {
            ApplyToFlaggedElements( fav, flag_output_,
                [val]( double& x ) { x = val; } );
        });
    }
    else
    {
        ApplyByType(
            [val]( double& sc )              { sc = val; },
            [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) = val; },
            [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                                   for ( uint32_t j=0; j<dim; ++j ) ts(i,j) = val; }
        );
    }
    IsWithinRange();
    return *this;
}



// ============================================================================
//  Scalar compound assignment  (+=, -=, *=, /=  with double)
// ============================================================================

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator+=( double val )
{
    if ( key_.type == ARRAY )
    {
        ApplyArray( [val]( ArrayVariable& av ) { av += val; });
    }
    else if ( key_.type == FLAGGEDARRAY )
    {
        ApplyFlaggedArray( [val, this]( FlaggedArrayVariable& fav )
        {
            ApplyToFlaggedElements( fav, flag_output_,
                [val]( double& x ) { x += val; } );
        });
    }
    else
    {
        ApplyByType(
            [val]( double& sc )              { sc += val; },
            [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) += val; },
            [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                                   for ( uint32_t j=0; j<dim; ++j ) ts(i,j) += val; }
        );
    }
    IsWithinRange();
    return *this;
}


template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator-=( double val )
{
    if ( key_.type == ARRAY )
    {
        ApplyArray( [val]( ArrayVariable& av ) { av -= val; });
    }
    else if ( key_.type == FLAGGEDARRAY )
    {
        ApplyFlaggedArray( [val, this]( FlaggedArrayVariable& fav )
        {
            ApplyToFlaggedElements( fav, flag_output_,
                [val]( double& x ) { x -= val; } );
        });
    }
    else
    {
      ApplyByType(
          [val]( double& sc )              { sc -= val; },
          [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) -= val; },
          [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                                 for ( uint32_t j=0; j<dim; ++j ) ts(i,j) -= val; }
      );
    }
    IsWithinRange();
    return *this;
}



template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator*=( double val )
{
    if ( key_.type == ARRAY )
    {
        ApplyArray( [val]( ArrayVariable& av ) { av *= val; });
    }
    else if ( key_.type == FLAGGEDARRAY )
    {
        ApplyFlaggedArray( [val, this]( FlaggedArrayVariable& fav )
        {
            ApplyToFlaggedElements( fav, flag_output_,
                [val]( double& x ) { x *= val; } );
        });
    }
    else
    {
      ApplyByType(
          [val]( double& sc )              { sc *= val; },
          [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) *= val; },
          [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                                 for ( uint32_t j=0; j<dim; ++j ) ts(i,j) *= val; }
      );
    }
    IsWithinRange();
    return *this;
}



template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator/=( double val )
{
    if ( val == 0.0 )
        throw csmp::Exception( ERROR,
            "PropertyHandle::operator/=(double)",
            var_name_.c_str(),
            "Division by zero." );

    if ( key_.type == ARRAY )
    {
        ApplyArray( [val]( ArrayVariable& av ) { av /= val; });
    }
    else if ( key_.type == FLAGGEDARRAY )
    {
        ApplyFlaggedArray( [val, this]( FlaggedArrayVariable& fav )
        {
            ApplyToFlaggedElements( fav, flag_output_,
                [val]( double& x ) { x /= val; } );
        });
    }
    else
    {
      ApplyByType(
          [val]( double& sc )              { sc /= val; },
          [val]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) vc(j) /= val; },
          [val]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                                 for ( uint32_t j=0; j<dim; ++j ) ts(i,j) /= val; }
      );
    }
    IsWithinRange();
    return *this;
}

// ============================================================================
//  PropertyHandle compound assignment  (+=, -=, *=, /=)
// ============================================================================

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator+=(
    const PropertyHandle<dim>& o )
{
    return ApplyBinaryOpAligned( o,
        []( double& a, double b )                                  { a += b; },
        []( VectorVariable<dim>& a, const VectorVariable<dim>& b ) { a += b; },
        []( TensorVariable<dim>& a, const TensorVariable<dim>& b ) { a += b; } );
}

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator-=(
    const PropertyHandle<dim>& o )
{
    return ApplyBinaryOpAligned( o,
        []( double& a, double b )                                  { a -= b; },
        []( VectorVariable<dim>& a, const VectorVariable<dim>& b ) { a -= b; },
        []( TensorVariable<dim>& a, const TensorVariable<dim>& b ) { a -= b; } );
}

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator*=(
    const PropertyHandle<dim>& o )
{
    return ApplyBinaryOpAligned( o,
        []( double& a, double b )                                  { a *= b; },
        []( VectorVariable<dim>& a, const VectorVariable<dim>& b ) { a *= b; },
        []( TensorVariable<dim>& a, const TensorVariable<dim>& b ) { a *= b; } );
}

template<uint32_t dim>
PropertyHandle<dim>& PropertyHandle<dim>::operator/=(
    const PropertyHandle<dim>& o )
{
    return ApplyBinaryOpAligned( o,
        []( double& a, double b )                                  { a /= b; },
        []( VectorVariable<dim>& a, const VectorVariable<dim>& b ) { a /= b; },
        []( TensorVariable<dim>& a, const TensorVariable<dim>& b ) { a /= b; } );
}

// ============================================================================
//  Output
// ============================================================================

template<uint32_t dim>
void PropertyHandle<dim>::Out() const
{
    cout << "\nPropertyHandle('" << subdomain_name_ << "')::Out:\n";
    cout << "  Variable : " << var_name_ << "\n";
    cout << "  Condition: " << parseStatus( flag_output_ ) << "\n";
    WithSubdomain( [&]( auto& sd )
    {
        sd.OutputVariableToScreen( var_name_.c_str() );
    });
}



template<uint32_t dim>
void PropertyHandle<dim>::Out( const char* filename ) const
{
    ofstream ofs( filename );
    if ( !ofs )
        throw csmp::Exception( ERROR,
            "PropertyHandle::Out(filename)",
            var_name_.c_str(),
            "Could not open output file." );

    ofs << "PropertyHandle('" << subdomain_name_ << "')::Out\n";
    ofs << "Variable : " << var_name_ << "\n";
    ofs << "Condition: " << parseStatus( flag_output_ ) << "\n";
    ofs.close();

    TextInterface().OutputDataAsTextColumns(
        subdomain_name_.c_str(), model_, filename, var_name_.c_str() );
}



// ============================================================================
//  Explicit instantiations
// ============================================================================

template class PropertyHandle<1U>;
template class PropertyHandle<2U>;
template class PropertyHandle<3U>;

} // namespace csmp

