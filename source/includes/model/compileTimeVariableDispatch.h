//
//  compileTimeVariableDispatch.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 5/7/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_COMPILE_TIME_VARIABLE_DISPATCH_H
#define CSMP_COMPILE_TIME_VARIABLE_DISPATCH_H

// static dispatch for LocalVariableStorage
#include "Index.h"
#include "LocalVariables.h"
#include "IntegrationPointVariables.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;

namespace lvsCompileTimeDispatch {

    //=====================================================================================================================================================
    // Assert Indexes of Integration Points

    template<uint32_t dim, template<uint32_t> class STOREE>
    void AssertFiniteVolumeIntegrationPointIndex(const STOREE<dim>* /* no IntegrationPointStoree */, size_t, size_t )
    {
        /* nothing to assert */
    }
    template<uint32_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::Element<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors() );
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }
    template<uint32_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::Face<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors());
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }
    template<uint32_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::InterFace<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors());
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }


    //=====================================================================================================================================================
    // New size of property data storage
    
    /**
        templatized function to associate LocalVariableStorage with Model, Region, Boundary and SplitBoundary classes which all have no integration points. 
    */
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<size_t, size_t> containerNewSize( const STOREE<dim>*, const LocalVariables& lv, const IntegrationPointVariables& )
      {
         const size_t dataDepth(lv.totalDataDepth);
         const size_t flagDepth(lv.totalFlagDepth);
         return std::make_pair(dataDepth, flagDepth);
      }
    
    /**
        templatized function for Element, Face and InterFace classes which have integration points.
    */
    template<uint32_t dim>
    std::pair<size_t, size_t> containerNewSize(const csmp::Element<dim>* e, const LocalVariables& lv, const IntegrationPointVariables& ipv)
    {
        const size_t dataDepth(lv.totalDataDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalDataDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalDataDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalDataDepth);
        const size_t flagDepth(lv.totalFlagDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalFlagDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalFlagDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerNewSize(const csmp::Face<dim>* e, const LocalVariables& lv, const IntegrationPointVariables& ipv)
    {
        const size_t dataDepth(lv.totalDataDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalDataDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalDataDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalDataDepth);
        const size_t flagDepth(lv.totalFlagDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalFlagDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalFlagDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerNewSize(const csmp::InterFace<dim>* e, const LocalVariables& lv, const IntegrationPointVariables& ipv)
    {
        const size_t dataDepth(lv.totalDataDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalDataDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalDataDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalDataDepth);
        const size_t flagDepth(lv.totalFlagDepth
            + e->IntegrationPoints()*ipv.ipvSimplex.totalFlagDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * ipv.ipvSector.totalFlagDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * ipv.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    //=====================================================================================================================================================
    // Total Depth of new data
    
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<size_t, size_t> containerTotalDataDepth( const STOREE<dim>* /* noIntegrationPointStoree */, const csmp::Index& idx )
    {
        const size_t dataDepth(idx.localVariables.totalDataDepth);
        const size_t flagDepth(idx.localVariables.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerTotalDataDepth(const csmp::Element<dim>* e, const csmp::Index& idx)
    {
        const size_t dataDepth(idx.localVariables.totalDataDepth
            + e->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);
        const size_t flagDepth(idx.localVariables.totalFlagDepth
            + e->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + e->IntegrationPointsPerFacet()  * e->Facets()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerTotalDataDepth(const csmp::Face<dim>* f, const csmp::Index& idx)
    {
        const size_t dataDepth(idx.localVariables.totalDataDepth
            + f->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + f->IntegrationPointsPerFacet()  * f->Facets()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);
        const size_t flagDepth(idx.localVariables.totalFlagDepth
            + f->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + f->IntegrationPointsPerFacet()  * f->Facets()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerTotalDataDepth(const csmp::InterFace<dim>* f, const csmp::Index& idx)
    {
        const size_t dataDepth(idx.localVariables.totalDataDepth
            + f->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + f->IntegrationPointsPerFacet()  * f->Facets()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);
        const size_t flagDepth(idx.localVariables.totalFlagDepth
            + f->IntegrationPoints()*idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + f->IntegrationPointsPerFacet()  * f->Facets()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    //=====================================================================================================================================================

    /** 
        Offset from beginning of data vector of the variable value record referred to by csmp::Index idx
        
                 Master template
    */
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<size_t, size_t> containerOffset( const STOREE<dim>* /* noIntegrationPointStoree */, const csmp::Index& idx )
    {
        const size_t dataOffset(idx.dataOffset);
        const size_t flagOffset(idx.flagOffset);
        return std::make_pair(dataOffset, flagOffset);
    }



    /**
        Specialization for Element class.
        
        Offset of sector or facet integration variable records from beginning of data vector in local variable storage associated with an element

        Function used only if index refers to a sector- or facet integration point variable.
    
        Method takes offset calculated in property database ignoring the number of integration points
        because these will vary from element to element type.
        This offset is now expanded for the element variable of interest 
        by multiplicating those other element integration point variables which precede the current 
        variable.
        
        @attention storage order is: all sector integration point variables for the first ip, then the second
        and so forth. Subsequently all facet integration point variables in the same fashion.

        @attention: this method is only used for the reading of variables that are stored in the LVS
        associated with Element objects
        
        @attention Index::offsetFactorSector is zero for sector integration point variables and 1 for facet integration point variables
    */
    template<uint32_t dim>
    std::pair<size_t, size_t> containerOffset( const csmp::Element<dim>* e, const csmp::Index& idx )
    {
        const size_t dataOffset(idx.dataOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalDataDepth
            + idx.offsetFactorSimplex * e->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.offsetFactorSector  * e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth);

        const size_t flagOffset(idx.flagOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalFlagDepth
            + idx.offsetFactorSimplex * e->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.offsetFactorSector  * e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth);

        return std::make_pair(dataOffset, flagOffset);
    }


    /// specialisation Face
    template<uint32_t dim>
    std::pair<size_t, size_t> containerOffset( const csmp::Face<dim>* f, const csmp::Index& idx )
    {
        const size_t dataOffset(idx.dataOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalDataDepth
            + idx.offsetFactorSimplex * f->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.offsetFactorSector  * f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth);
        const size_t flagOffset(idx.flagOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalFlagDepth
            + idx.offsetFactorSimplex * f->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.offsetFactorSector  * f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth);
        return std::make_pair(dataOffset, flagOffset);
    }



    /// specialisation InterFace
    template<uint32_t dim>
    std::pair<size_t, size_t> containerOffset(const csmp::InterFace<dim>* f, const csmp::Index& idx)
    {
        const size_t dataOffset(idx.dataOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalDataDepth
            + idx.offsetFactorSimplex * f->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.offsetFactorSector  * f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth);
        const size_t flagOffset(idx.flagOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalFlagDepth
            + idx.offsetFactorSimplex * f->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.offsetFactorSector  * f->IntegrationPointsPerSector() * f->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth);
        return std::make_pair(dataOffset, flagOffset);
    }

    //===================================================================================================================================================
    // Internal Cycles in new data container ( (1,1) for lvs,  (1,IPs) for ipvSimplex, ( Sectors, IPs ) for ipvSector and ( Facets, IPs ) for ipvFacet )
    
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<int, int> containerIPCycles(const STOREE<dim>* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(1, 1);
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycles(const csmp::Element<dim>* e, const csmp::Index& idx)
    {
        const long cycle1((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex /* 1 element */
            + idx.ipFactorSector  * e->Sectors()
            + idx.ipFactorFacet   * e->Facets());

        const long cycle2((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex * e->IntegrationPoints()
            + idx.ipFactorSector  * e->IntegrationPointsPerSector()
            + idx.ipFactorFacet   * e->IntegrationPointsPerFacet());

        return std::make_pair(cycle1, cycle2);
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycles(const csmp::Face<dim>* f, const csmp::Index& idx)
    {
        const long cycle1((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex /* 1 element */
            + idx.ipFactorSector  * f->Sectors()
            + idx.ipFactorFacet   * f->Facets());

        const long cycle2((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex * f->IntegrationPoints()
            + idx.ipFactorSector  * f->IntegrationPointsPerSector()
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet());

        return std::make_pair(cycle1, cycle2);
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycles(const csmp::InterFace<dim>* f, const csmp::Index& idx)
    {
        const long cycle1((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex /* 1 element */
            + idx.ipFactorSector  * f->Sectors()
            + idx.ipFactorFacet   * f->Facets());

        const long cycle2((1 - idx.ipFactorSimplex - idx.ipFactorSector - idx.ipFactorFacet)
            + idx.ipFactorSimplex * f->IntegrationPoints()
            + idx.ipFactorSector  * f->IntegrationPointsPerSector()
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet());

        return std::make_pair(cycle1, cycle2);
    }


    //=====================================================================================================================================================
    // Offset within first cycle ( 0 for lvs and ipvSimplex, TotalIpvDepth for ipvSector and ipvFacet )
    
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<int, int> containerIPCycle1Offset(const STOREE<dim>* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(static_cast<int>(0), static_cast<int>(0));
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycle1Offset(const csmp::Element<dim>* e, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSector  * e->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + idx.ipFactorFacet   * e->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);

        const long cycleFlagOffset(idx.ipFactorSector  * e->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + idx.ipFactorFacet   * e->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycle1Offset(const csmp::Face<dim>* f, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);

        const long cycleFlagOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<uint32_t dim>
    std::pair<int, int> containerIPCycle1Offset(const csmp::InterFace<dim>* f, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);

        const long cycleFlagOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    //=====================================================================================================================================================
    // Offset within second cycle ( 0 for lvs, TotalDepth for ipvSimplex, Depth for ipvSector and ipvFacet )
    
    template<uint32_t dim, template<uint32_t> class STOREE>
    std::pair<size_t, size_t> containerIPCycle2Offset(const STOREE<dim>* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(static_cast<int>(0), static_cast<int>(0));
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerIPCycle2Offset(const csmp::Element<dim>*, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.ipFactorSector  * idx.dataDepth
            + idx.ipFactorFacet   * idx.dataDepth);

        const long cycleFlagOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.ipFactorSector  * idx.flagDepth
            + idx.ipFactorFacet   * idx.flagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerIPCycle2Offset(const csmp::Face<dim>*, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.ipFactorSector  * idx.dataDepth
            + idx.ipFactorFacet   * idx.dataDepth);

        const long cycleFlagOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.ipFactorSector  * idx.flagDepth
            + idx.ipFactorFacet   * idx.flagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<uint32_t dim>
    std::pair<size_t, size_t> containerIPCycle2Offset(const csmp::InterFace<dim>*, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalDataDepth
            + idx.ipFactorSector  * idx.dataDepth
            + idx.ipFactorFacet   * idx.dataDepth);

        const long cycleFlagOffset(idx.ipFactorSimplex * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.ipFactorSector  * idx.flagDepth
            + idx.ipFactorFacet   * idx.flagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }
    
} // lvsCompileTimeDispatch





} // end csmp

#endif /* CSMP_COMPILE_TIME_VARIABLE_DISPATCH_H */
