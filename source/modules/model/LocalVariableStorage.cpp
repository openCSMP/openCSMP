#include "LocalVariableStorage.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"

#define VARIABLE_STORAGE_DEBUG

namespace csmp {

namespace lvsCompileTimeDispatch{

    //=====================================================================================================================================================
    // Assert Indexes of Integration Points

    template<class STOREE>
    void AssertFiniteVolumeIntegrationPointIndex(const STOREE* /* no IntegrationPointStoree */, size_t, size_t )
    {
        /* nothing to assert */
    }

    template<size_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::Element<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors() );
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }
    template<size_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::Face<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors());
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }
    template<size_t dim>
    void AssertFiniteVolumeIntegrationPointIndex(const csmp::InterFace<dim>* e, size_t sector_or_facet, size_t ip)
    {
        assert(sector_or_facet < e->Facets() || sector_or_facet < e->Sectors());
        assert(ip < e->IntegrationPointsPerSector() || ip < e->IntegrationPointsPerFacet());
    }


    //=====================================================================================================================================================
    // New size of property data storage

    template<class STOREE>
    std::pair<size_t, size_t> containerNewSize(const STOREE* /* noIntegrationPointStoree */, const LocalVariables& lv, const IntegrationPointVariables& )
    {
        const size_t dataDepth(lv.totalDataDepth);
        const size_t flagDepth(lv.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<size_t dim>
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

    template<size_t dim>
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

    template<size_t dim>
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
    template<class STOREE>
    std::pair<size_t, size_t> containerTotalDataDepth(const STOREE* /* noIntegrationPointStoree */, const csmp::Index& idx)
    {
        const size_t dataDepth(idx.localVariables.totalDataDepth);
        const size_t flagDepth(idx.localVariables.totalFlagDepth);
        return std::make_pair(dataDepth, flagDepth);
    }

    template<size_t dim>
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

    template<size_t dim>
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

    template<size_t dim>
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
    */
    template<class STOREE>
    std::pair<size_t, size_t> containerOffset( const STOREE* /* noIntegrationPointStoree */, const csmp::Index& idx )
    {
        const size_t dataOffset(idx.dataOffset);
        const size_t flagOffset(idx.flagOffset);
        return std::make_pair(dataOffset, flagOffset);
    }



    /**
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
    template<size_t dim>
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

/* ORIGINAL (worked for sector integration point variables only)
    template<size_t dim>
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
*/


/* SKM EXTENDED VERSION
    template<size_t dim>
    std::pair<size_t, size_t> containerOffset( const csmp::Element<dim>* e, const csmp::Index& idx )
    {
        const size_t dataOffset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                   idx.dataOffset + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalDataDepth
                                 + idx.offsetFactorSimplex * e->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalDataDepth
                                 + idx.offsetFactorSector  * e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth
                                                                          : // FACET_INTEGRATION_POINT
                                   idx.dataOffset + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalDataDepth
                                 + idx.offsetFactorSimplex * e->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalDataDepth
                                 + idx.offsetFactorSector  * e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalDataDepth;
                                 + idx.offsetFactorSector+1 * e->IntegrationPointsPerFacet() * e->Facets() * idx.integrationPointVariables.ipvFacet.totalDataDepth;

        const size_t flagOffset(idx.flagOffset
            + (idx.ipFactorSimplex + idx.ipFactorSector + idx.ipFactorFacet) * idx.localVariables.totalFlagDepth
            + idx.offsetFactorSimplex * e->IntegrationPoints() * idx.integrationPointVariables.ipvSimplex.totalFlagDepth
            + idx.offsetFactorSector  * e->IntegrationPointsPerSector() * e->Sectors() * idx.integrationPointVariables.ipvSector.totalFlagDepth);

        return std::make_pair(dataOffset, flagOffset);
    }
*/


    template<size_t dim>
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

    template<size_t dim>
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
    template<class STOREE>
    std::pair<int, int> containerIPCycles(const STOREE* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(1, 1);
    }

    template<size_t dim>
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

    template<size_t dim>
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

    template<size_t dim>
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
    template<class STOREE>
    std::pair<int, int> containerIPCycle1Offset(const STOREE* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(static_cast<int>(0), static_cast<int>(0));
    }

    template<size_t dim>
    std::pair<int, int> containerIPCycle1Offset(const csmp::Element<dim>* e, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSector  * e->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + idx.ipFactorFacet   * e->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);

        const long cycleFlagOffset(idx.ipFactorSector  * e->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + idx.ipFactorFacet   * e->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<size_t dim>
    std::pair<int, int> containerIPCycle1Offset(const csmp::Face<dim>* f, const csmp::Index& idx)
    {
        const long cycleDataOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalDataDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalDataDepth);

        const long cycleFlagOffset(idx.ipFactorSector  * f->IntegrationPointsPerSector() * idx.integrationPointVariables.ipvSector.totalFlagDepth
            + idx.ipFactorFacet   * f->IntegrationPointsPerFacet()  * idx.integrationPointVariables.ipvFacet.totalFlagDepth);

        return std::make_pair(cycleDataOffset, cycleFlagOffset);
    }

    template<size_t dim>
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
    template<class STOREE>
    std::pair<size_t, size_t> containerIPCycle2Offset(const STOREE* /* noIntegrationPointStoree */, const csmp::Index&)
    {
        return std::make_pair(static_cast<int>(0), static_cast<int>(0));
    }

    template<size_t dim>
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

    template<size_t dim>
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

    template<size_t dim>
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


















#ifndef NDEBUG
/// We check in debug only
template<size_t dim, class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertPlacement(const csmp::Index& idx) const
{
    assert(idx.place == static_cast<const STOREE*>(this)->Placement());
}


/// We check in debug only
template<size_t dim, class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertIntegrationPointPlacement(const csmp::Index& idx) const
{
    assert(idx.place == ELEMENT_INTEGRATION_POINT || idx.place == FACE_INTEGRATION_POINT || idx.place == INTER_FACE_INTEGRATION_POINT);
}

/// We check in debug only
template<size_t dim, class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertFiniteVolumeIntegrationPointPlacement(const csmp::Index& idx) const
{
    assert(idx.place == SECTOR_INTEGRATION_POINT || idx.place == FACET_INTEGRATION_POINT ||
        idx.place == FACE_SECTOR_INTEGRATION_POINT || idx.place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
        idx.place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT);
}

/// We store in debug only
template<size_t dim, class STOREE>
void LocalVariableStorage<dim, STOREE>::StoreLocalState(const LocalVariables& lv)
{
    data_.scalars = lv.scalars;
    data_.vectors = lv.vectors;
    data_.tensors = lv.tensors;
    data_.arrays = lv.arrayCount;
    data_.arrayLength = lv.arrayLength;
    data_.flaggedArrays = lv.flaggedArrayCount;
    data_.flaggedArrayLength = lv.flaggedArrayLength;
}
#endif


/// Ctor for plain(empty) variable storage
template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage()
    : data_()
{
}


/// Ctor for storees without IntegrationPointVariables
template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage(const LocalVariables& lv)
    : data_()
{
    ResizePropertyStorage(lv);
}


/// Goto overload for storees with IntegrationPointVariables @attention Will fail AT COMPILE TIME (very nice...) if used for an storee without integration points
template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage(const LocalVariables& lv, const IntegrationPointVariables& ipv)
    : data_()
{
    ResizePropertyStorage(lv, ipv);
}


template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>::~LocalVariableStorage()
{
}


template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage(const LocalVariableStorage<dim,STOREE>& ps )
    : data_(ps.data_)
{
}


template<size_t dim, class STOREE>
LocalVariableStorage<dim, STOREE>& LocalVariableStorage<dim, STOREE>::operator=(const LocalVariableStorage& ps)
{
    if (&ps != this)
        data_.data = ps.data_.data;
    return *this;
}



/// Goto overload for storees without IntegrationPointVariables
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv )
  {
    ResizePropertyStorage( lv.totalDataDepth , lv.totalFlagDepth );

#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }


/// Goto overload for storees with IntegrationPointVariables
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv, const IntegrationPointVariables& ipv )
  {
    const std::pair<size_t, size_t> newContainerSize = lvsCompileTimeDispatch::containerNewSize(static_cast<const STOREE*>(this), lv, ipv);
    ResizePropertyStorage(newContainerSize.first, newContainerSize.second);

#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }


/// Resizes the storage preserving original values if any (this is where all ResizePropertyStorage end up)
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( size_t newDataComponentCount, size_t newFlagComponentCount ) 
  {
    // needed?
    if( data_.flags.size() == newFlagComponentCount && data_.data.size() == newDataComponentCount )
      return;

    // resizing
    data_.flags.resize( newFlagComponentCount, ANY );
    data_.data.resize( newDataComponentCount, std::numeric_limits<double64>::quiet_NaN() );

    // trimming excess capacity
    std::vector<VARIABLE_FLAG>( data_.flags ).swap( data_.flags );
    std::vector<double64>( data_.data ).swap( data_.data );
  }


/// Creates an empty space in the storage preserving the original variable values. LocalVariables in Index have to reflect NEW state (incl added)
/**
@todo (2-F) Integration point properties not supported
*/
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::AddProperty( const csmp::Index& prop_key )
  {
    /// @todo (2-F) Asserts missing

    // Resize to new state (this already supports ip vars)
    const std::pair<size_t,size_t> newTotaDataDepth = lvsCompileTimeDispatch::containerTotalDataDepth( static_cast<const STOREE*>(this), prop_key );
    ResizePropertyStorage( newTotaDataDepth.first, newTotaDataDepth.second );
    const size_t dataSize( data_.data.size() );
    const size_t flagSize( data_.flags.size() );

    /// Roman, 2013: without ipvs support

//    // Offset, End and Depth
//    const size_t dataOffset (prop_key.dataOffset);
//    const size_t dataDepth  (prop_key.dataDepth);
//    const size_t dataEnd    (dataOffset+dataDepth-1);
//    // DataDepth
//    const size_t flagOffset (prop_key.flagOffset);
//    const size_t flagDepth  (prop_key.flagDepth);
//    const size_t flagEnd    (flagOffset+flagDepth-1);
//    // Moving data (same for all types)
//    for ( size_t i=dataSize-1U; i>dataEnd; --i )
//      data_.data[i] = data_.data[i-dataDepth];
//    // and flags
//    for ( size_t i=flagSize-1U; i>flagEnd; --i )
//      data_.flags[i] = data_.flags[i-flagDepth];
//    // Adding plain variable
//    for ( size_t i=dataOffset; i<=dataEnd; ++i )
//      data_.data[i] = std::numeric_limits<double>::quiet_NaN();
//    for ( size_t i=flagOffset; i<=flagEnd; ++i )
//      data_.flags[i] = ANY;

    /// Roman, 2013: with ipvs support

    // Offset
    const std::pair<size_t,size_t> newOffset = lvsCompileTimeDispatch::containerOffset( static_cast<const STOREE*>(this), prop_key );
    const size_t dataOffset         ( newOffset.first       );
    const size_t flagOffset         ( newOffset.second      );

    // DataDepth
    const size_t dataDepth          ( prop_key.dataDepth    );
    const size_t flagDepth          ( prop_key.flagDepth    );

    // Cycles
    const std::pair<int,int> newCycles = lvsCompileTimeDispatch::containerIPCycles( static_cast<const STOREE*>(this), prop_key );
    const int nipCycles1            ( newCycles.first       );
    const int nipCycles2            ( newCycles.second      );

    const std::pair<int,int> newCycles1Data = lvsCompileTimeDispatch::containerIPCycle1Offset( static_cast<const STOREE*>(this), prop_key );
    const int ipCycle1DataOffset    ( newCycles1Data.first  );
    const int ipCycle1FlagOffset    ( newCycles1Data.second );

    const std::pair<int,int> newCycles2Data = lvsCompileTimeDispatch::containerIPCycle2Offset( static_cast<const STOREE*>(this), prop_key );
    const int ipCycle2DataOffset    ( newCycles2Data.first  );
    const int ipCycle2FlagOffset    ( newCycles2Data.second );

    // Data and Flag Bounds
    std::vector<std::vector<std::pair<size_t,size_t> > > dataBounds( nipCycles1, std::vector<std::pair<size_t,size_t> >(nipCycles2+1,std::pair<size_t,size_t>(0,0)) );
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
            dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
        }
        dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
    }
    dataBounds[nipCycles1-1][nipCycles2].first  = dataSize;
    dataBounds[nipCycles1-1][nipCycles2].second = dataSize;

    std::vector<std::vector<std::pair<size_t,size_t> > > flagBounds( nipCycles1, std::vector<std::pair<size_t,size_t> >(nipCycles2+1,std::pair<size_t,size_t>(0,0)) );
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
            flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
        }
        flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
    }
    flagBounds[nipCycles1-1][nipCycles2].first  = flagSize;
    flagBounds[nipCycles1-1][nipCycles2].second = flagSize;

    // Moving data (same for all types)
    for (int cycle1=nipCycles1-1U; cycle1>=0; --cycle1  ){
        for (int cycle2=nipCycles2-1U; cycle2>=0; --cycle2  ){
            const size_t dataStart  ( dataBounds[ cycle1 ][ cycle2     ].second - 1 );
            const size_t dataEnd    ( dataBounds[ cycle1 ][ cycle2 + 1 ].second - 1 );
            for ( size_t i=dataEnd; i>dataStart; --i )
              data_.data[i] = data_.data[i-(cycle1+1)*(cycle2+1)*dataDepth];
        }
    }
    // Moving flags (same for all types)
    for (int cycle1=nipCycles1-1U; cycle1>=0; --cycle1  ){
        for (int cycle2=nipCycles2-1U; cycle2>=0; --cycle2  ){
            const size_t flagsStart ( flagBounds[ cycle1 ][ cycle2     ].second - 1 );
            const size_t flagsEnd   ( flagBounds[ cycle1 ][ cycle2 + 1 ].second - 1 );
            for ( size_t i=flagsEnd; i>flagsStart; --i )
                data_.flags[i] = data_.flags[i-(cycle1+1)*(cycle2+1)*flagDepth];
        }
    }

    // Adding non-initialized data of new variable
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( size_t i=dataBounds[cycle1][cycle2].first; i<dataBounds[cycle1][cycle2].second; ++i )
                data_.data[i] = std::numeric_limits<double64>::quiet_NaN();
    // Adding non-initialized flags of new variable
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( size_t i=flagBounds[cycle1][cycle2].first; i<flagBounds[cycle1][cycle2].second; ++i )
                data_.flags[i] = ANY;


#ifndef NDEBUG
    StoreLocalState(prop_key.localVariables);
#endif
  } // end AddProperty


/// Deletes data and flags of property LocalVariables in Index have to reflect CURRENT state (incl property to be removed)
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::DeleteProperty( const csmp::Index& prop_key )
  {
    /// @todo (2-F) Asserts missing
    const size_t flagSize( data_.flags.size() );
    const size_t dataSize( data_.data.size() );

    /// Roman, 2013: without ipvs support

//    // from here on we need to implement ip var support
//    // shifting
//    for ( size_t i(prop_key.flagOffset); i<(flagSize-prop_key.flagDepth); ++i )
//      data_.flags[i] = data_.flags[i+prop_key.flagDepth];
//    for ( size_t i(prop_key.dataOffset); i<(dataSize-prop_key.dataDepth); ++i )
//      data_.data[i]  = data_.data[i+prop_key.dataDepth];

//    // trim excessive
//    ResizePropertyStorage( dataSize-prop_key.dataDepth, flagSize-prop_key.flagDepth );

//#ifndef NDEBUG
//    if( prop_key.type       == SCALAR )
//        --data_.scalars;
//    else if( prop_key.type  == VECTOR )
//        --data_.vectors;
//    else if( prop_key.type  == TENSOR )
//        --data_.tensors;
//    else if( prop_key.type  == ARRAY )
//      {
//        --data_.arrays;
//        data_.arrayLength -= prop_key.dataDepth;
//      }
//    else if( prop_key.type  == FLAGGEDARRAY )
//      {
//        --data_.flaggedArrays;
//        data_.flaggedArrayLength -= prop_key.dataDepth;
//      }

//#endif

        /// Roman, 2013: with ipvs support

        // Offset
        const std::pair<size_t,size_t> newOffset = lvsCompileTimeDispatch::containerOffset( static_cast<const STOREE*>(this), prop_key );
        const size_t dataOffset         ( newOffset.first       );
        const size_t flagOffset         ( newOffset.second      );

        // DataDepth
        const size_t dataDepth          ( prop_key.dataDepth    );
        const size_t flagDepth          ( prop_key.flagDepth    );

        // Cycles
        const std::pair<int,int> newCycles = lvsCompileTimeDispatch::containerIPCycles( static_cast<const STOREE*>(this), prop_key );
        const int nipCycles1            ( newCycles.first       );
        const int nipCycles2            ( newCycles.second      );

        const std::pair<int,int> newCycles1Data = lvsCompileTimeDispatch::containerIPCycle1Offset( static_cast<const STOREE*>(this), prop_key );
        const int ipCycle1DataOffset    ( newCycles1Data.first  );
        const int ipCycle1FlagOffset    ( newCycles1Data.second );

        const std::pair<int,int> newCycles2Data = lvsCompileTimeDispatch::containerIPCycle2Offset( static_cast<const STOREE*>(this), prop_key );
        const int ipCycle2DataOffset    ( newCycles2Data.first  );
        const int ipCycle2FlagOffset    ( newCycles2Data.second );

        // Data and Flag Bounds
        std::vector<std::vector<std::pair<size_t,size_t> > > dataBounds( nipCycles1, std::vector<std::pair<size_t,size_t> >(nipCycles2+1,std::pair<size_t,size_t>(0,0)) );
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
                dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
            }
            dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
            dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        }
        dataBounds[nipCycles1-1][nipCycles2].first  = dataSize;
        dataBounds[nipCycles1-1][nipCycles2].second = dataSize;

        std::vector<std::vector<std::pair<size_t,size_t> > > flagBounds( nipCycles1, std::vector<std::pair<size_t,size_t> >(nipCycles2+1,std::pair<size_t,size_t>(0,0)) );
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
                flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
            }
            flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
            flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        }
        flagBounds[nipCycles1-1][nipCycles2].first  = flagSize;
        flagBounds[nipCycles1-1][nipCycles2].second = flagSize;


        // Shifting data (same for all types)
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const size_t dataStart  ( dataBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*dataDepth );
                const size_t dataEnd    ( dataBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*dataDepth);
                for ( size_t i=dataStart; i<dataEnd; ++i )
                  data_.data[i] = data_.data[i+(cycle1*nipCycles2+cycle2+1)*dataDepth];
            }
        }
        // Shifting flags (same for all types)
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const size_t flagsStart ( flagBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*flagDepth );
                const size_t flagsEnd   ( flagBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*flagDepth);
                for ( size_t i=flagsStart; i<flagsEnd; ++i )
                    data_.flags[i] = data_.flags[i+(cycle1*nipCycles2+cycle2+1)*flagDepth];
            }
        }

        // trim excessive
        const size_t newDataSize ( dataSize - nipCycles1*nipCycles2*dataDepth );
        const size_t newFlagSize ( flagSize - nipCycles1*nipCycles2*flagDepth );
        ResizePropertyStorage( newDataSize, newFlagSize );

     #ifndef NDEBUG
        if( ( prop_key.ipFactorSimplex + prop_key.ipFactorSector + prop_key.ipFactorFacet ) == 0 )
        {
            if( prop_key.type       == SCALAR )
                --data_.scalars;
            else if( prop_key.type  == VECTOR )
                --data_.vectors;
            else if( prop_key.type  == TENSOR )
                --data_.tensors;
            else if( prop_key.type  == ARRAY )
              {
                --data_.arrays;
                data_.arrayLength -= prop_key.dataDepth;
              }
            else if( prop_key.type  == FLAGGEDARRAY )
              {
                --data_.flaggedArrays;
                data_.flaggedArrayLength -= prop_key.dataDepth;
              }
        }
     #endif

  } // end DeleteProperty


template<size_t dim,class STOREE>
bool LocalVariableStorage<dim,STOREE>::EmptyLVS() const
  {
    return data_.data.empty();
  }


/// @todo (2-F) LVS-ARRAY, FLAGGEDARRAY, IPVs
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::OutLVS() const
  {
    #ifndef NDEBUG
    std::cout <<"\nLocalVariableStorage<" << dim << ">::Out: ";
    std::cout <<"\n\tstored scalar variables: ";
    if ( data_.scalars > 0U ) {
      csmp::Index  idx(SCALAR,MODEL,0U);
      while ( idx.index < data_.scalars ) {
        double64 sc = Read( idx );
        std::cout << std::endl <<"\t\t"<< sc;
        idx.index++;
        }
      }

    if ( data_.vectors > 0U ) {
      std::cout <<"\n\n\tstored vector variables: ";
      csmp::Index  idx(VECTOR,MODEL,0U);
      VectorVariable<dim>  vc;
      while ( idx.index < data_.vectors ) {
        Read( idx, vc );
        std::cout << std::endl <<"\t\t"<< vc;
        idx.index++;
        }
      }

    if ( data_.flags.size() > data_.scalars + data_.vectors * dim ) {
      std::cout <<"\n\n\tstored tensor variables: ";
      const size_t  tensors( data_.tensors );
      csmp::Index  idx(TENSOR,MODEL,0U);
      TensorVariable<dim>  ts;
      while ( idx.index < tensors ) {
        Read( idx, ts );
        std::cout << std::endl <<"\t\t"<< ts;
        idx.index++;
        }
      }
    std::cout << std::endl;
    #endif
  }






















// ===============
// LOCAL VARIABLES
// ===============

/// Scalar variable value
template<size_t dim,class STOREE>
double64 LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx ) const  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.dataOffset < data_.data.size() );
#endif
    return data_.data[idx.dataOffset];
 }
 
 
/// Scalar variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ScalarVariable& sc ) const  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.dataOffset < data_.data.size() );
 assert( idx.flagOffset < data_.flags.size() );
#endif
    sc.Flag() = data_.flags[idx.flagOffset];
    sc        = data_.data[idx.dataOffset];
 }


/// Scalar variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const ScalarVariable& sc )  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
 assert( idx.dataOffset < data_.data.size() );
#endif
    data_.flags[idx.flagOffset] = sc.Flag();
    data_.data[idx.dataOffset]  = sc();
 }


/// Scalar & Array variable flag
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx ) const 
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR || idx.type == ARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    return data_.flags[idx.flagOffset];
 }
 

/// Vector, Tensor, FlaggedArray variable flag
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, size_t i ) const  
 {
 #ifndef NDEBUG
  AssertPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    return data_.flags[ idx.flagOffset + i ];
 }
 
 
/// Scalar & Array variable flag
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, VARIABLE_FLAG flag ) 
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR  || idx.type == ARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    data_.flags[idx.flagOffset] = flag;
 }

 
/// Vector,Tensor, FlaggedArray variable flag
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, size_t i, VARIABLE_FLAG flag ) 
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    data_.flags[ idx.flagOffset + i ] = flag;
 }
 

/// Vector variable 
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const VectorVariable<dim>& vc )  
 {
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( size_t i(0); i<dim; ++i ) {
         data_.flags[ idx.flagOffset+i ] = vc.Flag(i);
         data_.data[ idx.dataOffset+i ]  = vc[i];
      }
 }


/// Vector variable 
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, VectorVariable<dim>& vc ) const  
 {
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( size_t i(0); i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ idx.flagOffset+i ];
         vc(i)      = data_.data[ idx.dataOffset+i ];
      }
 }
 

/// Tensor variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const TensorVariable<dim>& ts )
 {
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == TENSOR );
  assert( idx.index < data_.tensors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.flagOffset+dim-1) < data_.flags.size() );
  assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
    const size_t dataOffset(idx.dataOffset);
    const size_t flagOffset(idx.flagOffset);
    for ( size_t i=0U; i<dim; i++ ) 
      {
        data_.flags[ flagOffset+i ] = ts.Flag(i);
        for ( size_t j=0U; j<dim; j++ )
          data_.data[ dataOffset+i*dim+j ] = ts(i,j);
      }
 }


/// Tensor variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, TensorVariable<dim>& ts ) const  
 {
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == TENSOR );
 assert( idx.index < data_.tensors );
#endif  
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
   const size_t dataOffset(idx.dataOffset);
   const size_t flagOffset(idx.flagOffset);
   for ( size_t i=0U; i<dim; i++ ) 
     {
       ts.Flag(i) = data_.flags[ flagOffset+i ] ;
       for ( size_t j=0U; j<dim; j++ )
         ts(i,j) = data_.data[ dataOffset+i*dim+j ];
     }
 }


/// Array variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const ArrayVariable& av )
  {
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const size_t data_offset( idx.dataOffset );
    const size_t flags_offset( idx.flagOffset );
    const size_t arraySize( idx.dataDepth );
    for( size_t i(0); i < arraySize; ++i )
      data_.data[ data_offset   + i ] = av[i];
    data_.flags[ flags_offset] = av.Flag();
  }


/// Array variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ArrayVariable& av ) const
  {
    av.Resize( idx.dataDepth );
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const size_t data_offset( idx.dataOffset );
    const size_t flags_offset( idx.flagOffset );
    const size_t arraySize( idx.dataDepth );
    for( size_t i(0); i < arraySize; ++i )
      av(i)     = data_.data[ data_offset +  i ];
    av.Flag()= data_.flags[ flags_offset];
}

/// FlaggedArray variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );
#endif
    const size_t data_offset( idx.dataOffset );
    const size_t flags_offset( idx.flagOffset );
    const size_t arraySize( idx.dataDepth );
    for( size_t i(0); i < arraySize; ++i )
    {
      data_.data[ data_offset   + i ] = av[i];
      data_.flags[ flags_offset + i ] = av.Flag(i);
    }
  }


/// FlaggedArray variable
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
  av.Resize( idx.dataDepth );
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );
#endif
    const size_t data_offset( idx.dataOffset );
    const size_t flags_offset( idx.flagOffset );
    const size_t arraySize( idx.dataDepth );
    for( size_t i(0); i < arraySize; ++i )
    {
      av(i)     = data_.data[ data_offset +  i ];
      av.Flag(i)= data_.flags[ flags_offset + i ];
    }
}


/**
    @todo replace this super wasteful method
*/
template<size_t dim,class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( const csmp::Index& idx, 
                                                      double64 vmin, double64 vmax ) const
  {
#ifndef NDEBUG
  AssertPlacement( idx );
#endif
    if ( idx.type == SCALAR ) {
      const double64 val = Read( idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim>  vc;
      Read( idx, vc );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts;
      Read( idx, ts );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      ArrayVariable  av;
      Read( idx, av );
      return av.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }

    std::cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }
























// ===================================
// SIMPLEX INTEGRATION POINT VARIABLES
// ===================================

#include "LocalVariableStorageIndexArithmetic.h"


/// Scalar variable value at integration point
template<size_t dim,class STOREE>
double64 LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx ) const  
  {
    const size_t offset(DATA_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
#endif

    return data_.data[offset];
  }


/// Scalar variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx, ScalarVariable& sc ) const  
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif

    sc.Flag() = data_.flags[flagOffset];
    sc        = data_.data[offset];
  }


/// Scalar variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t ip, const csmp::Index& idx, const ScalarVariable& sc )  
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = sc.Flag();
    data_.data[offset] = sc();
  }


/// Scalar & Array variable flag at integration point
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( size_t ip, const csmp::Index& idx ) const 
  {
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif

    return data_.flags[flagOffset];
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( size_t ip, const csmp::Index& idx, size_t i ) const  
  {
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
    AssertIntegrationPointPlacement(idx);
    assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
    assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
    assert( flagOffset+i < data_.flags.size() );
#endif

    return data_.flags[ flagOffset+i ];
  }


/// Scalar & Array variable flag at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( size_t ip, const csmp::Index& idx, VARIABLE_FLAG flag ) 
  {
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR  || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = flag;
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( size_t ip, const csmp::Index& idx, size_t i, VARIABLE_FLAG flag ) 
  {
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif

    data_.flags[ flagOffset+i ] = flag;
  }


/// Vector variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t ip, const csmp::Index& idx, const VectorVariable<dim>& vc )  
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

  for ( size_t i(0); i<dim; ++i ) {
    data_.flags[ flagOffset+i ] = vc.Flag(i);
    data_.data[ offset+i ]  = vc[i];
    }
  }


/// Vector variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx, VectorVariable<dim>& vc ) const  
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

    for ( size_t i(0); i<dim; ++i ) {
      vc.Flag(i) = data_.flags[ flagOffset+i ];
      vc(i)      = data_.data[ offset+i ];
      }
  }


/// Tensor variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t ip, const csmp::Index& idx, const TensorVariable<dim>& ts )
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

  for ( size_t i=0U; i<dim; i++ ) 
    {
    data_.flags[ flagOffset+i ] = ts.Flag(i);
    for ( size_t j=0U; j<dim; j++ )
      data_.data[ offset+i*dim+j ] = ts(i,j);
    }
  }


/// Tensor variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx, TensorVariable<dim>& ts ) const  
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

  for ( size_t i=0U; i<dim; i++ ) 
    {
    ts.Flag(i) = data_.flags[flagOffset+i];
    for ( size_t j=0U; j<dim; j++ )
      ts(i,j) = data_.data[ offset+i*dim+j ];
    }
  }

/// Array variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t ip, const csmp::Index& idx, const ArrayVariable& av )
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);
    const size_t arraySize( idx.dataDepth );

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif

    for( size_t i(0); i < arraySize; ++i )
      data_.data[ offset     + i ] = av[i];
    data_.flags[flagOffset] = av.Flag();
  }


/// Array variables
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx, ArrayVariable& av ) const
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);
    const size_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif

    for( size_t i(0); i < arraySize; ++i )
      av(i) = data_.data[ offset+i ];

    av.Flag( data_.flags[flagOffset] );
  }

/// FlaggedArray variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t ip, const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);
    const size_t arraySize( idx.dataDepth );

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif

    for( size_t i(0); i < arraySize; ++i )
    {
      data_.data[ offset     + i ] = av[i];
      data_.flags[flagOffset + i ] = av.Flag(i);
    }
  }


/// FlaggedArray variables
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t ip, const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
    const size_t offset(DATA_OFFSET_IP);
    const size_t flagOffset(FLAG_OFFSET_IP);
    const size_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif

    for( size_t i(0); i < arraySize; ++i )
    {
      av(i) = data_.data[ offset+i ];
      av.Flag( i, data_.flags[flagOffset+i] );
    }
  }


template<size_t dim,class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( size_t ip, const csmp::Index& idx, 
                                                      double64 vmin, double64 vmax ) const
  {
#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
#endif

    if ( idx.type == SCALAR ) {
      const double64 val = Read( ip, idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim>  vc;
      Read( ip, idx, vc );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts;
      Read( ip, idx, ts );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      ArrayVariable  av;
      Read( ip, idx, av );
      return av.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( ip, idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }
    std::cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }


















// =========================================
// FINITE VOLUME INTEGRATION POINT VARIABLES
// =========================================


/** 
    Scalar variable value at sector or facet integration points only.
*/
template<size_t dim,class STOREE>
double64 LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx ) const  
  {
    // offset to first instance of idx variable in the data vector
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offsetData.first < data_.data.size() );
#endif
    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth; // for facet integration points only
    ///                                                                                    ^^^^^^^^
    return data_.data[offsetData.first + sector_ip_offset];
  }



/// Scalar variable at facet or sector integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx, ScalarVariable& sc ) const  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    sc.Flag() = data_.flags[flagOffset + sector_ip_flag_offset];
    sc        = data_.data[offset + sector_ip_offset];
  }




/// Scalar variable at facet or sector integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t sector_or_facet, size_t ip, const csmp::Index& idx, const ScalarVariable& sc )  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    data_.flags[flagOffset + sector_ip_flag_offset] = sc.Flag();
    data_.data[offset + sector_ip_offset]           = sc();
  }





/// Scalar & Array variable flag at facet or sector integration point
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( size_t sector_or_facet, size_t ip, const csmp::Index& idx ) const 
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR || idx.type == ARRAY );
  assert( offsetData.second < data_.flags.size() );
#endif
    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[offsetData.second + sector_ip_offset];
  }




/// Vector, Tensor, FlaggedArray variable flag at integration point
template<size_t dim,class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( size_t sector_or_facet, size_t ip, const csmp::Index& idx, size_t i ) const  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif
    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[ flagOffset + sector_ip_offset + i ];
  }




/// Scalar & Array variable flag at facet or sector integration point (note that the Array has only a single flag)
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( size_t sector_or_facet, size_t ip, const csmp::Index& idx, VARIABLE_FLAG flag ) 
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR  || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif
    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[flagOffset + sector_ip_offset] = flag;
  }



/// Vector, Tensor, FlaggedArray variable flag at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( size_t sector_or_facet, size_t ip, const csmp::Index& idx, size_t i, VARIABLE_FLAG flag ) 
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif
    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[ flagOffset + sector_ip_offset + i ] = flag;
  }




/// Vector variable at integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t sector_or_facet, size_t ip, const csmp::Index& idx, const VectorVariable<dim>& vc )  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i<dim; ++i ) {
         data_.flags[ flagOffset + sector_ip_flag_offset + i ] = vc.Flag(i);
         data_.data[ offset + sector_ip_offset + i ]            = vc[i];
      }
  }




/// Vector variable at facet or sector integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx, VectorVariable<dim>& vc ) const  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this),sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ flagOffset + sector_ip_flag_offset + i ];
         vc(i)      = data_.data[ offset + sector_ip_offset + i ];
      }
  }


/// Tensor variable at sector or facet integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t sector_or_facet, size_t ip, const csmp::Index& idx, const TensorVariable<dim>& ts )
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i=0U; i<dim; ++i ) {
          data_.flags[ flagOffset + sector_ip_flag_offset + i ] = ts.Flag(i);
          for ( size_t j=0U; j<dim; ++j )
            data_.data[ offset + sector_ip_offset + i*dim + j ] = ts(i,j);
      }
 }


/// Tensor variable at facet or sector integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx, TensorVariable<dim>& ts ) const  
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i=0U; i<dim; ++i ) {
          ts.Flag(i) = data_.flags[flagOffset + sector_ip_flag_offset + i];
          for ( size_t j=0U; j<dim; ++j )
            ts(i,j) = data_.data[ offset + sector_ip_offset + i*dim + j ];
      }
  }


/// Array variable at sector or facet integration poin
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t sector_or_facet, size_t ip, const csmp::Index& idx, const ArrayVariable& av )
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

    const size_t arraySize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i < arraySize; ++i )
      data_.data[ offset + sector_ip_offset + i ]     = av[i];
    data_.flags[ flagOffset + sector_ip_flag_offset ] = av.Flag();
  }


/// Array variables at sector or facet integration points
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx, ArrayVariable& av ) const
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

    const size_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i < arraySize; ++i )
      av(i) = data_.data[ offset + sector_ip_offset + i ];
    av.Flag( data_.flags[flagOffset + sector_ip_flag_offset] );
  }



/// FlaggedArray variable at facet or sector integration point
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( size_t sector_or_facet, size_t ip, const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

    const size_t arraySize( idx.dataDepth );

#ifndef NDEBUG  
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i < arraySize; ++i ) {
         data_.data[ offset     + sector_ip_offset + i ]      = av[i];
         data_.flags[flagOffset + sector_ip_flag_offset + i ] = av.Flag(i);
      }
  }



/// FlaggedArray variables
template<size_t dim,class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( size_t sector_or_facet, size_t ip, const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
    const std::pair<size_t, size_t> offsetData = lvsCompileTimeDispatch::containerOffset(static_cast<const STOREE*>(this), idx);
    const size_t offset(offsetData.first);
    const size_t flagOffset(offsetData.second);

    const size_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const size_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const size_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( size_t i(0); i < arraySize; ++i ) {
         av(i) = data_.data[ offset + sector_ip_offset + i ];
         av.Flag( i, data_.flags[flagOffset + sector_ip_flag_offset + i] );
      }
  }


/**
    checks whether any value of the target variable is within the range supplied via the arguments vmin and vmax.
    
     @attention for a vector, tensor, or array this range check is performed on all their elements.
     
     @todo refactor: for vectors, tensors and array variables this method is terribly inefficient as it creates temporaries for the checking.
*/
template<size_t dim,class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( size_t sector_or_facet, size_t ip,
                                                      const csmp::Index& idx,
                                                      double64 vmin, double64 vmax ) const
  {
#ifndef NDEBUG
   lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex(static_cast<const STOREE*>(this), sector_or_facet, ip);
   AssertFiniteVolumeIntegrationPointPlacement(idx);
#endif

    if ( idx.type == SCALAR ) {
         const double64 val = Read( sector_or_facet, ip, idx );
         return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
         VectorVariable<dim>  vc;
         Read( sector_or_facet, ip, idx, vc );
         return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
         TensorVariable<dim>  ts;
         Read( sector_or_facet, ip, idx, ts );
         return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
         ArrayVariable  av;
         Read( sector_or_facet, ip, idx, av );
         return av.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == FLAGGEDARRAY ) {
         FlaggedArrayVariable  fa;
         Read( sector_or_facet, ip, idx, fa );
         return fa.IsWithinRange( vmin, vmax );
      }

    std::cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }

template class LocalVariableStorage<1U,Node<1U> >;
template class LocalVariableStorage<1U,Element<1U> >;
template class LocalVariableStorage<1U,Face<1U> >;
template class LocalVariableStorage<1U,InterFace<1U> >;
template class LocalVariableStorage<1U,ModelSubDomain<1U,Element> >;
template class LocalVariableStorage<1U,ModelSubDomain<1U,Face> >;
template class LocalVariableStorage<1U,ModelSubDomain<1U,InterFace> >;
template class LocalVariableStorage<1U,Model<1U> >;


template class LocalVariableStorage<2U, Node<2U> >;
template class LocalVariableStorage<2U,Element<2U> >;
template class LocalVariableStorage<2U,Face<2U> >;
template class LocalVariableStorage<2U,InterFace<2U> >;
template class LocalVariableStorage<2U,ModelSubDomain<2U,Element> >;
template class LocalVariableStorage<2U,ModelSubDomain<2U,Face> >;
template class LocalVariableStorage<2U,ModelSubDomain<2U,InterFace> >;
template class LocalVariableStorage<2U,Model<2U> >;

template class LocalVariableStorage<3U, Node<3U> >;
template class LocalVariableStorage<3U,Element<3U> >;
template class LocalVariableStorage<3U,Face<3U> >;
template class LocalVariableStorage<3U,InterFace<3U> >;
template class LocalVariableStorage<3U,ModelSubDomain<3U,Element> >;
template class LocalVariableStorage<3U,ModelSubDomain<3U,Face> >;
template class LocalVariableStorage<3U,ModelSubDomain<3U,InterFace> >;
template class LocalVariableStorage<3U,Model<3U> >;

} // end csmp
