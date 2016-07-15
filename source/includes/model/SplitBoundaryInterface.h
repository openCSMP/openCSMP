#ifndef SPLIT_BOUNDARY_INTERFACE_H
#define SPLIT_BOUNDARY_INTERFACE_H

#include "SplitBoundary.h"
#include "Region.h"
#include "ErrorHandler.h"

#include <map>

namespace csmp
{

template<size_t dim, template<size_t> class SPLITBOUNDARY_COMPLEX>
class SplitBoundaryInterface
  {
  public:

    SplitBoundaryInterface();
    virtual ~SplitBoundaryInterface();

    csmp::SplitBoundary<dim>&             SplitBoundary( const std::string& spbName );
    const csmp::SplitBoundary<dim>&       SplitBoundary( const std::string& spbName ) const;
    bool                                  ContainsSplitBoundary( const std::string& bname ) const;


    // -----------------------------------------------
    // Binary input/output
    // -----------------------------------------------

    bool OutputSplitBoundariesToBinary( const char* fileName ) const;
    bool InputSplitBoundariesFromBinary( const char* fileName );


    // -------------------------------------------------------------------------
    // SplitBoundaries access and manipulations with the list of SplitBoundaries
    // -------------------------------------------------------------------------

    typedef typename std::map<std::string,csmp::SplitBoundary<dim> >::iterator         splitBoundaryIterator;
    typedef typename std::map<std::string,csmp::SplitBoundary<dim> >::const_iterator   splitBoundaryConstIterator;

    splitBoundaryIterator                 SplitBoundariesBegin();
    splitBoundaryIterator                 SplitBoundariesEnd();
    splitBoundaryConstIterator            SplitBoundariesBegin() const;
    splitBoundaryConstIterator            SplitBoundariesEnd() const;

    size_t                                SplitBoundaries() const;

    /// Removes splitboundary
    void RemoveSplitBoundary( csmp::SplitBoundary<dim>& splitboundary, bool deleteElements = false );

    // -----------------------------------------------
    // SplitBoundaries creation
    // -----------------------------------------------

    /// Creation of SplitBoundary around region
    bool InsertSplitBoundary( const std::string& group, bool deleteRegionAndItsElements = false );

    /// Creation of SplitBoundary between regions
    bool InsertSplitBoundary( const std::string& group1, const std::string& group2, bool createRegionBetween = false );

    // -----------------------------------------------
    // Manipulations with already existed SplitBoundaries
    // -----------------------------------------------
    bool IsEligibleSplitBoundaryRegionName( const std::string& regionName ) const;

  protected:

    std::map<std::string,csmp::SplitBoundary<dim> >   interFaceSplitBoundaryMap_;

  private:

    SplitBoundaryInterface( const SplitBoundaryInterface& bd );

  };

} // csmp

#endif
