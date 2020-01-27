#ifndef ANSYS_MODEL_QUADRATIC_H
#define ANSYS_MODEL_QUADRATIC_H

#include "ANSYS_Interface.h"
#include "ANSYS_ElementSpecifications.h"

namespace csmp {

  class FiniteElementManager;
  class LocalVariables;
  template<size_t> class Point;
  template<size_t> class Node;
  template<size_t> class Element;
  

  /** ANSYS interface for quadratic isoparametric elements
 
  Motivation behind creating this kind of redundant class is that ANSYS_Model2D is not able to read
  quadratic meshes from Icem input files and I gave up debugging. It performs nicely for smaller models,
  for larger models it's slower than ANSYS_Model2D or ANSYS_Model3D, mostly due to the limitation to ascii
  outputs. Instead of a '...-regions.txt' file it requires a '...-subdomains.txt' file which is supposed
  to look as such

  @code
  regions
  MATRIX
  boundaries
  LEFT
  RIGHT
  split
  FRACTURE
  eof
  @endcode

  where MATRIX, LEFT, RIGHT and FRACTURE are families listed in the '...asc' file. Here, 'MATRIX' becomes
  a Region, 'LEFT' and 'RIGHT' a Boundary and 'FRACTURE' a SplitBoundary.
  
  @author P. Lang
  @date 2012

  @attention This reads from ascii '.dat' files only and is hence slow for larger models.

  @todo  Relies on counter clockwise ordering in ANSYS
  @todo  Allows for single element type regions only
  @todo  Mesh either all linear, all quadratic or all cubic
  @todo  ASCII .dat files only
  @todo  'eof' tag
  */
  template<size_t dim>
  class ANSYS_Model_Quadratic : public Model<dim>
    {
      enum INTERPOLATION_ORDER{ LINEAR=1, QUADRATIC=2, CUBIC=3 };

    public:

      /// input from ANSYS *.asc, *.dat, *-variable.txt and *-regions.txt files.
      ANSYS_Model_Quadratic( const std::string& mesh_file_set,
                             const std::string& regions_file_prefix,
                             const std::string& variable_file );

      virtual ~ANSYS_Model_Quadratic();

    private:

      void Initialize();
      void AnnounceStart() const;
      void AnnounceEnd() const;
      bool ReadInData();
      void ReadInFamilies();
      void ReadInSubDomains();
      bool IsNot( std::ifstream& file,
                  std::string& line,
                  const std::string& isNot ) const;
      int  ExtractEntityCount( const std::string& line );
      void EraseNodes();
      void EraseElements();
      void InitializeLineElementNeighbors( size_t regionIndex );
      void ExtractRegionInfo( std::ifstream& file, std::string& line );
      bool RegionInfoExists( const std::string& regionName );
      void NeighborConnectivity();
      void CreateRegions();
      void CreateBoundaries();
      void CreateSplitBoundaries();
      void CreateRegion( const std::string& regionName,
                         typename std::vector<Element<dim>*>::const_iterator elementsBegin,
                         typename std::vector<Element<dim>*>::const_iterator elementsEnd );
      void DeleteRedundantElements();
      size_t InterpolationOrder( const std::vector<CSMP_FEM_TYPE>& femTypes ) const;
      void CheckSubDomainRequest() const;
      void AnnounceSubDomains() const;
      void AnnounceSubDomains( const std::vector<std::string>& subdomains,
                               const std::string& subdomain ) const;
      std::string ParseInterpolationOrder() const;

    private:

      std::string                       meshFileName_;
      std::string                       topologyFileName_;
      size_t                            interpolationOrder_;

      std::vector<Node<dim>*>           nodes_;
      std::vector<Element<dim>*>        elements_;

      std::set<size_t>                  elementsInUse_;
      std::set<std::string>             subdomains_;
      std::set<std::string>             families_;
      std::vector<std::pair<std::string,std::vector<size_t> > > regionElements_;

      std::vector<std::string>          regions_;
      std::vector<std::string>          boundaries_;
      std::vector<std::string>          splitBoundaries_;

    };

} // csmp

#endif         
