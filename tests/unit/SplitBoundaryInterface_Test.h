#ifndef SPLIT_BOUNDARY_INTERFACE_TEST_H
#define SPLIT_BOUNDARY_INTERFACE_TEST_H

#include "Test.h"
#include "SplitBoundaryInterface.h"

namespace csmp {

template<uint32_t> class Region;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Model;

// utility functions

template<uint32_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, INTERFACE_SIDE side, double xShift, double yShift, double zShift );

template<uint32_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, double shift );

template<uint32_t dim>
void shiftInterfaceTips( Region<dim>& region, double shift );

template<uint32_t dim>
void shiftRegion( Region<dim>& region, double xShift, double yShift, double zShift );

template<uint32_t dim>
void shiftRegionAboveLine( Region<dim>& region, size_t x_or_y_or_z, double line_coordinate, double shift, double eps );

template<uint32_t dim>
void shiftRegionBelowLine( Region<dim>& region, size_t x_or_y_or_z, double line_coordinate, double shift, double eps );

template<uint32_t dim>
void scaleRegionSymmetricOverZero( Region<dim>& region, double xScale, double yScale, double zScale );

template<uint32_t dim>
void scaleRegion( Region<dim>& region, double xScale, double yScale, double zScale );




/**
    SplitBoundaryInterface_Test unit test tests all functionality related to the creation and destruction of SplitBoundary objects.
    Refer to SplitBoundary_Test for tests of the functionality of the SplitBoundary itself.

     @author Revised testst originally conceived by Roman Manasipov in 2013,
     @author refactored by Junchul Kim (2019).
     @author refactored by Edoardo Pezzulli (2022)
     @date ported to 2016 version by SKM.
 */
 
template<uint32_t dim>
class SplitBoundaryInterface_Test : public Test
    {
    public:
      virtual void run();

      // CURRENTLY RUNNING:
      void Test_splitboundary_from_lower_dim_region(); //E.P Implemented and Tested - Contains RIGOROUS checking of nodes, elements, and interfaces correctly calibrated


      //RUNNING, BUT WITH NO DIAGNOSTICS: Just checks methods dont crash - TODO: Add quantitative checks/tests to each of these
      void Test_splitboundary_between_regions( const std::string& model_name );
      void Test_splitboundary_around_regions( const std::string& model_name );
      void Detect_and_create_splitboundaries( const std::string& model_name );
      void Detect_and_create_splitboundaries_from_constructor( const std::string& model_name );


      // NOT RUNNING/USED: These tests are currently not used/run
      void NodeParents( const csmp::Region<dim>& region, size_t minParentCount = 2 );
      void ElementNodes( const csmp::Region<dim>& region );
      void CheckRemovedLowDimParents( csmp::Boundary<dim>& boundary );
      bool NoNeighborNull( const csmp::InterFace<dim>& interFace );       /// tests whether all InterFace elements of the SplitBoundary have both higher-dimensional neighbors
      void TestSplitNodeAssignment( const csmp::Model<dim>& );      /// tests whether the splitted nodes share parent elements, located on different sides of interfaces
      void TestUnitNormals( csmp::Model<dim>&, const std::string& test_name  );

      //Visualisation functions
      void VisualiseSplitBoundaries( csmp::Model<dim>&, const std::string&test_name  );
      void PullApartSplitboundaries( Model<dim>& model, std::vector<std::string>& fractures, double displacement );
      void PullApartSplitboundaries( Model<dim>& model, double displacement );
      void LoadModel( const std::string& model_name );
      void LoadModel( const std::string& model_name, std::vector<std::string>& regions );
      void LoadContiguousModel( const std::string& model_name, std::vector<std::string>& fractures );
      void EstablishContiguousRegionsList( Model<dim>& model,
                                           const std::set<std::string>& fractures_basic_set,
                                           std::set<std::string>& fractures );

      /// reading fault modeling input data
      void InputFromFile( const char* file_name,
                          std::set<std::string>& fractures_basic_set );

      void InputFromFile( const char* file_name,
                          std::vector<std::string>& fractures );

      /// writing fault modeling input data into text file
      void OutputToFile( const char* file_name,
                         const std::vector<std::string>& fractures );
                         
      void OutputToFile( const char* file_name,
                         const std::set<std::string>& fractures );



protected:

      // CreateFromRegion Test Added by E.P
      bool Test_NodeCorrespondance_2D(const char* mesh_file );
      bool Test_NodeAndElementsCorrespondance_3D(const char* mesh_file );
      bool Test_NodeAndElementsCorrespondance_3D_X_Intersection( const char* mesh_file);
      bool Test_NodeAndElementsCorrespondance_3D_X_Intersection_Reverse( const char* mesh_file);


      private:
        static const bool verbose_ = true;
  };

  } // csmp

#endif /* SPLIT_BOUNDARY_INTERFACE_TEST_H */
