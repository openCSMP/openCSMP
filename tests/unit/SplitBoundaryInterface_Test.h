#ifndef SPLIT_BOUNDARY_INTERFACE_TEST_H
#define SPLIT_BOUNDARY_INTERFACE_TEST_H

#include "CSMP_definitions.h"
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
     
    @todo let's try to use a smaller test model that does not come from ANSYS.
 */
class SplitBoundaryInterface_Test : public Test
    {
   public:
      virtual void run();

      // CURRENTLY RUNNING:
      void Test_CreateSplitBoundaryFromLowerDimRegion();
      // calls
      void Test_ConversionOfNormalFault( const std::string& model_name="fault_boundary_test" );
      void Test_ConversionOfNormalFaultInSplitModel( const std::string& model_name="fault_boundary_test" );

      //RUNNING, BUT WITH NO DIAGNOSTICS: Just checks methods dont crash
      template<uint32_t dim>
      void Test_CreateSplitBoundaryBetween( const std::string& model_name );
      
      // DetectAndCreateSplitBoundaries() - is tested in ANSYS_SplitBoundaryMatch_Test
      
      // TODO: test the following
      //   size_t FormSplitBoundariesFrom( const ModelTopology& );
      //   size_t SeparateUniqueRegionsBySplitBoundaries();
      //   std::pair<std::string,bool>  InsertRegionIntoSplitBoundary( const char* split_boundary, int32_t material_id_for_new_elements );
      //   InsertLowerDimensionalRegionsIntoSplitBoundaries( int32_t material_id_for_new_elements );
      //   tested inside of CreateSplitBoundaryBetween: InsertLowerDimensionalRegionsIntoSplitBoundaries()
      //   bool  CreateSplitBoundaryFromInterfaces( const char* name_of_new_splitboundary ); -> uses AddSplitBoundary()
      //   void RemoveSplitBoundary( const char* split_boundary, bool erase_interfaces );
      //   std::string MergeSplitBoundaries( const char* new_sb_name, const std::set<std::string>& splitboundaries );
      //   void CreateNonUniqueSplitBoundaryGroup( const std::set<std::string>& input_split_boundaries, const char* exact_new_split_boundary_name );
      

      //Visualisation functions
      template<uint32_t dim>
      void VisualiseSplitBoundaries( csmp::Model<dim>&, const std::string&test_name );

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
