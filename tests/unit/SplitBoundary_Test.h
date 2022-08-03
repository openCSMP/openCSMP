#ifndef SPLIT_BOUNDARY_TEST_H
#define SPLIT_BOUNDARY_TEST_H

#include "Test.h"
#include "SplitBoundary.h"
#include "SplitBoundaryInterface_Test.h"

namespace csmp {

/**
    SplitBoundary_Test of the functionality of individual SpitBoundary objects
    as opposed to their creation (see SplitBoundaryInterface_Test) or that of the InterFace
    objects that the SplitBoundary objects consist of (see InterFace_Test).
    
    @attention this also performs test on SplitBoundary objects created in ANSYS
    
    @section Tests Performed Here (2D and 3D)
    
    - do all interfaces have higher dimensional neighbors
    - do correspoding nodes coincide
    - is there only a single node along the perimeter of a split boundary
    - are all of the above operations still valid for a SplitBoundary read from binary file?

     @author SKM
     @date 26/1/2020
 */
  class SplitBoundary_Test : public Test
    {
    public:
      SplitBoundary_Test() {}
    
      virtual void run();
      
      // TESTING of public methods - Top ones are implemented
      // ----------------------------------------------------------------------------
      // SplitBoundary creation gets tested by SplitBoundaryInterface_Test
       
      /// initialises 2022 SKUA model with properties from binary version of VSet
      void Test_InitialiseSKUA_Model( const std::string& file_name );
        
      /// tests that all properties are properly added to nodes on split boundary
      bool Test_InputNodePropertyValue(const char* mesh_file);

      /// integrates the property over the boundary line or surface
      bool Test_Area_and_SurfaceIntegral(const char* mesh_file);

      // TODO - ALl tests below are still to be done

      /// returns location of split boundary relative to adjacent region
      bool Test_FacingDirection();

      /// outputs length of perimeter curve of a 3D split boundary; no meaning in 1 or 2D models
      bool Test_PerimeterCalculation();
      
        /// integrates the property over the boundary line or surface
      bool Test_SurfaceIntegral();

      /// do all interfaces have higher dimensional neighbors
      bool Test_HigherDimensionalNeighbors();
    
      /// do correspoding nodes coincide
      bool Test_NodeCorrespondance();
      
      bool Test_Manifolds();
      
      /// is there only a single node along the perimeter of a split boundary
      bool Test_SingleParameterNode();
      
      /// are all of the above operations still valid for a SplitBoundary read from binary file?
      bool Test_IntegrityOfSplitBoundaryFromBinaryFile();
      
      // TEST SETUP
      // ----------------------------------------------------------------------------

      template<uint32_t dim>
      void LoadModel( const std::string& model_name );

      template<uint32_t dim>
      void LoadModel( const std::string& model_name, std::vector<std::string>& regions );

      template<uint32_t dim>
      void LoadContiguousModel( const std::string& model_name, std::vector<std::string>& fractures );

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
                         
      template<uint32_t dim>
      void EstablishContiguousRegionsList( Model<dim>& model, 
                                           const std::set<std::string>& interface_basic_set,
                                           std::set<std::string>& interface_sets );
      private:
        static const bool verbose_ = true;
  };

  } // csmp

#endif
