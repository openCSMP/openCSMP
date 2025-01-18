#ifndef NODE_MANIFOLD_MANAGER_TEST_H
#define NODE_MANIFOLD_MANAGER_TEST_H

#include "Test.h"
#include "NodeManifoldManager.h"
#include "Model.h"


namespace csmp {

  /**
    NodeManifoldManager_Test unit test tests all functionality related to the creation and destruction of NodeManifold objects.
  */

  class NodeManifoldManager_Test : public Test
  {
    public:
      virtual void run();

      //template<uint32_t dim>
      //void Create_splitboundary_between_regions( const std::string& model_name );

      template<uint32_t dim>
      void Test_nodemanifolds_created_from_splitboundaries_between_regions( const std::string& model_name );

      template<uint32_t dim>
      void Test_nodemanifolds_created_from_splitboundaries_around_regions( const std::string& model_name );      
      
      template<uint32_t dim>
      void Test_created_manifolds( Model<dim>& modelIN );

      template<uint32_t dim>
      void Create_splitboundary_between_regions( Model<dim>& modelIN );    

      template<uint32_t dim>
      void Create_splitboundary_around_regions( Model<dim>& model,
                                                std::vector<std::string>& interfaces );

      void InputFromFile( const char* file_name,
                          std::set<std::string>& fractures_basic_set );

      template<uint32_t dim>
      void EstablishContiguousRegionsList( Model<dim>& model,
                                           const std::set<std::string>& fractures_basic_set,
                                           std::set<std::string>& fractures );    

      void OutputToFile( const char* file_name,
                         const std::vector<std::string>& fractures );               
      

      private:
        static const bool verbose_ = true;
  };

  } // csmp

#endif /* NODE_MANIFOLD_MANAGER_TEST_H */
