#ifndef MODEL_COMPARATOR_H
#define MODEL_COMPARATOR_H

#include "Region.h"
#include "Point.h"
#include "PropertyAtPointVisitor.h"

namespace csmp {
  
  struct Index;
  class ScalarVariable;
  template<uint32_t> class Model;
  template<uint32_t> class VSet;

  /// allows for comparison between two csmp models
  template<uint32_t dim>
  class ModelComparator
  {
  public:

    /////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for identically equal meshes with identical numberiings of nodes//
    /////////////////////////////////////////////////////////////////////////////////////////

    double CompareModels( const Model<dim>& model1, const Model<dim>& model2,
                            const char* property1, const char* property2, const char* region = "Model" ) const;


    double CompareVSets( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                           const char* variablesFile1, const char* variablesFile2,
                           bool isoparametric = false, const char* region = "Model" ) const;
    
    double CompareVSets( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                           const char* variablesFile1, const char* variablesFile2,
                           bool isoparametric = false, const char* region = "Model" ) const;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for identically equal meshes but might be different numbering of nodes //
    ////////////////////////////////////////////////////////////////////////////////////////////////

    double CompareModelsRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                           const char* property1, const char* property2, const char* region = "Model" );

    double CompareVSetsRenumberedNodes( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                                          const char* variablesFile1, const char* variablesFile2,
                                          bool isoparametric = false, const char* region = "Model" );
    
    double CompareVSetsRenumberedNodes( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                           const char* variablesFile1, const char* variablesFile2,
                           bool isoparametric = false, const char* region = "Model" );

    /////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for meshes that can be in general different from each other     //
    /////////////////////////////////////////////////////////////////////////////////////////

    double CompareModelsAtPoints( Model<dim>& model1, Model<dim>& model2,
                            const char* property1, const char* property2, const char* region = "Model" );

    double CompareVSetsAtPoints( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                           const char* variablesFile1, const char* variablesFile2,
                           bool isoparametric = false, const char* region = "Model" );
    
    double CompareVSetsAtPoints( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                           const char* variablesFile1, const char* variablesFile2,
                           bool isoparametric = false, const char* region = "Model" );

  private:

    double CompareRegionScalarVariable( const Model<dim>& model1, const Model<dim>& model2,
                                          const Index model1Key, const Index model2Key, const char* region = "Model" ) const;

    double CompareScalarDequesL2( const std::deque<ScalarVariable>& deque1,
                                    const std::deque<ScalarVariable>& deque2 ) const;

    void ReadRegionElementScalarVariable( const Region<dim>& region,
                                          const Index propKey,
                                          std::deque<ScalarVariable>& scalarDeque ) const;

    void ReadRegionNodalScalarVariable( const Region<dim>& region,
                                        const Index propKey,
                                        std::deque<ScalarVariable>& scalarDeque ) const;

    double CompareRegionScalarVariableRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                          const char* property1, const char* property2, const char* region = "Model" );

    void ReadRegionNodalScalarVariableAndNodeCoordinates( const Region<dim>& region,
                                        const Index propKey,
                                        std::map<Point<dim>, ScalarVariable >& points_and_values);

    double CompareRegionScalarVariableAtPoints( Model<dim>& model1, Model<dim>& model2,
                                          const char* property1, const char* property2, const char* region = "Model" );

    void ReadRegionNodalScalarVariableAndNodeCoordinates( const Region<dim>& region,
                                        const Index propKey,
                                        std::vector<ScalarVariable> & values,
                                        std::map<size_t, std::vector<double> >& points);
  };

  /**
  @class ModelComparator
  @author P. Lang
  @date Aug 2011

  Compares Models, Vsets. For usage example, please refer to ModelComparator_Test
  for the time being.

  @todo (1) seperate comparisons of meshes from properties discretised on them
  @todo (3) Different error norms
  @todo (3) Different properties of same model
  @todo (3) Cor diff var types
  @todo (3) Cefactor(outsource) functions once diff var types implemented
  @todo (3) Cefactor static, no sense anymore - DONE(or not so much???)
  @todo (3) Cefaults for var file param
  @todo (3) Compare list of properties
  */

} // csmp

#endif
