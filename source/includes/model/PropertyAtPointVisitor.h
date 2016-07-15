#ifndef PROPERTY_AT_POINT_VISITOR_H
#define PROPERTY_AT_POINT_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "InterFace.h"
#include "Region.h"
#include "FiniteElementTraits.h"

namespace csmp {

template<size_t dim>
class PropertyAtPointVisitor : public Visitor<dim> {

public:
    /// attempts to retrieve the value of the property at the point in the mesh
    PropertyAtPointVisitor( const Model<dim>&,
                            const std::map<size_t,std::vector<double64> >& points_to_search,
                            const char* node_property_name );
  
    ~PropertyAtPointVisitor();


    virtual void Visit(Element<dim>*);
    virtual void Visit(Model<dim>*);

    std::map<size_t,size_t>& TargetFound();
  
    /// returns a pointer to the element that contains the point of interest
    Element<dim>*  ElementThatContains(size_t point) const;

    void PropertyValueAt( size_t point, ScalarVariable& outProperty );
    void PropertyValueAt( size_t point, VectorVariable<dim>& outProperty );
    void PropertyValueAt( size_t point, TensorVariable<dim>& outProperty );
    void PropertyValueAt( size_t point, ArrayVariable& outProperty );
    void PropertyValueAt( size_t point, FlaggedArrayVariable& outProperty );

    void SetPropertyKeys( const char* propS_ );
    void SetDebugOn();
    void SetContinueSearchOff();
    void SetBruteForceOff();
    void SetBruteForceOn();

    bool CheckResults();

private:
    const PropertyDatabase<dim>&                pref_;
    csmp::Index                                 prop_idx_;
    csmp::Index                                 bc_idx_;
    size_t                                      NumberOfElementNodes_;
    size_t                                      maxElementsInTheMesh_;
    size_t                                      numberOfPoints_;
    size_t                                      numberOfPointsFound_;

    std::map<size_t,size_t>                     elements_ids_found_;
    std::map<size_t, Element<dim>* >            elements_found_;
    std::vector<double64>                       props_;
    std::map<size_t, ScalarVariable >           propS_;
    std::map<size_t, VectorVariable<dim> >      propV_;
    std::map<size_t, TensorVariable<dim> >      propT_;
    std::map<size_t, ArrayVariable >            propA_;
    std::map<size_t, FlaggedArrayVariable >     propFA_;
    std::vector<ScalarVariable >                NPS_;       // scalar property at the nodes
    std::vector<VectorVariable<dim> >           NPV_;       // vector property at the nodes
    std::vector<TensorVariable<dim> >           NPT_;       // tensor property at the nodes
    std::vector<ArrayVariable >                 NPA_;       // array property at the nodes
    std::vector<FlaggedArrayVariable >          NPFA_;      // flagged array property at the nodes


    std::vector<double64>                       NI_;        // test-function vector
    std::map <size_t, std::vector<double64> >   xyz_;       // cloud of points


    bool                                        debug_;
    bool                                        continueSearch_;
    bool                                        bruteforce_;

    Element<dim>*                               TargetElement_;
    Element<dim>*                               CurrentElement_;
    size_t                                      startPoint_;
    size_t                                      maxIterations_;
    size_t                                      currIteration_;
    std::map <size_t, double64>                 rangingOfCloud_;
    double64                                    precision_;

    bool isCloseToBarycenter        ( const std::vector<double64> currXyz, Element<dim>* e, VectorVariable<dim> bc );
    bool FindPoint_BruteForceSearch ( const std::vector<double64> currXyz, Element<dim>* e, VectorVariable<dim> bc );
    bool FindPoint_NeighborSearch   ( const std::vector<double64> currXyz, Element<dim>* e, VectorVariable<dim> bc );

};


} // end csmp

#endif
