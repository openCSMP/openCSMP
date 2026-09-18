// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FRACTURE_H
#define FRACTURE_H

#include "Point.h"
#include "SplitBoundary.h"
#include "Boundary.h"
#include "FiniteElementPolicy.h"
#include "Element.h"

//#include "CohesiveTip.h"
//#include "DisplacementCorrelationTip.h"
//#include "HydraulicFractureTip.h"
//#include "J_IntegralTip.h"
//#include "J_IntegralViscousTip.h"


//#include "FractureTip.h"

template <uint32_t> class FractureTip;
template <uint32_t> class CohesiveTip;
template <uint32_t> class J_IntegralTip;
template <uint32_t> class HydraulicFractureTip;
template <uint32_t> class J_IntegralViscousTip;
template <uint32_t> class DisplacementCorrelationTip;


//Return to Fracture tip once code is fully ported
enum TIP_TYPE { DC_TIP = 0, COHESIVE_TIP = 1, HF_TIP = 2, J_DRY_TIP = 3, J_WET_TIP = 4 , JA_HFM_TIP = 5, J_HFM_TIP};

//Special library for analytical solutions of pressure
//#include "gsl/gsl_sf_hyperg.h"



namespace csmp {

#ifdef CSMP_WITH_SAMG_SOLVER
class SAMG_Settings;
#endif

/**

This class provides a tool to model fractures and accurately handle fracture mechanics.

@author M. Nejati
@date 2016

Extended by
@author E. Pezzulli
@date 2019

@brief The Fracture class is a FractureTip Manager. It is responsible for the creation of fracture tips relevant to a split boundary, and their updating when a FractureTip propagates.
   The Fracture class has 1 Fundamental requirement for creation
      1) A split boundary

   The Fracture class manages two fundamental data types
      1) A collection of TipNodes.
      2) A collection of FractureTips.

      Construction involves: Identifying Tip Nodes of the splitboundary & constructing an associated FractureTip.

      Propagation involves: A modification (without deletion/reconstruction) of the fracture tips and their associated tip nodes.

      @attention From a User Perspective! When advancing to a new time step, the user must set run SetOldCoordinatesToCurrentTip()
        in order to calibrate the old tip coordinate of all FractureTip objects (this ensures propagation velocity is calculated correctly every time step)
*/

// TODO: clean-up partially used macros etc.
template<uint32_t dim>
class Fracture
{
public:
    Fracture();
    Fracture(Model<dim>& model, std::string splitboundary, TIP_TYPE=DC_TIP);

    ~Fracture();


    //Fracture Names
    SplitBoundary<dim>&             SplitBoundaryRef(){return sb_ref_;}
    std::string                     MidRegionName(){ return midregion_->Name(); }
    std::string                     SplitBoundaryName(){return sb_ref_.Name(); }

    //data manipulation
    void                            AveragePropertyToMiddle(Index property, std::vector<InterFace<dim>*>& interfaces);
    void                            InputPropertyValue(const char* prop, ScalarVariable Sval, BOX_BOUNDARY bound_side);

    //Propagation
    //bool                            PropagationAlgorithm(Boundary<dim>& b_ref, double& dt, bool propagate_anyway = false);
    //void                            PropagateTip(FractureTip<dim>* f_tip, Boundary<dim>& b_ref);
    //std::vector<double>           EvaluateTips();
    std::vector<FractureTip<dim>*>  FractureTips() {return fracturetips_;}


    //Updating members
    void                                SetOldTipCoordinatesAsCurrent();
    //void UpdateTipAndBoxBoundaryNodes();
    //void UpdateFractureTips();

    //Access
    //std::vector<Node<dim>*>             TipNodes();
    std::vector<Node<dim>*>             BoxBoundaryNodes(){return boxboundaryNodes_;}         // Nodes of fluid region which intersect with boundary
    std::vector<Element<dim>*>          BoxBoundaryElements();
    std::vector<Element<dim>*>          BoxBoundaryElements( BOX_BOUNDARY side);
    std::map<InterFace<dim>*, size_t>   FindInterFaceOfNode( Node<dim>* n_ptr, INTERFACE_SIDE side, SUBDOMAIN_PART part = COMPLETE  );
    //Maps
    std::map<Point<dim>, Node<dim>*>    NodeMap(INTERFACE_SIDE side);
    std::map<double,   Node<dim>*>      NodeMap(INTERFACE_SIDE side, size_t xyz);
    std::map<Point<dim>, Element<dim>*> ElementMap(INTERFACE_SIDE side);


    //Querying
    double                            FractureLength(INTERFACE_SIDE side) ;
    double                            DistanceFromTip(std::map<Point<dim>,Node<dim>*>&, Node<dim>* N, BOX_BOUNDARY side = RIGHT);
    double                            MaxPropertyValue(const char* prop, INTERFACE_SIDE side = MIDDLE);
    double                            Volume(const char* aperture = "aperture");
    double                            SurfaceIntegral(const char* oper, INTERFACE_SIDE side);

    void                              ConfigureAnalyticalParametersAndSolutions(double constant_flux, double viscosity, double critSIF, double youngs_modulus, double poisson_ratio, bool plane_strain = true);
    double                            AnalyticalAperture(double t, double x = 0.0);
    double                            AnalyticalPressure(double t, double x = 0.0);
    double                            AnalyticalLength(double t);
    double                            DimensionlessToughness();
    double                            DimensionlessViscosity();

    //Updating Middle Mesh
    void                                StoreDisplacementDifferenceAsAperture(const std::string displacement, const std::string aperture, INTERFACE_SIDE side=MIDDLE);


    //Coupled HF Specific configurations
#ifdef CSMP_WITH_SAMG_SOLVER
    void                                SetSolverSettings( SAMG_Settings& settings);
#endif

    //Output
    void                                Out();
    void                                SolOut( INTERFACE_SIDE side = MIDDLE, const char* aperture = "aperture", const char* pressure = "fluid pressure");

    //Testing Methods
    //void                                TestSplitNodeAssignment( bool check_middle_elements );
    //void                                TestKGDMeshPropagation2D(Boundary<dim>& b_path, bool propagate_anyway = false);         //warning, this method splits the mesh!


protected:
    // Initialization
    void                                InitializeQuarterPointElements();
    void                                RestoreQuarterPointElementsToMidPoint();

    //used during construction
    void                                CreateFractureTips();

    //TODO: Deprecate
    void                                SetOldFractureLengthToCurrent(double dt = 0);

    //Subdomains
    Model<dim>*                     model_;
    csmp::SplitBoundary<dim>&       sb_ref_;                 // name of split boundary
    Region<dim>*                    midregion_;


    //Type of propagation algorithm
    TIP_TYPE                        tiptype_;
    //std::vector<Node<dim>*>         tipNodes_;            ///< vector of pointers to the nodes at tips (shared nodes of interface)
    std::vector<FractureTip<dim>*>  fracturetips_;        ///< Container for fracture tips
    std::vector<Node<dim>*>         boxboundaryNodes_;    ///< vector of pointers of middle regions nodes who are at boundary






    //Material keys
    Index           disp_key_;
    Index           E_key_;
    Index           nu_key_;

    bool            configured_, plane_strain_;
    double          Q_, mu_, Kc_, ym_, nu_;
    double          K_dominant, M_dominant, K_small, M_small;

    double          old_fracture_length_;
    double          accumulated_time_;

    //Bools
    bool            two_tips_;


};



}//csmp

#endif

