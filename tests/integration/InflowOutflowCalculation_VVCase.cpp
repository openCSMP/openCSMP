#include "InflowOutflowCalculation_VVCase.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "LinearSolver.h"
#include "PropertyHandle.h"
#include "PDE_Integrator.h"
#include "ScalarVariable.h"
#include "SteadyStateDiffusor.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "VelocityAndVolumeFlux.h"
#include "Timer.hpp"
#include "PL_Utilities.h"


using namespace std;

namespace csmp
{

template<uint32_t dim>
InflowOutflowCalculation_VVCase<dim>::InflowOutflowCalculation_VVCase()
  {
  }

template<uint32_t dim>
InflowOutflowCalculation_VVCase<dim>::InflowOutflowCalculation_VVCase( const char* prefix )
  {
    name_="InflowOutflowCalculation_VVCase";
    prefix_=prefix;
  }

template<uint32_t dim>
InflowOutflowCalculation_VVCase<dim>::~InflowOutflowCalculation_VVCase()
  {
  }

/*================================================================
  Calculation of inflow vs outflow rates between two boundaries.
  ================================================================
  Element type(s) : Linear Tests
  Target Mesh set : UnitCubeFine
  Test            : Steady State Pressure (Incompressible Single Phase Porous Media Flow)
                    Uses the NodeCenteredFiniteVolumeTransport class to calculate inflow vs outflow
                    of variable "concentration"

  =================================================================
*/
template<uint32_t dim>
void InflowOutflowCalculation_VVCase<dim>::run()
  {
    string input_file_name(prefix_);

    const ScalarVariable viscosityWater( PLAIN, 1.0E-03 );
    const ScalarVariable permeability( PLAIN, 1.0E-12 );
    const ScalarVariable porosity( PLAIN, 1.0E-01 );
    const ScalarVariable pressureLeft( DIRICH, 1.0E+05 );
    const ScalarVariable pressureRight( DIRICH, 0.0 );
    ScalarVariable concentration( DIRICH, 0.8 );

    // ----------------------------------
    // GEOMETRY SECTION
    // establishing model & output facility
    cout <<"Building Model..."<<endl;

    Model<dim>* model;
    switch(dim){
    case 1:
        cout<<"InflowOutflowCalculation_VVCase:: No calculation for 1D models is available yet"<<endl;
        exit(1);
    case 2:
        model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D(input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str()));
        break;
    case 3:
        model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D(input_file_name.data(),this->getName().c_str(),(this->getName()+".txt").c_str()));
        break;
    }

    //
    cout <<"Finished reading mesh..."<<endl;

    model->Region( "Model" ).InputPropertyValue( "fluid pressure",      makeScalar(PLAIN, 0.0) );
    model->Region( "Model" ).InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0.0 ) );
    model->Region( "Model" ).InputPropertyValue( "porosity",            porosity );
    model->Region( "Model" ).InputPropertyValue( "conductivity",        permeability/viscosityWater);

    //fluid pressure left right const dirich
    model->InputBoundaryValue( LEFT,  "fluid pressure", pressureLeft );
    model->InputBoundaryValue( RIGHT, "fluid pressure", pressureRight );

    model->InputBoundaryValue( LEFT,  "concentration", concentration );

    SteadyStateDiffusor<dim, Region> ssds( *model, "conductivity", "fluid pressure", "fluid volume source");
    VelocityAndVolumeFlux<dim,Element<dim> > veloandvflux( *model,
                                             "conductivity",               // conductivity (abs. perm./visc. for single phase)
                                             "porosity",                   // porosity
                                             "fluid pressure",             // fluid pressure
                                             false,
                                             "velocity");
    Boundary<dim>& modelBoundary( model->Boundary( "LEFT" ) );

    ssds.AddPostProcess(&veloandvflux);
    Timer timer;
    timer.Start();
    ssds.ComputeSteadyState( model->Region("Model") );
    cout << "\n\nSAMG took " << timer.Stop() << " sec\n";

    singlePhaseVelocity( *model, "Model" );

    model->InputPropertyValue( "concentration",                         concentration );
    model->InputPropertyValue( "nodal fluid volume source",             makeScalar(PLAIN,0.) );

    NodeCenteredFiniteVolumeTransport<dim>  advector( "Model", *model,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source", false, false );


    //Assuming incompressible flow, a constant concentration value everywhere (including LEFT boundary)
    //The total flux of concentration is estimated to be:
    const double q = concentration()*modelBoundary.Area()*permeability() * pressureLeft() / viscosityWater();
    //const double qNum( model->Region("Model").Average("volume flux") );

    /// @attention Due to comments in NCFVT code, I have changed this test to use Stephan's functionality.
    /// ComputeBoundaryFlux no longer exists. (Julian 12-10-2014)
    double inflow(0.);
    double outflow(0.);
    advector.BoundaryFluxes(inflow,             // inflow  ( magnitude )
                            outflow,            // outflow ( magnitude )
                            false,              // non-box shaped model
                            true,               // use advected variable
                            "concentration",    // advected variable
                            "fluid pressure",   // variable to determine no-flow conditions
                            ANY);

    cout<<endl;
    cout<<"Assuming incompressible flow, the estimated Concentration flux"<<endl;
    cout<<"based on input variables through LEFT boundary is: "<<q<<endl;
    cout<<"Obtained values of inflow: "<<fabs(inflow)<<" and outflow: "<<fabs(outflow)<<endl;
    cout<<endl;

    _equal( fabs(inflow),q, 1E-7 );
    _equal( fabs(inflow),fabs(outflow), 1E-7 );
    _equal( fabs(inflow)-fabs(outflow),0.0, 1E-7 );


    //    inflow=advector.ComputeBoundaryFlux("LEFT");
    //    outflow=advector.ComputeBoundaryFlux("RIGHT");
    //    _equal( fabs(inflow),fabs(outflow), 1E-7 );
    //    _equal( inflow+outflow,0.0, 1E-7 );

    //now with old functionality
    //    const double QNum = advector.ModelInflow();
    //    // unit area
    //    const double Q = q*concentration();
    //    cout << endl << endl << "Qnum: " << QNum << endl << "Hack: " << boundaryInflow( *model, model->Boundary("LEFT") )<< endl << "Analytic: " << Q << endl;
    //    //_equal( Q, QNum, 1.0E-07 );

  }


template class InflowOutflowCalculation_VVCase<2U>;
template class InflowOutflowCalculation_VVCase<3U>;

} //csmp
