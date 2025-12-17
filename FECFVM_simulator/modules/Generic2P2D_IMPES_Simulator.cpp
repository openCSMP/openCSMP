#include "Generic2P2D_IMPES_Simulator.h"
#include "IMPES_Setup.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "Model.h"
#include "Region.h"
#include "Node.h"
#include "Element.h"
#include "AuxillaryFunctions.h"
#include "PressureSolver.h"
#include "PropertyHandle.h"
#include "CSMP_highLevelUtilities.h"
#include "Standard_IO_Handler.h"
#include "FiniteElement.h"
#include "TransportModel.h"
#include "ModelTime.h"
#include "SignalHandler.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "samg.h"
#endif

// two phase flow models
#include "TwoPhaseModel.h"
#include "BrooksCorey.h"
#include "VanGenuchten.h"
#include "ExperimentalRT.h"
#include "HeterogeneityAndRateAwareModel.h"

// Transport classes
#include "TwoPhaseElementBasedTransport.h"
#include "ActiveElementTwoPhaseTransport.h"
#include "ActiveFaceTwoPhaseTransport.h"



struct pointcomp {
    bool operator() (csmp::Point<2U> lhs, csmp::Point<2U> rhs) const
    {
        if (lhs[0] != rhs[0])
            return lhs[0] < rhs[0];
        else
            return lhs[1] < rhs[1];
        
    }
};


using namespace std;

namespace csmp {

Generic2P2D_IMPES_Simulator::Generic2P2D_IMPES_Simulator( const char* model_name,
                                                          bool  ansys_model_true_or_csmp_binary_false,
                                                          IMPES_Setup<2>& setup)
 : modelName_     (model_name),
   setup_         (setup),
   reservoirModel_(NULL),
   relpermModel_  (NULL),
   pressureSolver_(NULL),
   advector_      (NULL),
   transportModel_(NULL),
   propKeys_      (NULL),
   vtuOutput_     (NULL)
 {
	std::string variables_file("IMPES-variables.txt");

    if( ansys_model_true_or_csmp_binary_false )
    {
        // input from ANSYS model
        const bool binary_file      ( true  );
        const bool use_regions_file ( true  );
        //const bool create_boundaries( true  );
        reservoirModel_ = new ANSYS_Model2D( model_name,
                                             variables_file.c_str(),
                                             binary_file,
                                             use_regions_file );
    }
    else
    {
        // binary input option
        reservoirModel_ = new Model<2U>( string(model_name) );
    }

    propKeys_  = new IMPES_SimulatorKeys<2U>( reservoirModel_->Database() );
    vtuOutput_ = new VTU_Interface<2U>( *reservoirModel_ );
 }


Generic2P2D_IMPES_Simulator::~Generic2P2D_IMPES_Simulator()
 {
    delete reservoirModel_;
    delete relpermModel_;
    delete pressureSolver_;
    delete advector_;
    delete transportModel_;
    delete propKeys_;
    delete vtuOutput_;
 }



void Generic2P2D_IMPES_Simulator::ConfigureModel()
{
    printModelDimensions( *reservoirModel_, false);

    setup_.ReadSettingsFromFile(reservoirModel_->Database());
    setup_.ReportSetting();

    InputDataManager<2U>  model_configuration;
    model_configuration.ConfigureFromFile(*reservoirModel_,
                                           modelName_.c_str(),
                                           false,
                                           true,   // 2) default prop.values
                                           true,   // 3) group prop.values
                                           true,   // 4) essential conditions
                                           true,   // 5) essential flags
                                           false,  // 6) arbitrary model shape boundary conditions
                                           runSettings_ );
    //WritePermAndVolmodValues();
    //ReadPermVolmod();


    // initializing uninitialized values for volume modifier
    Region<2U>&  gref(reservoirModel_->Region("Model"));
    ScalarVariable one_scalar(PLAIN, 1.);
    ScalarVariable temp;
    
    for ( auto it = gref.CellsBegin(); it != gref.CellsEnd(); ++it )
    {
        temp = (*it)->Read(propKeys_->vm_key);
        if (::isnan(temp()))
        {
            (*it)->Store(propKeys_->vm_key, one_scalar);
        }
        else if (temp() < 0.)
        {
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::ConfigureModel()",
                  "The 'volume modifier' variable must be greater than or equal to zero." );
        }
    
    }
    
    // creating transport model

    transportModel_ = new TransportModel<2U> (*reservoirModel_,
                                              "Model",
                                              "volume modifier",
                                              1U,
                                              1U,
                                              1U + (size_t) setup_.WithCapillaryForces(),
                                              2U);

/*
// for stephan
//  Region<2U,Element>&  gref(reservoirModel_->Region("Model"));
  csmp::Index sw_idx((reservoirModel_->Database().StorageKey("saturation water")));
  csmp::Index so_idx((reservoirModel_->Database().StorageKey("saturation oil")));


  for (vector<Element<2U>*>::iterator
       it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it)
  {
  
      if ((*it)->BaryCenter()[1U] < 400)
      {
          (*it)->Store(propKeys_->sw_key, ScalarVariable(PLAIN, 1.));
          (*it)->Store(propKeys_->so_key, ScalarVariable(PLAIN, 0.));
      }
  
  }
  
  

  Region<2U,Element>& well_ref(reservoirModel_->Region("WELL"));
  double perm((*well_ref.CellsBegin())->Read(propKeys_->k_key));
  double vol_mod((*well_ref.CellsBegin())->Read(propKeys_->vm_key));

  SeparateContinuousLineElementsToDifferentRegionAndAssignProperty();
  
  well_ref.InputPropertyValue("permeability", ScalarVariable(PLAIN, perm));
  well_ref.InputPropertyValue("volume modifier", ScalarVariable(PLAIN, vol_mod));

 

  delete transportModel_;
  
  transportModel_ = 0;
  
  transportModel_ = new TransportModel<2U> (reservoirModel_, "Model", "volume modifier", 1U, 1U, 1U + (size_t) setup_.WithCapillaryForces(), 2U);
// end for stephan
*/
/*    
    if (transportModel_->CVEs() > 30000)
    {
        cout << "\nThis version of simulator only works upto 30000 elements!\n";
        cin.get();
        exit(1);
        
    }
*/
    
    // setting two phase flow model
    switch (setup_.TwoPhaseFlowModel())
    {
        case BROOKSCOREY:

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING

            relpermModel_ = new BrooksCorey<2U>(reservoirModel_->Database(),
                                                "permeability",
                                                "viscosity oil", "viscosity water",
                                                "density oil",	 "density water",
                                                "brooks corey parameter",
                                                "entry pressure","pcmax",
                                                "saturation water",
                                                "residual saturation non-wetting phase", "residual saturation wetting phase",
                                                false);
#else
            relpermModel_ = new BrooksCorey<2U>(reservoirModel_->Database(), "brooks corey parameter", "entry pressure", false);
#endif
            break;

        case VANGENUCHTEN:
            relpermModel_ = new VanGenuchten<2U>(reservoirModel_->Database(), "van genuchten parameter", "van genuchten alpha", false, false);
            break;

        case EXPERIMENTAL:
            relpermModel_ = new ExperimentalRT<2U>(reservoirModel_->Database(), "rock type.txt", "rock type", false);
            break;

        case HETEROGENEITY_AWARE:
             relpermModel_ = new HeterogeneityAndRateAwareModel<2U>( reservoirModel_->Database(), false );
            break;

       default:

            std::cout << "Warning! Generic2P2D_IMPES_Simulator::ConfigureModel()" <<
                         " Two phase flow model was not recognized. It is set to Brooks Corey.\n";
                         
            relpermModel_ = new BrooksCorey<2U>(reservoirModel_->Database(), "brooks corey parameter", "entry pressure", false);
            break;
    }


    pressureSolver_ = new PressureSolver<2U,Element> (*reservoirModel_, setup_.WithGravitationalForces(), setup_.WithCapillaryForces());

    // creating the advector
    switch (setup_.Solution())
    {
    case IMPES:
        advector_ = new TwoPhaseElementBasedTransport<2U>(*reservoirModel_, *transportModel_, "Model", "porosity", "permeability",
                                                          "saturation oil", "saturation water", "divergence", "velocity", 
                                                          "oil volume source", "water volume source", "fluid volume source", setup_);
        break;

    case QEIMPES:
        advector_ = new ActiveElementTwoPhaseTransport<2U>(*reservoirModel_, *transportModel_, "Model", "porosity", "permeability",
                                                           "saturation oil", "saturation water", "divergence", "velocity", 
                                                           "oil volume source", "water volume source", "fluid volume source", setup_);
        break;

    case QFIMPES:
        advector_ = new ActiveFaceTwoPhaseTransport<2U>(*reservoirModel_, *transportModel_, "Model", "porosity", "permeability",
                                                        "saturation oil", "saturation water", "divergence", "velocity", 
                                                        "oil volume source", "water volume source", "fluid volume source", setup_);

        break;

    default:
        std::cout << "Warning! Generic2P2D_IMPES_Simulator::ConfigureModel()" <<
            " Solution method was not recognized. It is set to IMPES.\n";

        advector_ = new TwoPhaseElementBasedTransport<2U>(*reservoirModel_, *transportModel_, "Model", "porosity", "permeability",
                                                          "saturation oil", "saturation water", "divergence", "velocity", 
                                                          "oil volume source", "water volume source", "fluid volume source",setup_);
        break;
    }

} // end ConfigureModel





void Generic2P2D_IMPES_Simulator::InitializeProperties()
{
    cout << "\nGeneric2P2D_IMPES_Simulator::InitializeProperties()\n";
    
    double vmin, vmax;
    ScalarVariable zero_scalar(PLAIN, 0.);
    ScalarVariable temp;
    VectorVariable<2U> zero_vector(PLAIN, PLAIN, 0., 0.);
    Region<2U>&  gref(reservoirModel_->Region("Model"));
    
    reservoirModel_->MinMaxOf( "fluid volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'fluid volume source' has not been set for some or all elements. Zero will be assigned.\n";
    }

    reservoirModel_->MinMaxOf( "water volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'water volume source' has not been set for some or all elements. Zero will be assigned.\n";
    }

    reservoirModel_->MinMaxOf( "oil volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'oil volume source' has not been set for some or all elements. Zero will be assigned.\n";
    }

    reservoirModel_->MinMaxOf( "volume modifier", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'volume modifier' has not been set for some or all elements. One will be assigned.\n";
    }

    reservoirModel_->MinMaxOf( "nodal fluid volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'nodal fluid volume source' has not been set for some or all nodes. Zero will be assigned.\n";
    }

    // initializing uninitialized values
    for ( auto it = gref.CellsBegin(); it != gref.CellsEnd(); ++it)
    {
        // setting the capillary pressure term to zero
        (*it)->Store(propKeys_->ct_key, zero_vector);

        // setting fluid volume source
        temp = (*it)->Read(propKeys_->qf_key);
        if (::isnan(temp()))
        {
            (*it)->Store(propKeys_->qf_key, zero_scalar);
        }
        else if (temp() > 0.)
        {
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'fluid volume source' variable must be less than or equal to zero. It is used as a sink." );
        }
 
        // setting water volume source
        temp = (*it)->Read(propKeys_->qw_key);
        if (::isnan(temp()))
        {
            (*it)->Store(propKeys_->qw_key, zero_scalar);
        }
        else if (temp() < 0.)
        {
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'water volume source' variable must be greater than or equal to zero. It is used as a source." );
        }

        // setting oil volume source
        temp = (*it)->Read(propKeys_->qo_key);
        if (::isnan(temp()))
        {
            (*it)->Store(propKeys_->qo_key, zero_scalar);
        }
        else if (temp() < 0.)
        {
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'oil volume source' variable must be greater than or equal to zero. It is used as a source." );
        }
                   
    }
    
    // setting nodal fluid volume source
    for ( auto it = gref.NodesBegin(); it != gref.NodesEnd(); ++it)
    {
        temp = (*it)->Read(propKeys_->ns_key);
        if (::isnan(temp()))
        {
            (*it)->Store(propKeys_->ns_key, zero_scalar);
        }
    
    }

    // looping over interior point CVEs and initializing oil saturation using the average of their neighbors
    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin(); 
         cveit != transportModel_->PerimeterCVEsBegin(); ++cveit)
    {
        double so_k_phi_v(0.);
        double k_phi_v(0.);
        
        for (size_t n = 0U; n < cveit->Neighbors(); ++n)
        {
            k_phi_v += cveit->Neighbor(n)->E()->Read(propKeys_->k_key) * cveit->Neighbor(n)->E()->Read(propKeys_->phi_key) *
                       cveit->Neighbor(n)->Volume();
            so_k_phi_v += cveit->Neighbor(n)->E()->Read(propKeys_->so_key) * cveit->Neighbor(n)->E()->Read(propKeys_->k_key) *
                          cveit->Neighbor(n)->E()->Read(propKeys_->phi_key) * cveit->Neighbor(n)->Volume();
        }
        
        cveit->PropertyValue(1U, so_k_phi_v / k_phi_v);

    }


} // end InitializeProperties



/*
void Generic2P2D_IMPES_Simulator::InitializeProperties()
{
    cout << "\nGeneric2P2D_IMPES_Simulator::InitializeProperties()\n";
    
    double vmin, vmax;
    ScalarVariable  zero(PLAIN, 0.);
    VectorVariable<2U> zero_vec(PLAIN, PLAIN, 0., 0.);
    Region<2U,Element>&  gref(reservoirModel_->Region("Model"));
    
    // setting the capillary pressure term to zero
    for (vector<Element<2U>*>::iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); it++)
         (*it)->Store(propKeys_->ct_key, zero_vec);
    
    // looping over interior point CVEs and initializing oil saturation using the average of their neighbors
    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin(); 
         cveit != transportModel_->PerimeterCVEsBegin(); cveit++)
    {
        double so_k_phi_v(0.);
        double k_phi_v(0.);
        
        for (size_t n = 0U; n < cveit->Neighbors(); n++)
        {
            k_phi_v += cveit->Neighbor(n)->E()->Read(propKeys_->k_key) * cveit->Neighbor(n)->E()->Read(propKeys_->phi_key) *
                       cveit->Neighbor(n)->Volume();
            so_k_phi_v += cveit->Neighbor(n)->E()->Read(propKeys_->so_key) * cveit->Neighbor(n)->E()->Read(propKeys_->k_key) *
                          cveit->Neighbor(n)->E()->Read(propKeys_->phi_key) * cveit->Neighbor(n)->Volume();
        }
        
        cveit->PropertyValue(1U, so_k_phi_v / k_phi_v);

    }


    reservoirModel_->MinMaxOf( "fluid volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'fluid volume source' has not been set. Zero will be assigned.\n";

        for (vector<Element<2U>*>::iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); it++)
             (*it)->Store(propKeys_->qf_key, zero);
    }
    else
    {
        if (vmax > 0.)
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'fluid volume source' variable must be less than or equal to zero. It is used as a sink." );

    }
    
    reservoirModel_->MinMaxOf( "water volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'water volume source' has not been set. Zero will be assigned.\n";

        for (vector<Element<2U>*>::iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); it++)
             (*it)->Store(propKeys_->qw_key, zero);
    }
    else
    {
        if (vmin < 0.)
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'water volume source' variable must be greater than or equal to zero. It is used as a source." );

    }

    reservoirModel_->MinMaxOf( "oil volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'oil volume source' has not been set. Zero will be assigned.\n";

        for (vector<Element<2U>*>::iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); it++)
             (*it)->Store(propKeys_->qo_key, zero);
    }
    else
    {
        if (vmin < 0.)
            throw csmp::Exception( FATAL_ERROR, "Generic2P2D_IMPES_Simulator::InitializeProperties()",
                  "The 'oil volume source' variable must be greater than or equal to zero. It is used as a source." );

    }

    reservoirModel_->MinMaxOf( "nodal fluid volume source", vmin, vmax);
    if (::isnan(vmin) || ::isnan(vmax))
    {
        cout << "\nProperty 'nodal fluid volume source' has not been set. Zero will be assigned.\n";

        for (vector<Node<2U>*>::iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++)
             (*it)->Store(propKeys_->ns_key, zero);
    }

}

*/

/// SKM_FIX correct rounding of model time
void Generic2P2D_IMPES_Simulator::OutputToVTU( double time_unit )
{
  const double  n(ModelTime::Instance().modelTime/time_unit);
  const size_t  tstep(static_cast<long>(n >= 0. ? n + 0.5 : n - 0.5));
  static bool   first_time(true);

  if ( first_time )
  {
    vtuOutput_->OmitZeroInFileName(false);
    list<string> onceProps( setup_.OutputOnceBegin(), setup_.OutputOnceEnd() );
    vtuOutput_->OutputDataToVTU( "Initial Output", onceProps, "Model", tstep );
  }

  first_time = false;
  vtuOutput_->OmitZeroInFileName(false);
  list<string> alwaysProps( setup_.OutputAlwaysBegin(), setup_.OutputAlwaysEnd() );
  vtuOutput_->OutputDataToVTU( "Transient Output", alwaysProps, "Model", tstep );
}


/// SKM_FIX correct rounding of model time
void Generic2P2D_IMPES_Simulator::OutputToVTK( double time_unit )
{
  const double  n(ModelTime::Instance().modelTime/time_unit);
  const size_t  tstep(static_cast<long>(n >= 0. ? n + 0.5 : n - 0.5));
  static bool   first_time(true);

    if ( first_time )
    {
        for (std::set<std::string>::const_iterator sit = setup_.OutputOnceBegin(); sit != setup_.OutputOnceEnd(); ++sit)
        {
            vtkOutput_.OutputDataToVTK(*reservoirModel_, (*sit).c_str(), (*sit).c_str(), tstep);
        }
        
        first_time = false;
        
    }

    for (std::set<std::string>::const_iterator sit = setup_.OutputAlwaysBegin(); sit != setup_.OutputAlwaysEnd(); ++sit)
    {
        vtkOutput_.OutputDataToVTK(*reservoirModel_, (*sit).c_str(), (*sit).c_str(), tstep);
    }


/*

PropertyHandle<2U> active_element(reservoirModel_, "active element", SCALAR, ELEMENT);
Index active_element_key(reservoirModel_->Database().StorageKey("active element"));


for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->CVEsBegin();
    cveit != transportModel_->PointCVEsBegin(); cveit++)
{
    if (cveit->Active())
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 1.));
    }
    else
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 0.));
    }
}

for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterCVEsBegin();
    cveit != transportModel_->PerimeterPointCVEsBegin(); cveit++)
{
    if (cveit->Active())
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 1.));
    }
    else
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 0.));
    }
}

vtkOutput_.OutputDataToVTK(reservoirModel_, "active element", "active element", ModelTime::Instance().modelTime, false);

*/
/*
PropertyHandle<2U> active_node(reservoirModel_, "active node", SCALAR, NODE);
Index active_node_key(reservoirModel_->Database().StorageKey("active node"));
active_node = 0.;

for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin();
    cveit != transportModel_->PerimeterCVEsBegin(); cveit++)
{
    if (cveit->Active())
    {
        cveit->N(0U)->Store(active_node_key, ScalarVariable(PLAIN, 1.));
    }

}

vtkOutput_.OutputDataToVTK(reservoirModel_, "active node", "active node", ModelTime::Instance().modelTime, false);


std::set<ElementFace<2U>*> active_l_faces;
size_t counter_l(0U);
for (std::vector<ElementFace<2U>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->ActiveInteriorLineFacesBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->ActiveInteriorLineFacesEnd(); ++fit, ++counter_l)
{
    active_l_faces.insert(*fit);
}
for (std::set<ElementFace<2>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->MonitoringLineFacesBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->MonitoringLineFacesEnd(); ++fit, ++counter_l)
{
    active_l_faces.insert(*fit);
}

if (active_l_faces.size() != counter_l)
{
    std::cout << "line faces do not match " << active_l_faces.size() << " vs " << counter_l << "\n";
    system("pause");
}

std::set<ElementFace<2U>*> active_p_faces;
size_t counter_p(0U);
for (std::vector<ElementFace<2U>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->ActiveInteriorPointFacesBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->ActiveInteriorPointFacesEnd(); ++fit, ++counter_p, ++counter_l)
{
    active_p_faces.insert(*fit);
    active_l_faces.insert(*fit);
}
for (std::set<ElementFace<2>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->MonitoringPointFacesBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->MonitoringPointFacesEnd(); ++fit, ++counter_p, ++counter_l)
{
    active_p_faces.insert(*fit);
    active_l_faces.insert(*fit);
}

if (active_p_faces.size() != counter_p)
{
    std::cout << "point faces do not match " << active_p_faces.size() << " vs " << counter_p << "\n";
    system("pause");
}

if (active_l_faces.size() != counter_l)
{
    std::cout << "total faces do not match " << active_l_faces.size() << " vs " << counter_l << "\n";
    system("pause");
}


std::set<ControlVolumeElement<2U>*> active_2d_cve;
size_t counter_2d_cve(0U);
for (std::vector<ControlVolumeElement<2U>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->Active1D2D_CVEsBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->Active1D2D_CVEsEnd(); ++fit, ++counter_2d_cve)
{
    active_2d_cve.insert(*fit);
}
for (std::map<ControlVolumeElement<2>*, double>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->Monitoring1D2D_CVEsBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->Monitoring1D2D_CVEsEnd(); ++fit, ++counter_2d_cve)
{
    active_2d_cve.insert(fit->first);
}

if (active_2d_cve.size() != counter_2d_cve)
{
    std::cout << "2d cves do not match " << active_2d_cve.size() << " vs " << counter_2d_cve << "\n";
    system("pause");
}


std::set<ControlVolumeElement<2U>*> active_0d_cve;
size_t counter_0d_cve(0U);
for (std::vector<ControlVolumeElement<2U>*>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->Active0D_CVEsBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->Active0D_CVEsEnd(); ++fit, ++counter_0d_cve, ++counter_2d_cve)
{
    active_0d_cve.insert(*fit);
    active_2d_cve.insert(*fit);
}
for (std::map<ControlVolumeElement<2>*, double>::iterator fit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
     (advector_)->Monitoring0D_CVEsBegin(); fit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
     (advector_)->Monitoring0D_CVEsEnd(); ++fit, ++counter_0d_cve, ++counter_2d_cve)
{
    active_0d_cve.insert(fit->first);
    active_2d_cve.insert(fit->first);
}

if (active_0d_cve.size() != counter_0d_cve)
{
    std::cout << "2d cves do not match " << active_0d_cve.size() << " vs " << counter_0d_cve << "\n";
    system("pause");
}

if (active_2d_cve.size() != counter_2d_cve)
{
    std::cout << "total cves do not match " << active_2d_cve.size() << " vs " << counter_2d_cve << "\n";
    system("pause");
}
*/
 } // end OutputToVTK





 void Generic2P2D_IMPES_Simulator::ComputeFlowPropertiesAllRegions()
 {

    std::cout << "\nGeneric2P2D_IMPES_Simulator::ComputeFlowPropertiesAllRegions()\n";
    
    ScalarVariable  mob_t, sc;

    csmp::Region<2U>& sgref = reservoirModel_->Region("Model");

    // 0. Computing the saturation of water = 1 - So
    // ---------------------------------------------
    for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); ++it)
    {
        sc = 1. - (*it)->Read(propKeys_->so_key);
        (*it)->Store( propKeys_->sw_key, sc);
    }

    // 2. Computing the multiphase flow properties
    for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); ++it)
    {
        double vol_mod((*it)->Read(propKeys_->vm_key));
        
        // 0. setting up the relative permeability model
        // ---------------------------------------------
        relpermModel_->Initialize(*(*it));
        relpermModel_->EffectiveSaturation();

        // 1. k * total mobility
        // -------------------------------------------------------------------------------------------
        // in the NumIntegral_dNT_op_dN_dV operator there is a minus sign, so mob_t is multiplied by -1.
        mob_t = relpermModel_->TotalMobility() * vol_mod;
        (*it)->Store( propKeys_->tm_key, mob_t);

        // 2. gravity term: k g (lambda oil * rho oil + lambda water * rho water) for the PDE operator
        // -------------------------------------------------------------------------------------------
        if (setup_.WithGravitationalForces())
        {
            double temp(-9.8066 * relpermModel_->Permeability());
            temp *= (relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                     relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());
            temp *= vol_mod;
            VectorVariable<2U>  gt(PLAIN, PLAIN, 0., temp);

            // line elements
            if ( (*it)->FE()->IsLine() )
            {
                 // 1. find unit vector in direction of line element
                 Point<2U> evec((*it)->N(1U)->Coordinate() - (*it)->N(0U)->Coordinate()); 
                 // 2. check direction of evec and correct it to upward
                 if (evec[1U] < 0.) evec *= -1.;
                 evec /= evec.Length();
                 // 3. scaling element aligned unit vector by flow vector gravitational component
                 evec *= evec[1U] * temp;
                 // 4. add resulting velocity components to velocity vector
                 gt(0U) = evec[0U];
                 gt(1U) = evec[1U];
             }

            (*it)->Store(propKeys_->gt_key, gt);
        }

/*
if (setup_.WithCapillaryForces())
{
    VectorVariable<2U>  ct(PLAIN, PLAIN, 0., 0.);

    const double pc_e(relpermModel_->pc_Phase());
    const double perm_e(relpermModel_->Permeability());
    const double lambda_w_e(relpermModel_->MobilityPhase(1U));
    const Point<2U> bc_e((*it)->BaryCenter());
    

    for (size_t neighbor_no = 0U; neighbor_no < (*it)->Neighbors(); neighbor_no++)
    {
        
        if ((*it)->Neighbor(neighbor_no) != 0)
        {
            relpermModel_->Initialize(*((*it)->Neighbor(neighbor_no)));
            relpermModel_->EffectiveSaturation();

            const double pc_n(relpermModel_->pc_Phase());
            const double perm_n(relpermModel_->Permeability());
            const double lambda_w_n(relpermModel_->MobilityPhase(1U));
            const Point<2U> bc_n((*it)->Neighbor(neighbor_no)->BaryCenter());
            Point<2U> dpcdl(bc_e - bc_n);
            dpcdl.NormalizeLengthTo((pc_e - pc_n) / bc_e.DistanceTo(bc_n)); // capillary pressure derivative vector
            dpcdl *= (2. / (1. / perm_e + 1. / perm_n)) * 0.5 * (lambda_w_n + lambda_w_e);

            ct(0) += dpcdl[0U];
            ct(1) += dpcdl[1U];

        }
    
    }

    (*it)->Store( propKeys_->ct_key, ct );
}
*/        
    } // end looping over elements


        // for capillary term we loop over the faces
        
        // 3. capillary term: k * lambda water(upstream) * grad(pc) for the PDE operator
        // -------------------------------------------------------------------------------------------
        // for global formulation the capillary term doesn't appear in the pressure equation but for phase formulation
        // we should calculate capillary term using upstream water mobility and include in pressure equation
        if (setup_.WithCapillaryForces() && setup_.Formulation() == PHASE)
        {
          

for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); ++it)
    {
    VectorVariable<2U>  ct(PLAIN, PLAIN, 0., 0.);

    const double pc_e(relpermModel_->pc_Phase());
    const double perm_e(relpermModel_->Permeability());
    const double lambda_w_e(relpermModel_->MobilityPhase(1U));
    const Point<2U> bc_e((*it)->BaryCenter());

    for ( uint32_t neighbor_no = 0U; neighbor_no < (*it)->Neighbors(); ++neighbor_no)
    {
        
        if ((*it)->Neighbor(neighbor_no) != 0)
        {
            relpermModel_->Initialize(*((*it)->Neighbor(neighbor_no)));
            relpermModel_->EffectiveSaturation();

            const double pc_n(relpermModel_->pc_Phase());
            const double perm_n(relpermModel_->Permeability());
            const double lambda_w_n(relpermModel_->MobilityPhase(1U));
            const Point<2U> bc_n((*it)->Neighbor(neighbor_no)->BaryCenter());
            Point<2U> dpcdl(bc_e - bc_n);
            dpcdl.NormalizeLengthTo((pc_e - pc_n) / bc_e.DistanceTo(bc_n)); // capillary pressure derivative vector
            dpcdl *= (2. / (1. / perm_e + 1. / perm_n)) * 0.5 * (lambda_w_n + lambda_w_e);

            ct(0) += dpcdl[0U];
            ct(1) += dpcdl[1U];

        }

    }

    (*it)->Store( propKeys_->ct_key, ct );

}
                // solving the pressure equation without capillary term to find the upstream
                reservoirModel_->Apply(*pressureSolver_);
                // calculate total velocity
                ComputeTotalVelocity();
                // construct face fluxes to find upstream
                advector_->ConstructFluxes();
                
                
                VectorVariable<2U> zero_vec(PLAIN, PLAIN, 0., 0.);
                // setting the capillary pressure term to zero
                for ( auto it=sgref.CellsBegin(); it!=sgref.CellsEnd(); ++it)
                     (*it)->Store(propKeys_->ct_key, zero_vec);
                
                typedef std::vector<ElementFace<2U> >::iterator face_itr;
                
                for (face_itr fit = transportModel_->LineFacesBegin(); fit != transportModel_->PointFacesBegin(); ++fit)
                {                
                    relpermModel_->Initialize(*(fit->InsideCVE()->E()));
                    relpermModel_->EffectiveSaturation();
                    const double perm_in(relpermModel_->Permeability());
                    const double ro_w_in(relpermModel_->DensityWettingPhase());
                    const double ro_o_in(relpermModel_->DensityNonWettingPhase());
                    const double lambda_w_in(relpermModel_->MobilityPhase(1U));
                    const double lambda_o_in(relpermModel_->MobilityPhase(2U));
                    const double pc_in(relpermModel_->pc_Phase());
                    const double vol_mod_in(fit->InsideCVE()->E()->Read(propKeys_->vm_key));

                    relpermModel_->Initialize(*(fit->OutsideCVE()->E()));
                    relpermModel_->EffectiveSaturation();
                    const double perm_out(relpermModel_->Permeability());
                    const double ro_w_out(relpermModel_->DensityWettingPhase());
                    const double ro_o_out(relpermModel_->DensityNonWettingPhase());
                    const double lambda_w_out(relpermModel_->MobilityPhase(1U));
                    const double lambda_o_out(relpermModel_->MobilityPhase(2U));
                    const double pc_out(relpermModel_->pc_Phase());
                    const double vol_mod_out(fit->OutsideCVE()->E()->Read(propKeys_->vm_key));

                    const double capillary_grad(fit->PropertyValue(1U) * (pc_out - pc_in));
                    const double perm(2. / (1. / perm_out + 1. / perm_in));

                    double lambda_w, lambda_o;
                    double ro_w, ro_o;
                    advector_->FindUpstreamProperties(fit->PropertyValue(0U), fit->Area(), fit->UnitNormal(),
                                                      perm, lambda_w_out, lambda_o_out, lambda_w_in, lambda_o_in,
                                                      ro_w_out, ro_o_out, ro_w_in, ro_o_in, capillary_grad,
                                                      lambda_w, lambda_o, ro_w, ro_o);
                    
                   VectorVariable<2U> ct_in;
                   VectorVariable<2U> ct_out;
                   
                   fit->InsideCVE()->E()->Read(propKeys_->ct_key, ct_in);
                   fit->OutsideCVE()->E()->Read(propKeys_->ct_key, ct_out);

                   Point<2U> dpcdl(fit->OutsideCVE()->BaryCenter() - fit->InsideCVE()->BaryCenter());
                   dpcdl.NormalizeLengthTo((pc_out - pc_in) / fit->OutsideCVE()->BaryCenter().DistanceTo(fit->InsideCVE()->BaryCenter()));
                   dpcdl *= perm * lambda_w;
                   
                   ct_in(0) += dpcdl[0U] * vol_mod_in;
                   ct_in(1) += dpcdl[1U] * vol_mod_in;
                   ct_out(0) += dpcdl[0U] * vol_mod_out;
                   ct_out(1) += dpcdl[1U] * vol_mod_out;
                   
                   fit->InsideCVE()->E()->Store(propKeys_->ct_key, ct_in);
                   fit->OutsideCVE()->E()->Store(propKeys_->ct_key, ct_out);
                   
               }
            
//          }
            
//          first_time = false;
          
      }
      


} // end ComputeFlowPropertiesAllRegions





 void Generic2P2D_IMPES_Simulator::ComputeFlowPropertiesActiveRegions()
 {

    std::cout << "\nGeneric2P2D_IMPES_Simulator::ComputeFlowPropertiesActiveRegions()\n";
    
    ScalarVariable  mob_t, sc;

    if (setup_.Solution() == QEIMPES)
    {
        // Computing the multiphase flow properties active 1d2d cves
        for (std::vector<ControlVolumeElement<2U>*>::iterator cveit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
             (advector_)->Active1D2D_CVEsBegin(); cveit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
             (advector_)->Active1D2D_CVEsEnd(); ++cveit)
        {
            double vol_mod((*cveit)->E()->Read(propKeys_->vm_key));
            
            // 0. setting up the relative permeability model
            // ---------------------------------------------
            relpermModel_->Initialize(*((*cveit)->E()));
            relpermModel_->EffectiveSaturation();

            // 1. k * total mobility
            // -------------------------------------------------------------------------------------------
            // in the NumIntegral_dNT_op_dN_dV operator there is a minus sign, so mob_t is multiplied by -1.
            mob_t = relpermModel_->TotalMobility() * vol_mod;
            (*cveit)->E()->Store(propKeys_->tm_key, mob_t);

            // 2. gravity term: k g (lambda oil * rho oil + lambda water * rho water) for the PDE operator
            // -------------------------------------------------------------------------------------------
            if (setup_.WithGravitationalForces())
            {
                double temp(-9.8066 * relpermModel_->Permeability());
                temp *= (relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                         relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());
                temp *= vol_mod;
                VectorVariable<2U>  gt(PLAIN, PLAIN, 0., temp);

                // line elements
                if ((*cveit)->E()->FE()->IsLine()) 
                {
                     // 1. find unit vector in direction of line element
                     Point<2U> evec((*cveit)->E()->N(1U)->Coordinate() - (*cveit)->E()->N(0U)->Coordinate()); 
                     // 2. check direction of evec and correct it to upward
                     if (evec[1U] < 0.) evec *= -1.;
                     evec /= evec.Length();
                     // 3. scaling element aligned unit vector by flow vector gravitational component
                     evec *= evec[1U] * temp;
                     // 4. add resulting velocity components to velocity vector
                     gt(0U) = evec[0U];
                     gt(1U) = evec[1U];
                 }

                (*cveit)->E()->Store(propKeys_->gt_key, gt);
            }

        }

        // Computing the multiphase flow properties monitoring 1d2d cves
        for (std::map<ControlVolumeElement<2U>*, double>::iterator cveit = static_cast<ActiveElementTwoPhaseTransport<2U>*>
             (advector_)->Monitoring1D2D_CVEsBegin(); cveit != static_cast<ActiveElementTwoPhaseTransport<2U>*> 
             (advector_)->Monitoring1D2D_CVEsEnd(); ++cveit)
        {
            double vol_mod(cveit->first->E()->Read(propKeys_->vm_key));
            
            // 0. setting up the relative permeability model
            // ---------------------------------------------
            relpermModel_->Initialize(*(cveit->first->E()));
            relpermModel_->EffectiveSaturation();

            // 1. k * total mobility
            // -------------------------------------------------------------------------------------------
            // in the NumIntegral_dNT_op_dN_dV operator there is a minus sign, so mob_t is multiplied by -1.
            mob_t = relpermModel_->TotalMobility() * vol_mod;
            cveit->first->E()->Store(propKeys_->tm_key, mob_t);

            // 2. gravity term: k g (lambda oil * rho oil + lambda water * rho water) for the PDE operator
            // -------------------------------------------------------------------------------------------
            if (setup_.WithGravitationalForces())
            {
                double temp(-9.8066 * relpermModel_->Permeability());
                temp *= (relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                         relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());
                temp *= vol_mod;
                VectorVariable<2U>  gt(PLAIN, PLAIN, 0., temp);

                // line elements
                if (cveit->first->E()->FE()->IsLine())
                {
                     // 1. find unit vector in direction of line element
                     Point<2U> evec(cveit->first->E()->N(1U)->Coordinate() - cveit->first->E()->N(0U)->Coordinate()); 
                     // 2. check direction of evec and correct it to upward
                     if (evec[1U] < 0.) evec *= -1.;
                     evec /= evec.Length();
                     // 3. scaling element aligned unit vector by flow vector gravitational component
                     evec *= evec[1U] * temp;
                     // 4. add resulting velocity components to velocity vector
                     gt(0U) = evec[0U];
                     gt(1U) = evec[1U];
                 }

                cveit->first->E()->Store(propKeys_->gt_key, gt);
            }

        }
    
    } // if QEIMPES solution method
    

    if (setup_.Solution() == QFIMPES)
    {
        // Computing the multiphase flow properties active 1d2d cves
        for (std::vector<ControlVolumeElement<2U>*>::iterator cveit = static_cast<ActiveFaceTwoPhaseTransport<2U>*>
             (advector_)->Active1D2D_CVEsBegin(); cveit != static_cast<ActiveFaceTwoPhaseTransport<2U>*> 
             (advector_)->Active1D2D_CVEsEnd(); ++cveit)
        {
            double vol_mod((*cveit)->E()->Read(propKeys_->vm_key));
            
            // 0. setting up the relative permeability model
            // ---------------------------------------------
            relpermModel_->Initialize(*((*cveit)->E()));
            relpermModel_->EffectiveSaturation();

            // 1. k * total mobility
            // -------------------------------------------------------------------------------------------
            // in the NumIntegral_dNT_op_dN_dV operator there is a minus sign, so mob_t is multiplied by -1.
            mob_t = relpermModel_->TotalMobility() * vol_mod;
            (*cveit)->E()->Store(propKeys_->tm_key, mob_t);

            // 2. gravity term: k g (lambda oil * rho oil + lambda water * rho water) for the PDE operator
            // -------------------------------------------------------------------------------------------
            if (setup_.WithGravitationalForces())
            {
                double temp(-9.8066 * relpermModel_->Permeability());
                temp *= (relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                         relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());
                temp *= vol_mod;
                VectorVariable<2U>  gt(PLAIN, PLAIN, 0., temp);

                // line elements
                if ((*cveit)->E()->FE()->IsLine())
                {
                     // 1. find unit vector in direction of line element
                     Point<2U> evec((*cveit)->E()->N(1U)->Coordinate() - (*cveit)->E()->N(0U)->Coordinate()); 
                     // 2. check direction of evec and correct it to upward
                     if (evec[1U] < 0.) evec *= -1.;
                     evec /= evec.Length();
                     // 3. scaling element aligned unit vector by flow vector gravitational component
                     evec *= evec[1U] * temp;
                     // 4. add resulting velocity components to velocity vector
                     gt(0U) = evec[0U];
                     gt(1U) = evec[1U];
                 }

                (*cveit)->E()->Store(propKeys_->gt_key, gt);
            }

        }

    } // if QFIMPES solution method
    
} // end ComputeFlowPropertiesActiveRegions




void Generic2P2D_IMPES_Simulator::ComputeFlowProperties()
{
    static bool first_time(true);
     
    if (first_time)
    {
        ComputeFlowPropertiesAllRegions();
        first_time = false;
    }
    else
    {
        if (setup_.Solution() == IMPES)
        {
            ComputeFlowPropertiesAllRegions();
        }
        else
        {
            ComputeFlowPropertiesActiveRegions();
        }

    }
 
} // end ComputeFlowProperties





void Generic2P2D_IMPES_Simulator::ComputeTotalVelocity()
{
    VectorVariable<2U>  vt(PLAIN, PLAIN, 0., 0.);
    Region<2U>&  gref(reservoirModel_->Region("Model"));

    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it)
    {
         // total mobility at element barycenter
         double mob_t = (*it)->Read(propKeys_->tm_key);
         double vol_mod = (*it)->Read(propKeys_->vm_key);

         // interpolation function derivatives at element barycenter
         vt = 0.;
         //FiniteElementTraits<2U,Element> 
         (*(*it)).dN_AtBaryCenter(DERIV_, 1U);
         for (uint32_t i=0U; i<(*it)->Nodes(); ++i)
         {
             double pf = (*it)->N(i)->Read(propKeys_->pf_key);
	           vt(0) += -pf * DERIV_(0,i);
	           vt(1) += -pf * DERIV_(1,i);
         }
         
         vt *= mob_t / vol_mod;

         if (setup_.WithGravitationalForces())
         {
            VectorVariable<2U> gt;
            (*it)->Read(propKeys_->gt_key, gt);
            gt /= vol_mod;
            vt += gt;
 
         }

         if (setup_.WithCapillaryForces())
         {
            VectorVariable<2U> ct;
            (*it)->Read(propKeys_->ct_key, ct);
            ct /= vol_mod;
            vt += ct;
         }
 /*        
         VectorVariable<2U>  vtt(PLAIN, PLAIN, 0.01, 0.);
         vt = vtt;
*/
         // storing the computed velocity
        (*it)->Store(propKeys_->vt_key, vt);
     }


}



void Generic2P2D_IMPES_Simulator::ComputeOilVelocity()
{

    VectorVariable<2U>  vo(PLAIN, PLAIN, 0., 0.);
    Region<2U>&  gref(reservoirModel_->Region("Model"));


    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
    {
         relpermModel_->Initialize( *(*it) );
         relpermModel_->EffectiveSaturation();

         // total mobility at element barycenter
         double mob_o = -1. * relpermModel_->Permeability() * relpermModel_->MobilityPhase(2U);

         // interpolation function derivatives at element barycenter
         vo = 0.;
         //FiniteElementTraits<2U,Element> 
         (*(*it)).dN_AtBaryCenter(DERIV_, 1U);
         for ( auto i{0U}; i<(*it)->Nodes(); ++i)
         {
             double pf = (*it)->N(i)->Read( propKeys_->pf_key );
	           vo(0) += pf * DERIV_(0,i) * mob_o;
	           vo(1) += pf * DERIV_(1,i) * mob_o;
         }

         // for gravity
         if (setup_.WithGravitationalForces())
         {
            double temp(-9.8066 * relpermModel_->Permeability());
            temp *= (relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());

            // line elements
            if ( (*it)->FE()->IsLine() )
            {
                 // 1. find unit vector in direction of line element
                 Point<2U> evec((*it)->N(1U)->Coordinate() - (*it)->N(0U)->Coordinate()); 
                 // 2. check direction of evec and correct it to upward
                 if (evec[1U] < 0.) evec *= -1.;
                 evec /= evec.Length();
                 // 3. scaling element aligned unit vector by flow vector gravitational component
                 evec *= evec[1U] * temp;
                 // 4. add resulting velocity components to velocity vector
                 vo(0U) += evec[0U];
                 vo(1U) += evec[1U];
             }
             else
             {
                vo(1U) += temp;
             }
             

         }

         // no capillary term is needed for calculating oil velocity

         // storing the computed velocity
        (*it)->Store(propKeys_->vo_key, vo);
     }


}


void Generic2P2D_IMPES_Simulator::ComputeCapillaryPressure()
{

	std::cout << "\nGeneric2P2D_IMPES_Simulator::ComputeCapillaryPressure()\n";

	ScalarVariable  pc_e;

    csmp::Region<2U>& sgref = reservoirModel_->Region("Model");

	// Computing capillary pressure at each element
	// ---------------------------------------------
	for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); ++it)
	{
		relpermModel_->Initialize(*(*it));
		relpermModel_->EffectiveSaturation();
		pc_e = relpermModel_->pc_Phase();

        (*it)->Store(propKeys_->pc_key, pc_e );
	}
}




void Generic2P2D_IMPES_Simulator::CorrectPressureBoundaryConditionBasedonHydrostaticHead(bool top_as_reference)
{
    cout << "\nGeneric2P2D_IMPES_Simulator::CorrectPressureBoundaryConditionBasedonHydrostaticHead\n";
    
    static bool first_time(true);
    
    if (setup_.WithGravitationalForces())
    {

        static vector<size_t> left_boundary_nodes;
        static vector<size_t> right_boundary_nodes;
        static vector<pair<size_t, double> > left_boundary_elements;
        static vector<pair<size_t, double> > right_boundary_elements;
        
        csmp::Region<2U>& sgref = reservoirModel_->Region("Model");

        if (first_time)
        {
            size_t boundary_node(sgref.InteriorNodes());
            
            for ( auto nit=sgref.PerimeterNodesBegin(); nit!=sgref.NodesEnd(); ++nit, ++boundary_node)
            {
                if ((*nit)->AtBoundary() == LEFT || (*nit)->AtBoundary() == CNR1 || (*nit)->AtBoundary() == CNR4)
                {
                    if (left_boundary_nodes.empty())
                    {
                        left_boundary_nodes.push_back(boundary_node);
                        continue;
                    }
                    
                    if ((*nit)->Coordinate()[1U] < sgref.N(*left_boundary_nodes.begin())->Coordinate()[1U])
                    {
                        left_boundary_nodes.insert(left_boundary_nodes.begin(), boundary_node);
                        continue;
                    }
                    
                    for (vector<size_t>::iterator vit = left_boundary_nodes.begin(); vit < left_boundary_nodes.end(); ++vit)
                    {
                        if (vit != left_boundary_nodes.end() - 1)
                        {
                            if (((*nit)->Coordinate()[1U] > sgref.N((*vit))->Coordinate()[1U]) &&
                                ((*nit)->Coordinate()[1U] < sgref.N((*(vit + 1)))->Coordinate()[1U]))
                            {
                                left_boundary_nodes.insert(vit + 1, boundary_node);
                                break;
                            }
                            
                        }
                        else
                        {
                            left_boundary_nodes.insert(left_boundary_nodes.end(), boundary_node);
                            break;
                        }
                        
                    }
                    
                } // end LEFT boundary

                if ((*nit)->AtBoundary() == RIGHT || (*nit)->AtBoundary() == CNR2 || (*nit)->AtBoundary() == CNR3)
                {
                    if (right_boundary_nodes.empty())
                    {
                        right_boundary_nodes.push_back(boundary_node);
                        continue;
                    }
                    
                    if ((*nit)->Coordinate()[1U] < sgref.N(*right_boundary_nodes.begin())->Coordinate()[1U])
                    {
                        right_boundary_nodes.insert(right_boundary_nodes.begin(), boundary_node);
                        continue;
                    }
                    
                    for (vector<size_t>::iterator vit = right_boundary_nodes.begin(); vit < right_boundary_nodes.end(); ++vit)
                    {
                        if (vit != right_boundary_nodes.end() - 1)
                        {
                            if (((*nit)->Coordinate()[1U] > sgref.N((*vit))->Coordinate()[1U]) &&
                                ((*nit)->Coordinate()[1U] < sgref.N((*(vit + 1)))->Coordinate()[1U]))
                            {
                                right_boundary_nodes.insert(vit + 1, boundary_node);
                                break;
                            }
                            
                        }
                        else
                        {
                            right_boundary_nodes.insert(right_boundary_nodes.end(), boundary_node);
                            break;
                        }
                        
                    }
                    
                } // end RIGHT boundary
                
            } // end boundary nodes loop
            

            size_t boundary_element(sgref.InteriorCells());
            
            for ( auto eit=sgref.PerimeterCellsBegin(); eit!=sgref.CellsEnd(); ++eit, ++boundary_element)
            {
                Point<2U> N1, N2;
                BOX_BOUNDARY boundary(NOT);
                
                for ( uint32_t i = 0U; i < (*eit)->Nodes(); ++i)
                {                
                    if ((*eit)->N(i)->AtBoundary() == LEFT || (*eit)->N(i)->AtBoundary() == CNR1 || (*eit)->N(i)->AtBoundary() == CNR4)
                    {
                        N1 = (*eit)->N(i)->Coordinate();
                        boundary = LEFT;
                        break;
                    }

                    if ((*eit)->N(i)->AtBoundary() == RIGHT || (*eit)->N(i)->AtBoundary() == CNR2 || (*eit)->N(i)->AtBoundary() == CNR3)
                    {
                        N1 = (*eit)->N(i)->Coordinate();
                        boundary = RIGHT;
                        break;
                    }
                    
                }
                
                if (boundary == LEFT || boundary == RIGHT)
                {
                    for (uint32_t i = (*eit)->Nodes() - 1U; i >= 0U; i--)
                    {                
                        if ((*eit)->N(i)->AtBoundary() == LEFT || (*eit)->N(i)->AtBoundary() == CNR1 || (*eit)->N(i)->AtBoundary() == CNR4)
                        {
                            N2 = (*eit)->N(i)->Coordinate();
                            break;
                        }

                        if ((*eit)->N(i)->AtBoundary() == RIGHT || (*eit)->N(i)->AtBoundary() == CNR2 || (*eit)->N(i)->AtBoundary() == CNR3)
                        {
                            N2 = (*eit)->N(i)->Coordinate();
                            break;
                        }
                        
                    }
                
                }
                
                if (boundary == LEFT)
                {
                    if (left_boundary_elements.empty())
                    {
                        left_boundary_elements.push_back(make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                        continue;
                    }
                    
                    if (0.5*(N1[1U]+N2[1U]) < (*left_boundary_elements.begin()).second)
                    {
                        left_boundary_elements.insert(left_boundary_elements.begin(), make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                        continue;
                    }
                    
                    for (vector<pair<size_t, double> >::iterator vit = left_boundary_elements.begin(); vit < left_boundary_elements.end(); ++vit)
                    {
                        if (vit != left_boundary_elements.end() - 1)
                        {
                            if ((0.5*(N1[1U]+N2[1U]) > (*vit).second) &&
                                (0.5*(N1[1U]+N2[1U]) < (*(vit + 1)).second))
                            {
                                left_boundary_elements.insert(vit + 1, make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                                break;
                            }
                            
                        }
                        else
                        {
                            left_boundary_elements.insert(left_boundary_elements.end(), make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                            break;
                        }
                        
                    }
                
                }
                else if (boundary == RIGHT)
                {
                    if (right_boundary_elements.empty())
                    {
                        right_boundary_elements.push_back(make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                        continue;
                    }
                    
                    if (0.5*(N1[1U]+N2[1U]) < (*right_boundary_elements.begin()).second)
                    {
                        right_boundary_elements.insert(right_boundary_elements.begin(), make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                        continue;
                    }
                    
                    for (vector<pair<size_t, double> >::iterator vit = right_boundary_elements.begin(); vit < right_boundary_elements.end(); ++vit)
                    {
                        if (vit != right_boundary_elements.end() - 1)
                        {
                            if ((0.5*(N1[1U]+N2[1U]) > (*vit).second) &&
                                (0.5*(N1[1U]+N2[1U]) < (*(vit + 1)).second))
                            {
                                right_boundary_elements.insert(vit + 1, make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                                break;
                            }
                            
                        }
                        else
                        {
                            right_boundary_elements.insert(right_boundary_elements.end(), make_pair(boundary_element, 0.5*(N1[1U]+N2[1U])));
                            break;
                        }
                        
                    }
                
                }
                

            }
            
            first_time = false;
            
        } // end if first_time
        
        static const VARIABLE_FLAG left_flag(sgref.N(left_boundary_nodes[1U])->Status(propKeys_->pf_key));
        static const VARIABLE_FLAG right_flag(sgref.N(right_boundary_nodes[1U])->Status(propKeys_->pf_key));
        
        if (left_flag == DIRICH)
        {
            static const double pl(sgref.N(left_boundary_nodes[1U])->Read(propKeys_->pf_key));
                        
            if (top_as_reference)
            {
                sgref.N(left_boundary_nodes[left_boundary_nodes.size() - 1U])->Store(propKeys_->pf_key, ScalarVariable(left_flag, pl));

                for (size_t i = left_boundary_nodes.size() - 2U; i >= 0U; i--)
                {
                    if (i > left_boundary_nodes.size()) break;
                    
                    relpermModel_->Initialize(*sgref.E(left_boundary_elements[i].first));
                    relpermModel_->EffectiveSaturation();
                    
                    const double gradient((relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                                            relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase()) / 
                                            (relpermModel_->MobilityPhase(1U) + relpermModel_->MobilityPhase(2U)) * 9.8066);
                                            
                    const double p(sgref.N(left_boundary_nodes[i + 1U])->Read(propKeys_->pf_key) +
                                     (sgref.N(left_boundary_nodes[i + 1U])->Coordinate()[1U] - 
                                     sgref.N(left_boundary_nodes[i])->Coordinate()[1U]) * gradient);
                    
                    sgref.N(left_boundary_nodes[i])->Store(propKeys_->pf_key, ScalarVariable(left_flag, p));


                }
                
            }
            else
            {
                sgref.N(left_boundary_nodes[0U])->Store(propKeys_->pf_key, ScalarVariable(left_flag, pl));

                for (size_t i = 1U; i < left_boundary_nodes.size(); ++i)
                {
                    relpermModel_->Initialize(*sgref.E(left_boundary_elements[i - 1U].first));
                    relpermModel_->EffectiveSaturation();
                    
                    const double gradient((relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                                            relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase()) / 
                                            (relpermModel_->MobilityPhase(1U) + relpermModel_->MobilityPhase(2U)) * 9.8066);
                                            
                    const double p(sgref.N(left_boundary_nodes[i - 1U])->Read(propKeys_->pf_key) +
                                     (sgref.N(left_boundary_nodes[i - 1U])->Coordinate()[1U] - 
                                     sgref.N(left_boundary_nodes[i])->Coordinate()[1U]) * gradient);
                    
                    sgref.N(left_boundary_nodes[i])->Store(propKeys_->pf_key, ScalarVariable(left_flag, p));

                }

            }

        }

        if (right_flag == DIRICH)
        {
            static const double pr(sgref.N(right_boundary_nodes[1U])->Read(propKeys_->pf_key));
                        
            if (top_as_reference)
            {
                sgref.N(right_boundary_nodes[right_boundary_nodes.size() - 1U])->Store(propKeys_->pf_key, ScalarVariable(right_flag, pr));

                for (size_t i = right_boundary_nodes.size() - 2U; i >= 0U; i--)
                {
                    if (i > right_boundary_nodes.size()) break;
                    
                    relpermModel_->Initialize(*sgref.E(right_boundary_elements[i].first));
                    relpermModel_->EffectiveSaturation();
                    
                    const double gradient((relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                                            relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase()) / 
                                            (relpermModel_->MobilityPhase(1U) + relpermModel_->MobilityPhase(2U)) * 9.8066);
                                            
                    const double p(sgref.N(right_boundary_nodes[i + 1U])->Read(propKeys_->pf_key) +
                                     (sgref.N(right_boundary_nodes[i + 1U])->Coordinate()[1U] - 
                                     sgref.N(right_boundary_nodes[i])->Coordinate()[1U]) * gradient);
                    
                    sgref.N(right_boundary_nodes[i])->Store(propKeys_->pf_key, ScalarVariable(right_flag, p));


                }
                
            }
            else
            {
                sgref.N(right_boundary_nodes[0U])->Store(propKeys_->pf_key, ScalarVariable(right_flag, pr));

                for (size_t i = 1U; i < right_boundary_nodes.size(); ++i)
                {
                    relpermModel_->Initialize(*sgref.E(right_boundary_elements[i - 1U].first));
                    relpermModel_->EffectiveSaturation();
                    
                    const double gradient((relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase() +
                                            relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase()) / 
                                            (relpermModel_->MobilityPhase(1U) + relpermModel_->MobilityPhase(2U)) * 9.8066);
                                            
                    const double p(sgref.N(right_boundary_nodes[i - 1U])->Read(propKeys_->pf_key) +
                                     (sgref.N(right_boundary_nodes[i - 1U])->Coordinate()[1U] - 
                                     sgref.N(right_boundary_nodes[i])->Coordinate()[1U]) * gradient);
                    
                    sgref.N(right_boundary_nodes[i])->Store(propKeys_->pf_key, ScalarVariable(right_flag, p));

                }

            }
        
        }

    } // end if setup_.WithGravitationalForces()

}



void Generic2P2D_IMPES_Simulator::WriteRestartFile()
{
    
    std::string restart_file_name(modelName_);
    restart_file_name += ".rsd";

    std::string last_backup(modelName_);
    last_backup += ".bak";
    
    cout << "\nGeneric2P2D_IMPES_Simulator::WriteRestartFile() Writing restart file: " << restart_file_name << "\n";
    
    static bool first_time(true);
    
    if (first_time)
    {
        map<Point<2U>, size_t, pointcomp> barycenter_map;
        std::string restart_geo_file_name(modelName_);
        restart_geo_file_name += ".rsg";
        std::ofstream rsg_file(restart_geo_file_name.c_str());
        size_t counter(0U);
        
        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->CVEsBegin();
               cveit != transportModel_->PointCVEsBegin(); ++cveit, ++counter)
        {
            barycenter_map[cveit->E()->BaryCenter()] = counter;
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterCVEsBegin();
               cveit != transportModel_->PerimeterPointCVEsBegin(); ++cveit, ++counter)
        {
            barycenter_map[cveit->E()->BaryCenter()] = counter;
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin();
               cveit != transportModel_->PerimeterCVEsBegin(); ++cveit, ++counter)
        {
            barycenter_map[cveit->N(0U)->Coordinate()] = counter;
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterPointCVEsBegin();
               cveit != transportModel_->CVEsEnd(); ++cveit, ++counter)
        {
            barycenter_map[cveit->N(0U)->Coordinate()] = counter;
        }        
        
        // writing rsg file
        for (map<csmp::Point<2U>, size_t, pointcomp>::iterator it = barycenter_map.begin(); it != barycenter_map.end(); ++it)
        {
            rsg_file << it->second << "\n";
        }

        rsg_file.close();
        
        first_time = false;
        
    } // end if first_time

    std::ifstream restart_file_check;
    restart_file_check.open(restart_file_name.c_str());
    
    // deleting the backup if restart data exist
    if (restart_file_check)
    {
        remove(last_backup.c_str());
    }
    
    restart_file_check.close();    

    // renaming the last output to have a backup if something happened during writing new restart file
    rename(restart_file_name.c_str(), last_backup.c_str());
    
    std::ofstream restart_file(restart_file_name.c_str());
    
    restart_file << ModelTime::Instance().modelTime << "\n";
    
    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->CVEsBegin();
           cveit != transportModel_->PointCVEsBegin(); ++cveit)
    {
        double saturation(cveit->E()->Read(propKeys_->so_key));
        restart_file << saturation << "\n";
    }

    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterCVEsBegin();
           cveit != transportModel_->PerimeterPointCVEsBegin(); ++cveit)
    {
        double saturation(cveit->E()->Read(propKeys_->so_key));
        restart_file << saturation << "\n";
    }

    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin();
           cveit != transportModel_->PerimeterCVEsBegin(); ++cveit)
    {
        double saturation(cveit->PropertyValue(1U));
        restart_file << saturation << "\n";
    }

    for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterPointCVEsBegin();
           cveit != transportModel_->CVEsEnd(); ++cveit)
    {
        double saturation(cveit->PropertyValue(1U));
        restart_file << saturation << "\n";
    }
    
    restart_file.flush();
    
    restart_file.close();
    
    // deleting the backup
    //remove(last_backup.c_str());

} // end WriteRestartFile




void Generic2P2D_IMPES_Simulator::ReadRestartFile()
{
    
    std::string restart_geo_file_name(modelName_);
    restart_geo_file_name += ".rsg";

    std::string restart_dat_file_name(modelName_);
    restart_dat_file_name += ".rsd";
    
    cout << "\nGeneric2P2D_IMPES_Simulator::ReadRestartFile() Reading restart files: " << 
            restart_geo_file_name << " and " << restart_dat_file_name << "\n";

    std::ifstream restart_g_file;
    std::ifstream restart_d_file;
    
    restart_g_file.open(restart_geo_file_name.c_str());
    restart_d_file.open(restart_dat_file_name.c_str());
    
    if (!restart_g_file || !restart_d_file)
    {
        cout << "\nGeneric2P2D_IMPES_Simulator::ReadRestartFile() Can not read restart files!\n";
        ModelTime::Instance().modelTime = 0.;
    }
    else
    {
        restart_d_file >> ModelTime::Instance().modelTime;
        
        double data;
        map<size_t, double> data_map;
        
        for (size_t i = 0U; i < transportModel_->CVEs(); ++i)
        {
            restart_d_file >> data;
            data_map[i] = data;
        }

        restart_d_file.close();
                
        size_t point;
        map<size_t, double> converted_data_map;
        
        for (size_t i = 0U; i < transportModel_->CVEs(); ++i)
        {
            restart_g_file >> point;
            converted_data_map[i] = data_map[point];
        }
   
        restart_g_file.close();
        
        map<Point<2U>, ControlVolumeElement<2U>*, pointcomp> point_cve_map;
        
        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->CVEsBegin();
               cveit != transportModel_->PointCVEsBegin(); ++cveit)
        {
            point_cve_map[cveit->E()->BaryCenter()] = &(*cveit);
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterCVEsBegin();
               cveit != transportModel_->PerimeterPointCVEsBegin(); ++cveit)
        {
            point_cve_map[cveit->E()->BaryCenter()] = &(*cveit);
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin();
               cveit != transportModel_->PerimeterCVEsBegin(); ++cveit)
        {
            point_cve_map[cveit->N(0U)->Coordinate()] = &(*cveit);
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterPointCVEsBegin();
               cveit != transportModel_->CVEsEnd(); ++cveit)
        {
            point_cve_map[cveit->N(0U)->Coordinate()] = &(*cveit);
        }        
        
        // assigning the values to the corresponding cves
        map<size_t, double>::iterator data_itr = converted_data_map.begin();
        for (map<Point<2U>, ControlVolumeElement<2U>*, pointcomp>::iterator cveit = point_cve_map.begin(); 
             cveit != point_cve_map.end(); ++cveit, ++data_itr)
        {
            if (cveit->second->Dimension() != 0U)
            {
                cveit->second->E()->Store(propKeys_->so_key, ScalarVariable(PLAIN, data_itr->second));
                cveit->second->E()->Store(propKeys_->sw_key, ScalarVariable(PLAIN, 1. - data_itr->second));
            }
            else
            {
                cveit->second->PropertyValue(1U, data_itr->second);
            }
            
        }

        
        
/*        
        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->CVEsBegin();
               cveit != transportModel_->PointCVEsBegin(); ++cveit)
        {
            double saturation;
            restart_file >> saturation;
            cveit->E()->Store(propKeys_->so_key, ScalarVariable(PLAIN, saturation));
            cveit->E()->Store(propKeys_->sw_key, ScalarVariable(PLAIN, 1. - saturation));
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterCVEsBegin();
               cveit != transportModel_->PerimeterPointCVEsBegin(); ++cveit)
        {
            double saturation;
            restart_file >> saturation;
            cveit->E()->Store(propKeys_->so_key, ScalarVariable(PLAIN, saturation));
            cveit->E()->Store(propKeys_->sw_key, ScalarVariable(PLAIN, 1. - saturation));
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PointCVEsBegin();
               cveit != transportModel_->PerimeterCVEsBegin(); ++cveit)
        {
            double saturation;
            restart_file >> saturation;
            cveit->PropertyValue(1U, saturation);
        }

        for (std::vector<ControlVolumeElement<2U> >::iterator cveit = transportModel_->PerimeterPointCVEsBegin();
               cveit != transportModel_->CVEsEnd(); ++cveit)
        {
            double saturation;
            restart_file >> saturation;
            cveit->PropertyValue(1U, saturation);
        }
       
        
        restart_file.close();
 */        
    }
    
}



void Generic2P2D_IMPES_Simulator::WritePermAndVolmodValues()
{

    map<Point<2U>, size_t, pointcomp> barycenter_map;
    std::string geo_file_name("geo_file.txt");
    std::ofstream geo_file(geo_file_name.c_str());
    size_t counter(0U);

    Region<2U>&  gref(reservoirModel_->Region("Model")); //-------------------->

    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it, ++counter)
    {
        barycenter_map[(*it)->BaryCenter()] = counter;
    }
     

    // writing geo file
    for (map<csmp::Point<2U>, size_t, pointcomp>::iterator it = barycenter_map.begin(); it != barycenter_map.end(); ++it)
    {
        geo_file << it->second << "\n";
    }

    geo_file.close();

    std::ofstream permeability_file("permeabilities.txt");
    std::ofstream volume_modifier_file("volume modifiers.txt");


    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it)
    {
        double permeability((*it)->Read(propKeys_->k_key)); //------------------>
        double volume_modifier((*it)->Read(propKeys_->vm_key)); //------------------>
        permeability_file << permeability << "\n";
        volume_modifier_file << volume_modifier<< "\n";
    }


    permeability_file.flush();
    volume_modifier_file.flush();    
    permeability_file.close();
    volume_modifier_file.close();


} // end WritePermAndVolmodValues



void Generic2P2D_IMPES_Simulator::ReadPermVolmod()
{
    
    std::string geo_file_name("geo_file.txt");

    std::string perm_file_name("permeabilities.txt");
    std::string volmod_file_name("volume modifiers.txt");
    
    
    cout << "\nGeneric2P2D_IMPES_Simulator::ReadPermVolmod() Reading data files: " << 
            geo_file_name << ", " << perm_file_name << " and " << volmod_file_name << "\n";

    std::ifstream geo_file;
    std::ifstream perm_file;
    std::ifstream volmod_file;
    
    geo_file.open(geo_file_name.c_str());
    perm_file.open(perm_file_name.c_str());
    volmod_file.open(volmod_file_name.c_str());
    
    if (!geo_file || !perm_file || !volmod_file)
    {
        cout << "\nGeneric2P2D_IMPES_Simulator::ReadPermVolmod() Can not read files!\n";
        exit(1);

    }
    else
    {
       
        double data;
        map<size_t, double> perm_data_map, volmod_data_map;
        
        const size_t n_elmts = reservoirModel_->Region("Model").Cells();
        for (size_t i = 0U; i < n_elmts; ++i)
          {
              perm_file >> data;
              perm_data_map[i] = data;
              
              volmod_file >> data;
              volmod_data_map[i] = data;
          }

        perm_file.close();
        volmod_file.close();
                
        size_t point;
        map<size_t, double> converted_perm_data_map;
        map<size_t, double> converted_volmod_data_map;
        
        for (size_t i = 0U; i < n_elmts; ++i)
        {
            geo_file >> point;
            converted_perm_data_map[i] = perm_data_map[point];
            converted_volmod_data_map[i] = volmod_data_map[point];
        }
   
        geo_file.close();
        
        map<Point<2U>, Element<2U>*, pointcomp> point_ele_map;
        
        Region<2U>&  gref(reservoirModel_->Region("Model"));
        
        for ( auto it = gref.CellsBegin(); it != gref.CellsEnd(); ++it )
        {
            point_ele_map[(*it)->BaryCenter()] = *it;
        }
        
        // assigning the values to the corresponding cves
        map<size_t, double>::iterator perm_data_itr = converted_perm_data_map.begin();
        map<size_t, double>::iterator volmod_data_itr = converted_volmod_data_map.begin();
        
        for (map<Point<2U>, Element<2U>*, pointcomp>::iterator eit = point_ele_map.begin(); 
             eit != point_ele_map.end(); ++eit, ++perm_data_itr, ++volmod_data_itr)
        {
            eit->second->Store(propKeys_->k_key, ScalarVariable(PLAIN, perm_data_itr->second));
            eit->second->Store(propKeys_->vm_key, ScalarVariable(PLAIN, volmod_data_itr->second));
            
        }

       
    }
    
}



void Generic2P2D_IMPES_Simulator::ProcessPotentialSignals()
{
   SignalHandler sig;
      
   if (sig.OutputSignal())
   {
      OutputToVTU(1.);
      OutputFluidVolumes(ModelTime::Instance().modelTime, false);
      sig.OutputSignal(false);
   }
   
   if (sig.RestartSignal())
   {
      WriteRestartFile();
      sig.RestartSignal(false);
   }

   if (sig.QuitSignal())
   {
       // releasing the memory occupied by samg
#ifdef CSMP_WITH_SAMG_SOLVER
       int ierr;
       SAMG_REFRESH(&ierr);
#endif
       exit(1);
   }
   
}



void Generic2P2D_IMPES_Simulator::OutputVelocityValues(double time)
{
    
    std::cout << "\nGeneric2P2D_IMPES_Simulator::OutputVelocityValues()\n";
    
    map<std::string,Region<2U> >::const_iterator  grit=reservoirModel_->RegionsBegin();

    std::string file_name("velocities-");
    
    std::ostringstream strs;
    strs << time;
    file_name += strs.str();

    //file_name += time;
    file_name += ".txt";
    
    std::ofstream out_file(file_name.c_str(), ios::out | ios::app);
    
    VectorVariable<2U>  vt(PLAIN, PLAIN, 0., 0.);
    
	out_file << "vt\tvo\tvw\n";
    
    // go through all groups and output values
    // -------------------------------------------
    for ( grit=reservoirModel_->UniqueRegionsBegin(); grit!=reservoirModel_->UniqueRegionsEnd(); grit++ )
    {
    
        out_file << (*grit).first << "\n";
        
        for (vector<Element<2U>*>::const_iterator
             it=(*grit).second.CellsBegin(); it!=(*grit).second.CellsEnd(); it++)
        {

            // total mobility at element barycenter
            double mob_t = (*it)->Read(propKeys_->tm_key);
            double vol_mod = (*it)->Read(propKeys_->vm_key);

            // interpolation function derivatives at element barycenter
            vt = 0.;
            //FiniteElementTraits<2U,Element> 
            (*(*it)).dN_AtBaryCenter(DERIV_, 1U);
            for ( uint32_t i=0U; i<(*it)->Nodes(); ++i)
            {
                double pf = (*it)->N(i)->Read(propKeys_->pf_key);
                vt(0) += -pf * DERIV_(0,i);
                vt(1) += -pf * DERIV_(1,i);
            }

            vt *= mob_t / vol_mod;

            if (setup_.WithGravitationalForces())
            {
                VectorVariable<2U> gt;
                (*it)->Read(propKeys_->gt_key, gt);
                gt /= vol_mod;
                vt += gt;

            }

            if (setup_.WithCapillaryForces())
            {
                VectorVariable<2U> ct;
                (*it)->Read(propKeys_->ct_key, ct);
                ct /= vol_mod;
                vt += ct;
            }
            
            // output total velocity into the text file
            out_file << vt.Length() << "\t";


            relpermModel_->Initialize( *(*it) );
            relpermModel_->EffectiveSaturation();

            // oil mobility at element barycenter
            double mob_o = -1. * relpermModel_->Permeability() * relpermModel_->MobilityPhase(2U);

            // interpolation function derivatives at element barycenter
            vt = 0.;
            //FiniteElementTraits<2U,Element> 
            (*(*it)).dN_AtBaryCenter(DERIV_, 1U);
            for (uint32_t i=0U; i<(*it)->Nodes(); ++i)
            {
                double pf = (*it)->N(i)->Read( propKeys_->pf_key );
  	            vt(0) += pf * DERIV_(0,i) * mob_o;
  	            vt(1) += pf * DERIV_(1,i) * mob_o;
            }

            // for gravity
            if (setup_.WithGravitationalForces())
            {
               double temp(-9.8066 * relpermModel_->Permeability());
               temp *= (relpermModel_->MobilityPhase(2U) * relpermModel_->DensityNonWettingPhase());

               // line elements
               if ( (*it)->IsLine() )
               {
                    // 1. find unit vector in direction of line element
                    Point<2U> evec((*it)->N(1U)->Coordinate() - (*it)->N(0U)->Coordinate()); 
                    // 2. check direction of evec and correct it to upward
                    if (evec[1U] < 0.) evec *= -1.;
                    evec /= evec.Length();
                    // 3. scaling element aligned unit vector by flow vector gravitational component
                    evec *= evec[1U] * temp;
                    // 4. add resulting velocity components to velocity vector
                    vt(0U) += evec[0U];
                    vt(1U) += evec[1U];
                }
                else
                {
                   vt(1U) += temp;
                }
            
            }
            
            // output oil velocity into the text file
            out_file << vt.Length() << "\t";

            // water mobility at element barycenter
            double mob_w = -1. * relpermModel_->Permeability() * relpermModel_->MobilityPhase(1U);

            // interpolation function derivatives at element barycenter
            vt = 0.;
            //FiniteElementTraits<2U,Element> 
            (*(*it)).dN_AtBaryCenter(DERIV_, 1U);
            for ( uint32_t i=0U; i<(*it)->Nodes(); ++i)
            {
                double pf = (*it)->N(i)->Read( propKeys_->pf_key );
  	            vt(0) += pf * DERIV_(0,i) * mob_w;
  	            vt(1) += pf * DERIV_(1,i) * mob_w;
            }

            // for gravity
            if (setup_.WithGravitationalForces())
            {
               double temp(-9.8066 * relpermModel_->Permeability());
               temp *= (relpermModel_->MobilityPhase(1U) * relpermModel_->DensityWettingPhase());

               // line elements
               if ( (*it)->IsLine() )
               {
                    // 1. find unit vector in direction of line element
                    Point<2U> evec((*it)->N(1U)->Coordinate() - (*it)->N(0U)->Coordinate()); 
                    // 2. check direction of evec and correct it to upward
                    if (evec[1U] < 0.) evec *= -1.;
                    evec /= evec.Length();
                    // 3. scaling element aligned unit vector by flow vector gravitational component
                    evec *= evec[1U] * temp;
                    // 4. add resulting velocity components to velocity vector
                    vt(0U) += evec[0U];
                    vt(1U) += evec[1U];
                }
                else
                {
                   vt(1U) += temp;
                }
            
            }
            
            
            // output water velocity into the text file
            out_file << vt.Length() << "\n";
        
        }
         
    }
    
    out_file.close();

} // end OutputVelocityValues



void Generic2P2D_IMPES_Simulator::OutputInflowOutflow(double time, bool restart)
{
	static bool first_time(!restart);
	std::string file_name(modelName_);
	file_name += "-average-velocities.txt";

	std::ofstream out_file(file_name.c_str(), ios::out | ios::app);

	if (first_time)
	{
		out_file.close();
		remove(file_name.c_str());
		out_file.open(file_name.c_str(), ios::out | ios::app);

		out_file << "Time (s)";

		out_file << "\tOil Inflow Flux(m3/s)\tOil Inflow Area(m2)\tOil Inflow Velocity(m/s)"
				 << "\tOil Outflow Flux(m3/s)\tOil Outflow Area(m2)\tOil Outflow Velocity(m/s)"
				 << "\tWater Inflow Flux(m3/s)\tWater Inflow Area(m2)\tWater Inflow Velocity(m/s)"
				 << "\tWater Outflow Flux(m3/s)\tWater Outflow Area(m2)\tWater Outflow Velocity(m/s)"
				 << "\tTotal Inflow Flux(m3/s)\tTotal Inflow Area(m2)\tTotal Inflow Velocity(m/s)"
				 << "\tTotal Outflow Flux(m3/s)\tTotal Outflow Area(m2)\tTotal Outflow Velocity(m/s)";

		out_file << "\n";

		first_time = false;
	}

	out_file << time;

  Region<2U>& gref(reservoirModel_->Region("Model"));
	
	uint32_t inside_node;
	uint32_t outside_node;
	double   facet_area;

	enum{ OIL = 0, WATER = 1, TOTAL = 2 };

	std::vector<double> in_flow( 3, 0.0 );
	std::vector<double> out_flow(3, 0.0);
	std::vector<double> in_flow_region_area(3, 0.0);
	std::vector<double> out_flow_region_area(3, 0.0);
	std::vector<double> in_flow_velocity(3, 0.0);
	std::vector<double> out_flow_velocity(3, 0.0);
	std::vector<double> facet_flux(3, 0.0);

	VectorVariable<2U>   oil_velo;
	VectorVariable<2U>   water_velo;
	VectorVariable<2U>   total_velo;
  ScalarVariable       volume_modifier;

	vector< std::pair<double, double> > OIL_BFLUXES  (gref.Nodes(), std::pair<double, double>(0.0, 0.0));
	vector< std::pair<double, double> > WATER_BFLUXES(gref.Nodes(), std::pair<double, double>(0.0, 0.0));
	vector< std::pair<double, double> > TOTAL_BFLUXES(gref.Nodes(), std::pair<double, double>(0.0, 0.0));

	gref.UpdateMemberIndexes();

	bool perimeter_element(false);
	bool perimeter_node(false);

	// computing total flux balance by looping over all the elements in the Region
	for ( size_t n = 0; n<gref.Cells(); n++)
	{
		perimeter_element = gref.IsPerimeterCell(gref.E(n));
		perimeter_node = false;
		if (!perimeter_element)
		{
			for (uint32_t ni = 0; ni< gref.E(n)->Nodes(); ni++)
			{
				if (gref.IsPerimeterNode(gref.E(n)->N(ni)))
				{
					perimeter_node = true;
					break;
				}
			}
		}
		if (perimeter_element || perimeter_node)
		{
			// velocities
            gref.E(n)->Read(propKeys_->vo_key, oil_velo);
            gref.E(n)->Read(propKeys_->vt_key, total_velo);
            gref.E(n)->Read(propKeys_->vm_key, volume_modifier );
			water_velo = total_velo - oil_velo;

			// finite-volume facet related information
			for (uint32_t i = 0U; i<gref.E(n)->FV()->Facets(); i++)
			{
				// facet data
				gref.E(n)->FV()->FacetEdgeNodes(i, inside_node, outside_node);

          //if (!gref.IsPerimeterNode(gref.E(n)->N(inside_node)) || !gref.IsPerimeterNode(gref.E(n)->N(outside_node)))
          //{
          const auto idx_inside_node = gref.E(n)->N(inside_node)->Idx();
          const auto idx_outside_node = gref.E(n)->N(outside_node)->Idx();

          facet_area = gref.E(n)->FacetArea(i)*volume_modifier();

					// oil flux
					facet_flux[OIL] = gref.E(n)->ProjectionOnFacetNormal(i,oil_velo);
					facet_flux[OIL] *= facet_area;
					OIL_BFLUXES[idx_inside_node].first   += facet_flux[OIL];
					OIL_BFLUXES[outside_node].first  -= facet_flux[OIL];
				
					// water flux
					facet_flux[WATER] = gref.E(n)->ProjectionOnFacetNormal(i, water_velo);
					facet_flux[WATER] *= facet_area;
					WATER_BFLUXES[idx_inside_node].first   += facet_flux[WATER];
					WATER_BFLUXES[idx_outside_node].first  -= facet_flux[WATER];

					// total flux
					facet_flux[TOTAL] = gref.E(n)->ProjectionOnFacetNormal(i, total_velo);
					facet_flux[TOTAL] *= facet_area;
					TOTAL_BFLUXES[idx_inside_node].first   += facet_flux[TOTAL];
					TOTAL_BFLUXES[idx_outside_node].first  -= facet_flux[TOTAL];

					OIL_BFLUXES[idx_inside_node].second    += facet_area;
					OIL_BFLUXES[idx_outside_node].second   += facet_area;
					WATER_BFLUXES[idx_inside_node].second  += facet_area;
					WATER_BFLUXES[idx_outside_node].second += facet_area;
					TOTAL_BFLUXES[idx_inside_node].second  += facet_area;
					TOTAL_BFLUXES[idx_outside_node].second += facet_area;
                //}
			}
		}
	}

	//  calculating the in and outfluxes through model boundaries
	fill( in_flow.begin(), in_flow.end(), 0.0 );
	fill( in_flow_region_area.begin(), in_flow_region_area.end(), 0.0);
	fill( in_flow_velocity.begin(), in_flow_velocity.end(), 0.0);

    const double zero( numeric_limits<double>::epsilon() );

	for (size_t i = gref.InteriorNodes(); i<gref.Nodes(); i++)
	{
		if ( OIL_BFLUXES[i].first > zero )
		{
			in_flow[OIL]			  += OIL_BFLUXES[i].first;
			in_flow_region_area[OIL]  += OIL_BFLUXES[i].second;
		}
		else if ( OIL_BFLUXES[i].first < zero )
		{
			out_flow[OIL]			  -= OIL_BFLUXES[i].first;
			out_flow_region_area[OIL] += OIL_BFLUXES[i].second;
		}

		if ( WATER_BFLUXES[i].first > zero)
		{
			in_flow[WATER]             += WATER_BFLUXES[i].first;
			in_flow_region_area[WATER] += WATER_BFLUXES[i].second;
		}
		else if ( WATER_BFLUXES[i].first < zero)
		{
			out_flow[WATER]             -= WATER_BFLUXES[i].first;
			out_flow_region_area[WATER] += WATER_BFLUXES[i].second;
		}

		if (TOTAL_BFLUXES[i].first > zero)
		{
			in_flow[TOTAL]             += TOTAL_BFLUXES[i].first;
			in_flow_region_area[TOTAL] += TOTAL_BFLUXES[i].second;
		}
		else if (TOTAL_BFLUXES[i].first < zero)
		{
			out_flow[TOTAL]             -= TOTAL_BFLUXES[i].first;
			out_flow_region_area[TOTAL] += TOTAL_BFLUXES[i].second;
		}
	}

	for ( size_t i = 0; i < 3; i++ )
	{
		in_flow_velocity[i]  = in_flow[i] / in_flow_region_area[i];
		out_flow_velocity[i] = out_flow[i] / out_flow_region_area[i];

		out_file << "\t" << in_flow[i]  << "\t" << in_flow_region_area[i]  << "\t" << in_flow_velocity[i]
			     << "\t" << out_flow[i] << "\t" << out_flow_region_area[i] << "\t" << out_flow_velocity[i];
	}

	out_file << "\n";

	out_file.close();


} // end OutputInflowOutflow

void Generic2P2D_IMPES_Simulator::OutputFluidVolumes(double time, bool restart)
{
    
    std::cout << "\nGeneric2P2D_IMPES_Simulator::OutputFluidVolumes()\n";
    
    static bool first_time(!restart);
    std::string file_name(modelName_);
    file_name += "-FIP.txt";
    
    std::ofstream out_file(file_name.c_str(), ios::out | ios::app);
    
    if (first_time)
    {
        out_file.close();
        remove (file_name.c_str());
        out_file.open(file_name.c_str(), ios::out | ios::app);
        
        out_file << "Time (s)\t";
        
        for (std::map<std::string,csmp::Region<2U> >::const_iterator rit = reservoirModel_->UniqueRegionsBegin();
             rit != reservoirModel_->UniqueRegionsEnd(); ++rit)
        {
            out_file << rit->first << "'s OIP (m3)\t" << rit->first << "'s WIP (m3)\t";
        }

//out_file << "pressure min\tpressure max";
        
        out_file << "\n";
        
        first_time = false;
        
    }

    out_file << time << "\t";
    
    for ( auto rit = reservoirModel_->UniqueRegionsBegin();
         rit != reservoirModel_->UniqueRegionsEnd(); ++rit)
    {
        double oil_volume(0.);
        double water_volume(0.);
        
        for ( auto eit = rit->second.CellsBegin(); eit != rit->second.CellsEnd(); ++eit)
        {
            const double volume((*eit)->Volume());
            const double volume_modifier((*eit)->Read(propKeys_->vm_key));
            const double water_saturation((*eit)->Read(propKeys_->sw_key));
            const double porosity((*eit)->Read(propKeys_->phi_key));
            const double oil_saturation((*eit)->Read(propKeys_->so_key));
            
            oil_volume += volume * volume_modifier * porosity * oil_saturation;
            water_volume += volume * volume_modifier * porosity * water_saturation;
        }
        
        out_file << oil_volume << "\t" << water_volume << "\t";
    }
    
    out_file << "\n";
  
    out_file.close();

} // end OutputFluidVolumes


/**
Writes the fluid pressure from all regions into seperate texfile with extension -fluid-pressure.txt
*/
void Generic2P2D_IMPES_Simulator::OutputFluidPressure(double time, bool restart) const
{

	std::cout << "\nGeneric2P2D_IMPES_Simulator::OutputFluidPressure()\n";
	static bool first_time(!restart);
	std::string file_name(modelName_);
	file_name += "-fluid-pressure.txt";

	std::ofstream out_file(file_name.c_str(), ios::out | ios::app);

	if (first_time)
	{
		out_file.close();
		remove(file_name.c_str());
		out_file.open(file_name.c_str(), ios::out | ios::app);

		out_file << "Time (s)\t";

        for (std::map<std::string, csmp::Region<2U> >::const_iterator rit = reservoirModel_->UniqueRegionsBegin();
            rit != reservoirModel_->UniqueRegionsEnd(); ++rit)
		{
			out_file << rit->first << "'s min fluid pressure (Pa)\t" << rit->first << "'s max fluid pressure (Pa)\t";
		}

		out_file << "\n";

		first_time = false;

	}

	out_file << time << "\t";

	vector<Node<2U>*>::const_iterator it;

    for (std::map<std::string, csmp::Region<2U> >::const_iterator rit = reservoirModel_->UniqueRegionsBegin();
        rit != reservoirModel_->UniqueRegionsEnd(); ++rit)
	{
		it = rit->second.NodesBegin();
        double min_p((*it)->Read(propKeys_->pf_key));
        double max_p((*it)->Read(propKeys_->pf_key));

		while (it != rit->second.NodesEnd())
		{
            double pf((*it)->Read(propKeys_->pf_key));
			if (pf > max_p) max_p = pf;
			if (pf < min_p) min_p = pf;
			it++;
		}

		out_file << min_p << "\t" << max_p << "\t";
	}

	out_file << "\n";

	out_file.close();

} // end OutputFluidPressures


/** 
     Writes the fluid pressure in the WELLF region into seperate texfile with extension -WELLF-fluid-pressure.txt
*/
void Generic2P2D_IMPES_Simulator::OutputWellF_FluidPressure( double time, bool restart ) const
{
    
    std::cout << "\nGeneric2P2D_IMPES_Simulator::OutputWellF_FluidPressure\n";
    static bool first_time(!restart);
    std::string file_name(modelName_);
    file_name += "-WELLF-fluid-pressure.txt";
    
    std::ofstream out_file(file_name.c_str(), ios::out | ios::app);
    
    if (first_time)
    {
        out_file.close();
        remove (file_name.c_str());
        out_file.open(file_name.c_str(), ios::out | ios::app);
        
        out_file << "Time (s)\t" << "min fluid pressure (Pa)\t" << "'max fluid pressure (Pa)\t in region WELLF\n";
      
        first_time = false;
    }


    // write out pressure extrema in mandatory WELLF region
    const Region<2U>&  gref(reservoirModel_->Region("WELLF"));

     vector<Node<2U>*>::const_iterator it=gref.NodesBegin();
     double min_p((*it)->Read(propKeys_->pf_key)), max_p((*it)->Read(propKeys_->pf_key));

    while ( it != gref.NodesEnd() )
    {
        double pf((*it)->Read(propKeys_->pf_key));
        if (pf > max_p) max_p = pf;
        if (pf < min_p) min_p = pf;
        it++;
    }

    out_file << time <<"\t"<< min_p << "\t" << max_p << "\n";
  
    out_file.close();

} // end OutputWellF_FluidPressure








void Generic2P2D_IMPES_Simulator::OutputCapillaryPressureGradient(double time)
{
    
    std::cout << "\nGeneric2P2D_IMPES_Simulator::OutputCapillaryPressureGradient()\n";
    
    std::string file_name("CapillaryPressureGradient-");
    
    std::ostringstream strs;
    strs << time;
    file_name += strs.str();

    file_name += ".txt";
    
    std::ofstream out_file(file_name.c_str(), ios::out);
    
    out_file << "x1\ty1\tx2\ty2\tpc1\tpc2\tsw1\tsw2\n";
    
    
    for (std::vector<ElementFace<2U> >::iterator fit = transportModel_->LineFacesBegin();
         fit != transportModel_->PointFacesBegin(); ++fit)
    {

        relpermModel_->Initialize(*(fit->InsideCVE()->E()));
        relpermModel_->EffectiveSaturation();
        const double pc_i(relpermModel_->pc_Phase());
        const double sw_i(relpermModel_->Saturation(1U));
        const Point<2U> bc_i(fit->InsideCVE()->E()->BaryCenter());

        relpermModel_->Initialize(*(fit->OutsideCVE()->E()));
        relpermModel_->EffectiveSaturation();
        const double pc_o(relpermModel_->pc_Phase());
        const double sw_o(relpermModel_->Saturation(1U));
        const Point<2U> bc_o(fit->OutsideCVE()->E()->BaryCenter());

        out_file << bc_i[0] << "\t" << bc_i[1] << "\t" 
                 << bc_o[0] << "\t" << bc_o[1] << "\t"
                 << pc_i << "\t" << pc_o << "\t"
                 << sw_i << "\t" << sw_o << "\n";
    
    }
    
    out_file.close();

    
} // end OutputCapillaryPressureGradient
    
    
    
    


void Generic2P2D_IMPES_Simulator::SolvePressure()
{
    static bool   first_time{true};
    static size_t counter{0};

    if (first_time)
      {
#ifdef CSMP_WITH_SAMG_SOLVER
          if (setup_.SAMG_Reconstruction() != 0U)
            {
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_itypu(1);
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(4);
                pressureSolver_->IntegrateOver(reservoirModel_->Region("Model"));
                // using the previous solution as initial guess for succeeding samg calls
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_itypu(0);
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(3);
                ++counter;
                first_time = false;
            }
          else {
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_itypu(1);
                static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(5);
                pressureSolver_->IntegrateOver(reservoirModel_->Region("Model"));
                first_time = false;
            }
#else
        pressureSolver_->IntegrateOver(reservoirModel_->Region("Model") );
        first_time = false;
#endif
    } 
   else // once the first solution has been obtained
    {
#ifdef CSMP_WITH_SAMG_SOLVER
        if (setup_.SAMG_Reconstruction() != 0U)
          {
              if (counter % setup_.SAMG_Reconstruction() == 0U )
                {
                    static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(4);
                    pressureSolver_->IntegrateOver(reservoirModel_->Region("Model"));
                    static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(3);
                    counter = 1U;
                }
              else {
                    pressureSolver_->IntegrateOver(reservoirModel_->Region("Model"));
                    ++counter;
                }
          }
        else {
            static_cast<PressureSolver<2U,Element>*> (pressureSolver_)->PressureSolverSettings().Set_iswit(5);
            pressureSolver_->IntegrateOver(reservoirModel_->Region("Model"));
        }
#else
        pressureSolver_->IntegrateOver(reservoirModel_->Region("Model") );
#endif
    }

} // end SolvePressure




void Generic2P2D_IMPES_Simulator::Run(bool restart)
{
  ConfigureModel();
  InitializeProperties();
  
  if (restart)
    {
        ReadRestartFile();
        // GIVE USER THE OPTION TO READ IN A DIFFERENT CONFIGURATION FILE
        Standard_IO_Handler stdio;
        if ( stdio.YesNo("Do you want to read a '*-restart-configuration.txt' file") ) {
        const std::string restart_config_file( modelName_ + "-restart" );
        InputDataManager<2U>  model_configuration;
        // TODO: this config call is different to the one in the main file
        model_configuration.ConfigureFromFile(*reservoirModel_,
                                               restart_config_file.c_str(),
                                               false,
                                               true,   // 2) default prop.values
                                               true,   // 3) group prop.values
                                               true,   // 4) essential conditions
                                               true,   // 5) essential flags
                                               false,  // 6) arbitrary model shape boundary conditions
                                               runSettings_);
         }
    }
  else ModelTime::Instance().modelTime = 0.;
 
  reservoirModel_->InstantiateFiniteVolumes();

  setPropertyToZero( *reservoirModel_, "nodal fluid volume source" );
  
  volumeWeightedDistributionOfElementPropertyToNodeAndAdd(*reservoirModel_,
                                                          "fluid volume source",
                                                          "nodal fluid volume source",
                                                          "volume modifier",
                                                          "Model");
  volumeWeightedDistributionOfElementPropertyToNodeAndAdd(*reservoirModel_,
                                                          "water volume source",
                                                          "nodal fluid volume source",
                                                          "volume modifier",
                                                          "Model");
  volumeWeightedDistributionOfElementPropertyToNodeAndAdd(*reservoirModel_,
                                                          "oil volume source",
                                                          "nodal fluid volume source",
                                                          "volume modifier",
                                                          "Model");
  // adding boundary fluxes to nodes
  areaWeightedDistributionOfNodePropertyToNodeAndAdd(*reservoirModel_,
                                                     "boundary source",
                                                     "nodal fluid volume source");
  
  if (setup_.WithGravitationalForces()) CorrectPressureBoundaryConditionBasedonHydrostaticHead(true);

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING

  ComputeCapillaryPressure();

#endif

  ComputeFlowProperties();

  SolvePressure();
  
  ComputeTotalVelocity();
  ComputeOilVelocity();

  OutputToVTU(1.);
  OutputFluidVolumes(ModelTime::Instance().modelTime, restart);
  OutputFluidPressure(ModelTime::Instance().modelTime, restart);
  //OutputWellF_FluidPressure( ModelTime::Instance().modelTime, restart );

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING

  OutputVelocityValues(ModelTime::Instance().modelTime); // ------------------------------->
  OutputCapillaryPressureGradient(ModelTime::Instance().modelTime); // ------------------------------->

#endif

#ifdef FECFV_SIMULATOR_INFLOW_OUTFLOW

  OutputInflowOutflow(ModelTime::Instance().modelTime,restart);

#endif

  // TwoPhaseElementBasedTransport
  advector_->AssignTransportBoundaryCondition("saturation");

  double truncated_timestep;
  
  while (ModelTime::Instance().modelTime < runSettings_.Duration())
  {
      const double output_time(runSettings_.NearestOutputTime(ModelTime::Instance().modelTime) );
      
      while (ModelTime::Instance().modelTime < output_time)
      {
          truncated_timestep = runSettings_.TimeToNearestOutputTime(ModelTime::Instance().modelTime);

          std::cout << "\nModel Time: " << ModelTime::Instance().modelTime << "\n";

          ModelTime::Instance().modelTime += advector_->Transport(truncated_timestep, *relpermModel_);

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING

          ComputeCapillaryPressure();

#endif
          // if model has received a signal
          ProcessPotentialSignals();

          if (setup_.WithGravitationalForces()) CorrectPressureBoundaryConditionBasedonHydrostaticHead(true);
         
          ComputeFlowProperties();

          SolvePressure();

          ComputeTotalVelocity();
          ComputeOilVelocity();
        
      }

      OutputToVTU(1);
      OutputFluidVolumes(ModelTime::Instance().modelTime, restart);
      OutputFluidPressure(ModelTime::Instance().modelTime, restart);
      //OutputWellF_FluidPressure( ModelTime::Instance().modelTime, restart );

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING
      
      OutputVelocityValues(ModelTime::Instance().modelTime); // ------------------------------->
      OutputCapillaryPressureGradient(ModelTime::Instance().modelTime); // ------------------------------->

#endif

#ifdef FECFV_SIMULATOR_INFLOW_OUTFLOW

      OutputInflowOutflow(ModelTime::Instance().modelTime,restart);

#endif

      WriteRestartFile();

  }

  // releasing the memory occupied by samg
#ifdef CSMP_WITH_SAMG_SOLVER
  int ierr;
  SAMG_REFRESH(&ierr);
#endif

} // end Run








void Generic2P2D_IMPES_Simulator::CalculateEffectivePermeability()
{

    ConfigureModel();
    InitializeProperties();
    
    setup_.Solution(IMPES);
#ifdef CSMP_WITH_SAMG_SOLVER
    setup_.SAMG_Reconstruction(0U);
#endif
    setup_.WithCapillaryForces(false);
    setup_.WithGravitationalForces(false);

    Point<2U> xyz_min, xyz_max; 
    reservoirModel_->MinMaxCoordinates( xyz_min, xyz_max );
    double Lx(std::abs(xyz_max[0] - xyz_min[0]));
    double Ly(std::abs(xyz_max[1] - xyz_min[1]));
    const double dp(1.0e+5);

    reservoirModel_->Region("Model").InputPropertyValue("viscosity oil", ScalarVariable(PLAIN,1.0e-3));
    reservoirModel_->Region("Model").InputPropertyValue("viscosity water", ScalarVariable(PLAIN,1.0e-3));
    reservoirModel_->Region("Model").InputPropertyValue("saturation oil", ScalarVariable(PLAIN,1.));
    reservoirModel_->Region("Model").InputPropertyValue("saturation water", ScalarVariable(PLAIN,0.));

    reservoirModel_->InputBoundaryValue(RIGHT, "fluid pressure", ScalarVariable(DIRICH, 1.0e+6));
    reservoirModel_->InputBoundaryValue(LEFT, "fluid pressure", ScalarVariable(DIRICH, 1.0e+6 + Lx * dp));
    reservoirModel_->InputBoundaryValue(TOP, "fluid pressure", ScalarVariable(PLAIN, 0.));
    reservoirModel_->InputBoundaryValue(BOTTOM, "fluid pressure", ScalarVariable(PLAIN, 0.));

    ComputeFlowProperties();
    SolvePressure();
    ComputeTotalVelocity();

    advector_->ConstructFluxes();

    double left_flux(0.);
    double right_flux(0.);

    for (std::vector<ElementFace<2U> >::iterator fit = transportModel_->PerimeterLineFacesBegin(); 
         fit != transportModel_->FacesEnd(); ++fit)
    {
        if (fit->Placement() == LEFT)
        {
            left_flux += fit->PropertyValue(0U);
        }

        if (fit->Placement() == RIGHT)
        {
            right_flux += fit->PropertyValue(0U);
        }

    }

    reservoirModel_->InputBoundaryValue(TOP, "fluid pressure", ScalarVariable(DIRICH, 1.0e+6));
    reservoirModel_->InputBoundaryValue(BOTTOM, "fluid pressure", ScalarVariable(DIRICH, 1.0e+6 + Ly * dp));
    reservoirModel_->InputBoundaryValue(LEFT, "fluid pressure", ScalarVariable(PLAIN, 0.));
    reservoirModel_->InputBoundaryValue(RIGHT, "fluid pressure", ScalarVariable(PLAIN, 0.));

    SolvePressure();
    ComputeTotalVelocity();

    advector_->ConstructFluxes();

    double top_flux(0.);
    double bottom_flux(0.);

    for ( auto fit = transportModel_->PerimeterLineFacesBegin();
          fit != transportModel_->FacesEnd(); ++fit)
    {
        if (fit->Placement() == TOP)
        {
            top_flux += fit->PropertyValue(0U);
        }

        if (fit->Placement() == BOTTOM)
        {
            bottom_flux += fit->PropertyValue(0U);
        }

    }

    cout << "\n\nHorizontal Effective Permeability Calculation:\n";
    cout << "----------------------------------------------\n\n";
    cout << "Model Length:          " << Ly << "\n";
    cout << "Pressure Difference:   " << Lx * dp << "\n";
    cout << "Viscosity:             " << 1.e-3 << "\n";
    cout << "LEFT Flux:             " << std::abs(left_flux) << "\n";
    cout << "RIGHT Flux:            " << std::abs(right_flux) << "\n";
    cout << "Average Flux:          " << 0.5 * (std::abs(left_flux) + std::abs(right_flux)) << "\n";
    cout << "KH Effective:          " << 0.5 * (std::abs(left_flux) + std::abs(right_flux)) * 1.e-3 / (Ly * dp) << "\n";
    cout << "\n\n";


    cout << "Vertical Effective Permeability Calculation:\n";
    cout << "--------------------------------------------\n\n";
    cout << "Model Height:          " << Lx << "\n";
    cout << "Pressure Difference:   " << Ly * dp << "\n";
    cout << "Viscosity:             " << 1.e-3 << "\n";
    cout << "TOP Flux:              " << std::abs(top_flux) << "\n";
    cout << "BOTTOM Flux:           " << std::abs(bottom_flux) << "\n";
    cout << "Average Flux:          " << 0.5 * (std::abs(top_flux) + std::abs(bottom_flux)) << "\n";
    cout << "KV Effective:          " << 0.5 * (std::abs(top_flux) + std::abs(bottom_flux)) * 1.e-3 / (Lx * dp) << "\n";


}



void Generic2P2D_IMPES_Simulator::FindContinuousLineRegionContainsCVE(ControlVolumeElement<2U>* elm, 
                                                                      std::set<ControlVolumeElement<2U>* >& visited,
                                                                      std::vector<ControlVolumeElement<2U>* >& discovered)
{
    for (size_t n = 0; n < elm->Neighbors(); ++n)
    {
        if ((elm->Neighbor(n)->Dimension() == 1) && (visited.find(elm->Neighbor(n)) == visited.end()))
        {
            visited.insert(elm->Neighbor(n));
            discovered.push_back(elm->Neighbor(n));
            FindContinuousLineRegionContainsCVE(elm->Neighbor(n), visited, discovered);
        }
        
    }
    
}



static double normalDistributionGenerator( double mean, double sd, double minimum, double maximum )
{
    static unsigned int counter( 1U );
    double result( 0. );
   
    do
    {
        double x1 ( (double)rand() / (double)RAND_MAX );
        double x2 ( (double)rand() / (double)RAND_MAX );
       
        if ( ( counter % 2 ) == 0 )
            result = sqrt(-2. * log(x1)) * sin(2 * 3.14159265 * x2);
        else
            result = sqrt(-2. * log(x1)) * cos(2 * 3.14159265 * x2);
           
        result = mean + result * sd;
        counter++;
   
    } while ( ( result < minimum ) || ( result > maximum ) );
   
    return result;
}



void Generic2P2D_IMPES_Simulator::SeparateContinuousLineElementsToDifferentRegionAndAssignProperty()
{

    cout << "\nGeneric2P2D_IMPES_Simulator::SeparateContinuousLineElementsToDifferentRegionAndAssignProperty()\n";
    
    std::set<ControlVolumeElement<2U>* > visited;
    std::vector<std::vector<ControlVolumeElement<2U>* > > discovered_regions;
    
    
    for (std::vector<ControlVolumeElement<2U> >::iterator lcveitr = transportModel_->LineCVEsBegin();
         lcveitr != transportModel_->PointCVEsBegin(); ++lcveitr)
    {
        if (visited.find(&(*lcveitr)) == visited.end())
        {
            std::vector<ControlVolumeElement<2U>* > new_region;
            FindContinuousLineRegionContainsCVE(&(*lcveitr), visited, new_region);
            discovered_regions.push_back(new_region);
        }
        
    }

    for (std::vector<ControlVolumeElement<2U> >::iterator lcveitr = transportModel_->PerimeterLineCVEsBegin();
         lcveitr != transportModel_->PerimeterPointCVEsBegin(); ++lcveitr)
    {
        if (visited.find(&(*lcveitr)) == visited.end())
        {
            std::vector<ControlVolumeElement<2U>* > new_region;
            FindContinuousLineRegionContainsCVE(&(*lcveitr), visited, new_region);
            discovered_regions.push_back(new_region);
        }
        
    }
    
    Point<2U> xyz_min, xyz_max; 
    reservoirModel_->MinMaxCoordinates( xyz_min, xyz_max );

    double flank_ap(0.0005);
    double crest_ap(0.005);
    double min_ap(0.0001);
    double max_ap(0.001);
    double deviation(0.3);
    double crest_position(0.5 * (xyz_max[0] + xyz_min[0]));
   
    for (std::vector<std::vector<ControlVolumeElement<2U>* > >::iterator reg_itr = discovered_regions.begin();
         reg_itr != discovered_regions.end(); ++reg_itr)
    {
        
        double aperture;
        double perm;
        
        if ((((*reg_itr)[0U])->E()->BaryCenter())[0] < crest_position)
        {
            aperture = flank_ap + (crest_ap - flank_ap) / (crest_position - xyz_min[0]) * ((((*reg_itr)[0U])->E()->BaryCenter())[0] - xyz_min[0]);
        }
        else
        {
            aperture = crest_ap + (flank_ap - crest_ap) / (xyz_max[0] - crest_position) * ((((*reg_itr)[0U])->E()->BaryCenter())[0] - crest_position);
        }
        
        aperture = normalDistributionGenerator(aperture, aperture * deviation, min_ap, max_ap);
        perm = aperture * aperture / 12;
        
        for (std::vector<ControlVolumeElement<2U>* >::iterator subreg_itr = reg_itr->begin();
             subreg_itr != reg_itr->end(); ++subreg_itr)
        {
            (*subreg_itr)->E()->Store(propKeys_->k_key, ScalarVariable(PLAIN, perm));
            (*subreg_itr)->E()->Store(propKeys_->vm_key, ScalarVariable(PLAIN, aperture * 10));
        }
    
    }
    
    //transportModel_->CorrectVolumes("volume modifier");
    //transportModel_->CorrectSectorsVolume("volume modifier");

}


} // end csmp
