//
//  GenericFiniteVolumeTransport_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/01/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

// generic transport scheme
#include "GenericFiniteVolumeTransport_Test.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "FluxEvaluator.h"
#include "TimeStepEvaluator.h"
#include "ExplicitTransport.h"

// model
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"

// interfaces
#include "Standard_IO_Handler.h"
#include "ComputationalSettings.h"
#include "InputDataManager.h"
#include "ANSYS_Model3D.h"
#include "ModelTopology.h"
#include "VTK_Interface.h"

// integration od PDEs and post-processing
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"

// interrelations
#include "ConstantFactor.h"

// legacy finite-volume transport scheme
#include "StencilProcessor.h"
#include "ExplicitStencilProcessor.h"
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"
#include "NodeCenteredFiniteVolumeTransport.h"

#ifdef FV_STENCIL_TESTING
#include "FiniteVolumeStencil_refactored.h"
#else
#include "FiniteVolumeStencil.h"
#endif

using namespace std;


namespace csmp {

    namespace {
        const double s_internal_flux_rel_err = 1.0e-10;
        const double s_flux_balance_thresh = 1.0e-13;


        double rel_error(double x, double y)
        {
            auto mag = std::max(std::abs(x),std::abs(y));
            return std::abs(x - y) / mag;
        }
        
        double max_abs(double x, double y)
        {
            return std::max(std::abs(x),std::abs(y));
        }

        
        template<typename It, typename P>
        void
        jacobian_at_point(It& it, const P& p)
        {
            std::vector<double> DNR, DNS, DNT;
            (*it)->FE()->dNr( p[0], p[1], p[2], DNR );
            (*it)->FE()->dNs( p[0], p[1], p[2], DNS );
            (*it)->FE()->dNt( p[0], p[1], p[2], DNT );
            (*it)->CoordinateMatrix();
            (*it)->FE()->Jacobian( DNR, DNS, DNT );
        }
    }
    


// *************************************************************************************************
//
// definitions of auxiliary functions
//
// *************************************************************************************************

/**
    Tests performed:
    
    Speed comparison between projections made in parametric versus physical space while checking accuracy
    at same time.
*/
void GenericFiniteVolumeTransport_Test::run()
 {
    TestBasics();
    BenchmarkGlobalVersusParametricIntegration();
   
 } // end run





/**
    Testing finite volume projection calculations for the underpinning stencils.
    
    1. get model and compute pressure gradients
    
    2. verify (for the interior elements only) that the flux balance is indeed zero using established approach
    
    3. verify face by face computations
    
    4. verify flux balance cell-by-cell
    
    5. verify overall flux balances for hybrid element mesh
    
    6. compare computation times for flux balances
 
*/
void GenericFiniteVolumeTransport_Test::TestBasics()
 {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
     string  model_name("prism_test");
     // string  model_name("fracs4");

//      ANSYS_Model3D  model( model_name.c_str(), "CSMP-2phase-variables.txt");
//      ANSYS_Model3D  model( model_name.c_str(), "ExplicitTransport_Test-variables.txt");
      ANSYS_Model3D  model( model_name.c_str(), "VariableSet_TracerTransfer-variables.txt");
      printModelDimensions( model, true );

     // ------------------------------------------------------------
     // 2. configuring the model
     // ------------------------------------------------------------
      InputDataManager<3U>   model_configuration;
      ComputationalSettings  run_settings;
      model_configuration.ConfigureFromFile( model,
                                             model_name.c_str(),
                                             false, 
                                             true,   // 2) default prop.values
                                             true,   // 3) region prop.values
                                             true,   // 4) essential box-boundary conditions
                                             true,   // 5) essential flags
                                             true,   // 6) boundary conditions
                                             run_settings );
      Standard_IO_Handler  stdio;
      printRangeOfVariable( model, stdio, "permeability" );

     // ------------------------------------------------------------
     // 3. Building the transport scheme
     // ------------------------------------------------------------
     ExplicitTransport<3>  transport_scheme( model, "Model" );


      // optional visualization of the input permeability and boundary conditions
      VTK_Interface<3U>  vtk_output;
      vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );


     // -----------------------------------------------------------------------
     // 3. hydraulic conductivity computation
     // -----------------------------------------------------------------------
      const double  fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model.Apply( conductivity );
      printRangeOfVariable( model, "conductivity" );
      vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0 );


     // -----------------------------------------------------------------------
     // 4. computing a steady-state fluid pressure distribution in the model
     // -----------------------------------------------------------------------
      SteadyStateDiffusor<3U,Element> steady_state_pressure( model,
                                                            "conductivity", "fluid pressure",
                                                            "fluid volume source" );
    
      // postprocessing of pressure gradients and flow velocities
      VelocityAndVolumeFlux<3U>  postpro0( model, "conductivity", "porosity", "fluid pressure" );
      steady_state_pressure.AddPostProcess( &postpro0 );

      // the calculation of fluid pressure
      steady_state_pressure.ComputeSteadyState( model.Region("Model") );

      // results: the pore velocity is the Darcy velocity divided by the porosity
      printRangeOfVariable( model, stdio, "fluid pressure" );
      printRangeOfVariable( model, stdio, "velocity" );
      printRangeOfVariable( model, stdio, "pore velocity" );

      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
      vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1 );
   

     // -----------------------------------------------------------------------
     // 4. stepping over the model comparing facet by facet flux calculations
     // -----------------------------------------------------------------------
     const Region<3>& model_domain(model.Region("Model"));
     const csmp::Index p_key(model.Database().StorageKey("fluid pressure")),
                       K_key(model.Database().StorageKey("conductivity")),
                       v_key(model.Database().StorageKey("velocity"));
   
     std::vector<double> DNR, DNS, DNT;
     for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it )
       {
          const size_t nodes((*it)->Nodes());
          // 1. computing facet velocity in parametric space
          // -----------------------------------------------
          // 1.1 getting ipol-functions at barycentre and computing the pressure gradient in parametric space
          Point<3U> bctr = (*it)->FV()->Barycenter();
          (*it)->FE()->dNr( bctr[0], bctr[1], bctr[2], DNR );
          (*it)->FE()->dNs( bctr[0], bctr[1], bctr[2], DNS );
          (*it)->FE()->dNt( bctr[0], bctr[1], bctr[2], DNT );
          // pressure gradients / velocities
          const double K((*it)->Read(K_key));
          double  dpdr(0.), dpds(0.), dpdt(0.);
          Point<3U> vD(0.);
          for ( uint32_t i{0u}; i<nodes; ++i ) {
               const double p_node((*it)->N(i)->Read(p_key));
               dpdr   = DNR[i] * p_node;
               dpds   = DNS[i] * p_node;
               dpdt   = DNT[i] * p_node;
               vD[0] += -K * dpdr;
               vD[1] += -K * dpds;
               vD[2] += -K * dpdt;
            }
         
           double total_flux_physical = 0.0;
           double total_flux_parametric = 0.0;

          // 2. facet projections
          // --------------------
          const uint32_t facets=(*it)->Facets();
          for ( uint32_t i{0u}; i<facets; ++i ) {
              assert( (*it)->IsVolume() );
              
               // 2.1 classic way of calculating facet fluxes in physical space
               // ------------------------------------------------------------
              const double projected_velocity_physical = (*it)->ProjectionOnFacetNormal( i, v_key );
              // Darcy velocity computation
              const double flux_physical = (*it)->FacetArea(i) * projected_velocity_physical;
              total_flux_physical += flux_physical;
            
               // 2.2 parametric space computation
               // --------------------------------
               /* Requirements
                  - velocity does not depend on position in element. Only one Jacobian transformation is required
                  - mapping of facet area to physical space (sqrt(J^T J) ), non-square Jacobian squared by multiplication with its transpose
                  - mapping of pressure gradient to physical space
                  
                  Learnings
                  - the velocity (computed in parametric space), vD is transformed into physical space by pre-multiplication with JINV at barycentre
                  - this applies to all velocities projected onto facet normals in parametric space (this operation involves renormalisation 
                    of parametric space velocities to vD projected length
               */
              const Point<3u> parametric_facet_normal((*it)->ParametricFacetNormal(i));
              
              jacobian_at_point(it, (*it)->FacetPoint(i,0).Coordinates());
              
//              double jinvdet_ip = (*it)->FE()->JacobianInverse();
              const DenseMatrix<DM_MIN>& jinv_ip = (*it)->FE()->JINV;
              
              Point<3u> facet_normal_remapped(jinv_ip * parametric_facet_normal.Coordinates());
//              double fnrlen = facet_normal_remapped.Length();
              facet_normal_remapped.NormalizeLengthTo(1.0);
              Point<3u> vDlocal(jinv_ip * vD.Coordinates());
              double projected_velocity_parametric = dotProduct(facet_normal_remapped, vDlocal);
              double flux_parametric = projected_velocity_parametric * (*it)->FacetAreaMapped(i);
              
              total_flux_parametric += flux_parametric;
              
              // 3. testing that the fluxes are the same
              // ---------------------------------------

              
              std::cerr << std::setprecision(15);
              std::cerr << "p.v. physical   " << projected_velocity_physical << '\n';
              std::cerr << "p.v. parametric " << projected_velocity_parametric << '\n';
              std::cerr << "flux physical   " << flux_physical << '\n';
              std::cerr << "flux parametric " << flux_parametric << '\n';
              
              double relative_error = rel_error(flux_physical, flux_parametric);
              if (relative_error > 1e-10) {
#if 0
                  std::cerr << "Element type: " << parseFiniteElementType((*it)->FE_Type()) << '\n';
                  std::cerr << "Facet: " << i << ' ' << parseFacetType(((*it)->FV()->FacetType(i))) << '\n';
                  std::cerr << "Flux physical: " << flux_physical << '\n';
                  std::cerr << "Flux from parametric: " << flux_parametric << '\n';
                  std::cerr << "Rel error: 10^" << std::log(relative_error) << '\n';
#endif
              }
              else {
                  _equal( flux_physical, flux_parametric, 1e-10 );
              }
          }
           std::cerr << "Element type: " << parseFiniteElementType((*it)->FE_Type()) << '\n';
           std::cerr << "Total flux physical: " << total_flux_physical << '\n';
           std::cerr << "Total flux parametric: " << total_flux_parametric << '\n';
           std::cerr << "Total facets: " << (*it)->Facets() << '\n';
       }

 } // end TestBasics






/**
    Testing finite volume projection calculations for the underpinning stencils.
    
    1. get model and compute pressure gradients
    
    2. verify (for the interior elements only) that the flux balance is indeed zero using established approach
    
    3. verify face by face computations
    
    4. verify flux balance cell-by-cell
    
    5. verify overall flux balances for hybrid element mesh
    
    6. compare computation times for flux balances
 
*/
void GenericFiniteVolumeTransport_Test::BenchmarkGlobalVersusParametricIntegration()
 {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
     string  model_name("prism_test");
     // string  model_name("fracs4");

      ANSYS_Model3D  model( model_name.c_str(), "example25.txt");
      printModelDimensions( model, true );
      model.InstantiateFiniteVolumes();

     // ------------------------------------------------------------
     // 2. configuring the model
     // ------------------------------------------------------------
      InputDataManager<3U>   model_configuration;
      ComputationalSettings  run_settings;
      model_configuration.ConfigureFromFile( model,
                                             model_name.c_str(),
                                             false, 
                                             true,   // 2) default prop.values
                                             true,   // 3) region prop.values
                                             true,   // 4) essential box-boundary conditions
                                             true,   // 5) essential flags
                                             true,   // 6) boundary conditions
                                             run_settings );
      Standard_IO_Handler  stdio;
      printRangeOfVariable( model, stdio, "permeability" );

      // optional visualization of the input permeability and boundary conditions
      VTK_Interface<3U>  vtk_output;
      vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );


     // -----------------------------------------------------------------------
     // 3. hydraulic conductivity computation
     // -----------------------------------------------------------------------
      const double  fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model.Apply( conductivity );
      printRangeOfVariable( model, "conductivity" );
      vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0 );


     // -----------------------------------------------------------------------
     // 4. computing a steady-state fluid pressure distribution in the model
     // -----------------------------------------------------------------------
      SteadyStateDiffusor<3U,Element> steady_state_pressure( model,
                                                            "conductivity", "fluid pressure",
                                                            "fluid volume source" );
    
      // postprocessing of pressure gradients and flow velocities
      VelocityAndVolumeFlux<3U>  postpro0( model, "conductivity", "porosity", "fluid pressure" );
      steady_state_pressure.AddPostProcess( &postpro0 );

      // the calculation of fluid pressure
      steady_state_pressure.ComputeSteadyState( model.Region("Model") );

      // results: the pore velocity is the Darcy velocity divided by the porosity
      printRangeOfVariable( model, stdio, "fluid pressure" );
      printRangeOfVariable( model, stdio, "velocity" );
      printRangeOfVariable( model, stdio, "pore velocity" );

      vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
      vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1 );

      Region<3>& model_domain(model.Region("Model"));
      const csmp::Index p_key(model.Database().StorageKey("fluid pressure")),
                        K_key(model.Database().StorageKey("conductivity")),
                        v_key(model.Database().StorageKey("velocity"));
   
      std::vector<double> DNR, DNS, DNT;
      std::vector<double> cross_section(model_domain.Nodes());
      std::vector<double> velocity_magnitude(model_domain.Nodes());

     // -----------------------------------------------------------------------
     // 4. stepping over the model comparing facet by facet flux calculations
     // -----------------------------------------------------------------------
     for (auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it) {
         double csa = 0.0;
         double surface_area = 0.0;
         VectorVariable<3U> vc(PLAIN,PLAIN,PLAIN,0.,0.,0.);

         double vDmax(0.0);
         Point<3u> vDavg(0.0);

         for ( uint32_t j{0u}; j < (*it)->Parents(); ++j) {
             auto e = (*it)->Parent(j);
             e->Read(v_key, vc);
             if (vDmax < vc.P().Length()) {
                 vDavg = vc.P();
                 vDmax = vc.P().Length();
             }
         }
         velocity_magnitude[(*it)->Idx()] = vDmax;
         vDavg.NormalizeLengthTo(1.0);
         for ( uint32_t j{0u}; j < (*it)->Parents(); ++j) {
             auto e = (*it)->Parent(j);

             const auto child = (*it)->ParentNodeNumber(j);
             
             for ( uint32_t k{0u}; k < e->FV()->FacetsPerSector(child); ++k) {
                 const auto facet = e->FV()->FacetSurroundingSector(child, k);
                 double costheta = std::abs(dotProduct(e->FacetNormal(facet), vDavg));
                 csa += costheta * e->FacetArea(facet);
                 surface_area += std::abs(e->FacetArea(facet));
             }

         }
         cross_section[(*it)->Idx()] = csa * 0.5;
     }
 
     for ( auto pit = model_domain.PerimeterNodesBegin(); pit != model_domain.NodesEnd(); ++pit )
     {
         _test( (*pit)->Idx() < model_domain.Nodes() );
         cross_section[(*pit)->Idx()] = 0.;
         velocity_magnitude[(*pit)->Idx()] = 0.;
     }

     std::vector<double> flux_balance_parametric(model_domain.Nodes());
     std::vector<double> flux_balance_physical(model_domain.Nodes());
     std::vector<Point<3u>> directed_area_para(model_domain.Nodes());
     std::vector<Point<3u>> directed_area_phys(model_domain.Nodes());

     std::vector<std::pair<std::set<Point<3u>>,Element<3u>*> > elements;
     elements.reserve(model_domain.Cells());
     for (auto it = model_domain.CellsBegin(); it != model_domain.CellsEnd(); ++it) {
         std::set<Point<3u>> nodes;
         for (auto itn = (*it)->NodesBegin(); itn != (*it)->NodesEnd(); ++itn) {
             nodes.insert((*itn)->Coordinate());
         }
         elements.push_back(make_pair(nodes, *it));
     }
     std::sort(elements.begin(), elements.end());
     double maxtheta = 0;
     

     for (auto& element_key : elements)
       {
           auto iti = element_key.second;
           auto it = &iti;
       
           // std::cerr << "Element type: " << parseFiniteElementType((*it)->FE_Type()) << '\n';
           assert( (*it)->IsVolume() );


          const uint32_t nodes((*it)->Nodes());
          // 1. computing facet velocity in parametric space
          // -----------------------------------------------
          // 1.1 getting ipol-functions at barycentre and computing the pressure gradient in parametric space
           const Point<3U> bctr = (*it)->FV()->Barycenter();
           (*it)->FE()->dNr( bctr[0], bctr[1], bctr[2], DNR );
           (*it)->FE()->dNs( bctr[0], bctr[1], bctr[2], DNS );
           (*it)->FE()->dNt( bctr[0], bctr[1], bctr[2], DNT );
           // pressure gradients / velocities
           const double K((*it)->Read(K_key));
           Point<3U> vD(0.);
           for ( uint32_t i{0u}; i<nodes; ++i ) {
               const double p_node((*it)->N(i)->Read(p_key));
               vD[0] += -K * DNR[i] * p_node;
               vD[1] += -K * DNS[i] * p_node;
               vD[2] += -K * DNT[i] * p_node;
           }
//           double vDlength = vD.Length();
           
           jacobian_at_point(it, bctr.Coordinates());
           const DenseMatrix<DM_MIN> jac_bctr((*it)->FE()->JAC);
//           double jinvdet_bctr = (*it)->FE()->JacobianInverse();
           const DenseMatrix<DM_MIN> jinv_bctr((*it)->FE()->JINV);
           
           const Point<3U> vDproj(jinv_bctr * vD.Coordinates());

           // XXX check logic
           bool at_boundary = model_domain.IsPerimeterCell((*it)->Idx());

           // 2. facet projections
           // --------------------
           const uint32_t facets=(*it)->Facets();
           
//           auto fetype = (*it)->FE_Type();
           
           for ( uint32_t i{0u}; i<facets; ++i ) {
               assert( (*it)->IsVolume() );
               
               // 2.1 classic way of calculating facet fluxes in physical space
               // ------------------------------------------------------------
               const double projected_velocity_physical = (*it)->ProjectionOnFacetNormal( i, v_key );
               // Darcy velocity computation
               const double flux_physical = (*it)->FacetArea(i) * projected_velocity_physical;

               // 2.2 parametric space computation
               // --------------------------------
               /* Requirements
                - velocity does not depend on position in element. Only one Jacobian transformation is required
                - mapping of facet area to physical space (sqrt(J^T J) ), non-square Jacobian squared by multiplication with its transpose
                - mapping of pressure gradient to physical space
                
                Learnings
                - the velocity (computed in parametric space), vD is transformed into physical space by pre-multiplication with JINV at barycentre
                - this applies to all velocities projected onto facet normals in parametric space (this operation involves renormalisation
                of parametric space velocities to vD projected length
                */

               Point<3u> v0(0.0);
               Point<3u> v1(0.0);
               for ( uint32_t nn{0u}; nn < (*it)->Nodes(); ++nn) {
                   auto xform_weights = (*it)->FV()->FacetNormalTransformationNodeWeights(i, nn);
                   const Point<3u> n((*it)->N(nn)->Coordinate());
                   v0 += xform_weights.first * n;
                   v1 += xform_weights.second * n;
               }
               Point<3u> parametric_normal_remapped(crossProduct(v1,v0));
               parametric_normal_remapped.NormalizeLengthTo(1.0);
               double projected_velocity_parametric = dotProduct(parametric_normal_remapped, vDproj);
               double flux_parametric = projected_velocity_parametric * (*it)->FacetArea(i);

               uint32_t inside_node, outside_node;
               (*it)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
               
               _test( (*it)->N(inside_node)->Idx() < model_domain.Nodes() );
               _test( (*it)->N(outside_node)->Idx() < model_domain.Nodes() );

               const double scale = std::max(cross_section[inside_node], cross_section[outside_node])
                   * std::max(velocity_magnitude[inside_node], velocity_magnitude[outside_node]);

               const double abserr = std::abs(flux_parametric - flux_physical);
               const double relerr = rel_error(flux_parametric, flux_physical);
               std::cerr << "scale = " << scale << "\n";
               std::cerr << "abserr = " << abserr << "\n";
               std::cerr << "relerr = " << relerr << "\n";
               std::cerr << "abserr = 10^" << std::log10(abserr) << "\n";

               // Any discrepancy should be explainable by plain old numerical error,
               // or by integration error.
               // _test(abserr < std::max(1.0e-10, scale * 1.0e-10));
               
               flux_balance_parametric[(*it)->N(inside_node)->Idx()] += flux_parametric;
               flux_balance_parametric[(*it)->N(outside_node)->Idx()] -= flux_parametric;
               flux_balance_physical[(*it)->N(inside_node)->Idx()] += flux_physical;
               flux_balance_physical[(*it)->N(outside_node)->Idx()] -= flux_physical;
               directed_area_para[(*it)->N(inside_node)->Idx()] += parametric_normal_remapped * (*it)->FacetArea(i);
               directed_area_para[(*it)->N(outside_node)->Idx()] -= parametric_normal_remapped * (*it)->FacetArea(i);
               directed_area_phys[(*it)->N(inside_node)->Idx()] += (*it)->FacetNormal(i) * (*it)->FacetArea(i);
               directed_area_phys[(*it)->N(outside_node)->Idx()] -= (*it)->FacetNormal(i) * (*it)->FacetArea(i);

               // 3. testing that the fluxes are the same
               // ---------------------------------------
               if (!at_boundary) {
                  _equal( flux_physical, flux_parametric, s_internal_flux_rel_err );
               }
           }
       }
     
     std::cerr << "maxtheta = " << maxtheta << '\n';
     
     for ( auto pit = model_domain.PerimeterNodesBegin(); pit != model_domain.NodesEnd(); ++pit )
     {
         _test( (*pit)->Idx() < model_domain.Nodes() );
         flux_balance_physical[(*pit)->Idx()] = 0.;
         flux_balance_parametric[(*pit)->Idx()] = 0.;
         directed_area_phys[(*pit)->Idx()] = 0.;
         directed_area_para[(*pit)->Idx()] = 0.;
     }
     
     double fmin_phys = +std::numeric_limits<double>::max();
     double fmax_phys = -std::numeric_limits<double>::max();
     double fmin_para = +std::numeric_limits<double>::max();
     double fmax_para = -std::numeric_limits<double>::max();

     size_t weird_nodes = 0;
     double maxdelta = 0;
     for ( uint32_t i{0u}; i < flux_balance_physical.size(); ++i) {
         _test(directed_area_para[i].Length() < 1.0e-14);
         _test(directed_area_phys[i].Length() < 1.0e-14);

         double scale = cross_section[i] * velocity_magnitude[i];
         double phys = flux_balance_physical[i];
         double para = flux_balance_parametric[i];

         fmin_phys = std::min(fmin_phys, phys);
         fmax_phys = std::max(fmin_phys, phys);
         fmin_para = std::min(fmin_para, para);
         fmax_para = std::max(fmin_para, para);
         double delta = std::abs(phys - para);
         std::cerr << "scale[" << i << "] = " << scale << '\n';
         std::cerr << "flux_balance_phys[" << i << "] = " << phys << '\n';
         std::cerr << "flux_balance_para[" << i << "] = " << para << '\n';
         std::cerr << "flux_balance_delta[" << i << "] = " << delta << '\n';
         maxdelta = std::max(maxdelta, delta);
         if (delta > s_flux_balance_thresh) {
             ++weird_nodes;
             std::cerr << "For node " << i << "\n";
             std::cerr << "Flux balance phys " << phys << "\n";
             std::cerr << "Flux balance para " << para << "\n";
             std::cerr << "delta " << delta << " = 10^" << std::log10(delta) << "\n";
         }
         // _test(std::abs(flux_balance_physical[i] - flux_balance_parametric[i]) < s_flux_balance_thresh);
     }
     std::cerr << "Flux balance phys: (" << fmin_phys << ", "  << fmax_phys << ")\n";
     std::cerr << "Flux balance para: (" << fmin_para << ", " << fmax_para << ")\n";
     std::cerr << "maxdelta: " << maxdelta << '\n';

     std::cerr << "Weird nodes: " << weird_nodes << " / " << model_domain.Nodes() << '\n';
 } // end





/**
    Executing and testing all finite-volume related methods
*/
static bool test_NCFVT_methods( Model<3U>& model3D, NodeCenteredFiniteVolumeTransport<3U>& advector3D, VTK_Interface<3U>& vtk_output )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  //----------------------------------------------------------------------
  //------------------------- TESTING METHODS ----------------------------
  //----------------------------------------------------------------------

          /*
          * -- TESTING OF NCFVT METHODS-----------------------
          *
          * -- 1.)AdvectVariable -----------------------------
          */
          const double timeInterval(1.e3);
          const double courantMultiplier(1.e5);
          double courantIncrement;
          cout << "\n\n\n\n\n";
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 1.)AdvectVariable -----------------------------\n";
          cout << "\nInputAguments - timeInterval: " << timeInterval <<
                  ", courantMultiplier: " << courantMultiplier;

          courantIncrement = advector3D.AdvectVariable(timeInterval,courantMultiplier);
          //                            ^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;



          /*
          * -----------   2.)TransportPhase    -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 2.)TransportPhase 1D ONLY ---------------------\n";
          //cout << "\nInputAguments - timeInterval: " << timeInterval;

          //advector1D.TransportPhase(relperms, timeInterval);
          //         ^^^^^^^^^^^^^^




          /*
          * -----------  3.)CourantIncrement() -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 3.)CourantIncrement() 1D ONLY -----------------\n";

          //courantIncrement = advector1D.CourantIncrement();
          //                            ^^^^^^^^^^^^^^^^

          //cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * -----------   4.)CourantIncrement    -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 4.)CourantIncrement 1D ONLY -------------------\n";

          //courantIncrement = advector1D.CourantIncrement(relperms);
          //                            ^^^^^^^^^^^^^^^^

          //cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * ----------- 5.)AnisotropicCourantIncrement() -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 5.)AnisotropicCourantIncrement() --------------\n";

          courantIncrement = advector3D.AnisotropicCourantIncrement();
          //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * ----------- 6.)AnisotropicCourantIncrement -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 6.)AnisotropicCourantIncrement ----------------\n";

          //courantIncrement = advector3D.AnisotropicCourantIncrement(relperms);
          //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nReturns courantIncrement of: " << courantIncrement;


          /*
          * -----------    7.)CFL_Multiplier   -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 7.)CFL_Multiplier -----------------------------\n";

          double CFLcheck(1.e2);
          advector3D.CFL_Multiplier(CFLcheck);
          //         ^^^^^^^^^^^^^^
          assert(CFLcheck == advector3D.CFL_Multiplier());
          cout << "\nManual CFL input value: " << CFLcheck ;
          cout << "\nReturn function for CFL multiplier gives: " << advector3D.CFL_Multiplier();
          //                                                                   ^^^^^^^^^^^^^^


          /*
          * -----------    8.)Inflow/Outflow   -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 8.)Model Inflow/Outflow -----------------------\n";

          double outflow(advector3D.ModelOutflow());
          //                          ^^^^^^^^^^^^
          double inflow(advector3D.ModelInflow());
          //                         ^^^^^^^^^^^^

          cout << "\nadvector.ModelOutflow(): " << outflow ;
          cout << "\nadvector.ModelInflow(): " << inflow ;


          /*
          * -----------  9.)BoundaryFluxes  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 9.)BoundaryFluxes -----------------------------\n";

          double boundaryFluxes(advector3D.BoundaryFluxes(inflow, outflow));
          //                                 ^^^^^^^^^^^^^^
          cout << "\nBoundary Fluxes with previous as Input";
          cout << " arguments returns: " << boundaryFluxes ;


          /*
          * -----------  10.)FluxBalance  -----------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 10.)FluxBalance -------------------------------\n";

          double fmin, fmax;
          advector3D.FluxBalance(fmin, fmax);
          //         ^^^^^^^^^^^
          cout << "\nFluxBalance returns " << fmin << " as minimum and ";
          cout << fmax << " as maximum FV flux balance.";


          /*
          * ----- 11.)MultiplyScalarNodePropertyByFiniteVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-----------------------";
          cout << "\n\n/** -- 11.)MultiplyScalarNodePropertyByFiniteVolume --\n";

          advector3D.MultiplyScalarNodePropertyByFiniteVolume("fluid pressure");
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


          /*
          * ----- 12.)VolumeIntegrateScalarFiniteElementVariable  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS------------------------";
          cout << "\n\n/** -- 12.)VolumeIntegrateScalarFiniteVolumeVariable --\n";

          double poreVolume(0.);
          poreVolume =
          advector3D.VolumeIntegrateScalarFiniteElementVariable("porosity", true);
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          cout << "\nIntegration of porosity  returns " << poreVolume << " for FE Integration.";



          /*
          * ----- 13.)FiniteVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 13.)FiniteVolumes -------------------------------\n";

          poreVolume =
          advector3D.FiniteVolume("concentration");
          //         ^^^^^^^^^^^^
          vtk_output.OutputDataToVTK( model3D, "finite-volume", "concentration", 0 );

          cout << "\nFinite Volumes of porosity  returns " << poreVolume;

          /*
          * ----- 14.)VolumeIntegrate  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 14.)VolumeIntegrate -----------------------------\n";

          advector3D.VolumeIntegrate("porosity", "poreVolume");
          //         ^^^^^^^^^^^^^^^
          advector3D.VolumeIntegrate("porosity", "saturation oil", "oilVolume");
          //         ^^^^^^^^^^^^^^^

          /*
          * ----- 15.)AssignScalarBoundaryValues  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 15.)AssignScalarBoundaryValues ------------------\n";
          double customPressure(2.e7);
          advector3D.AssignScalarBoundaryValues(LEFT, "fluid pressure", DIRICH, customPressure, true);
          //         ^^^^^^^^^^^^^^^^^^^^^^^^^^


          /*
          * ----- 16.)PoreVolume  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 16.)PoreVolume -----------------------------\n";

          poreVolume =
          advector3D.PoreVolume(50);
          //         ^^^^^^^^^^
          cout << "\nPoreVolume returns " << poreVolume;


          /*
          * ----- 17.)Out  ------
          */
          cout << "\n\n\n\n\n/** -- TESTING OF NCFVT METHODS-------------------------";
          cout << "\n\n/** -- 17.)Out -----------------------------------------\n";

          advector3D.Out();
          //         ^^^

  return true;

} // test_NCFVT_methods





static void testNodeCenteredFiniteVolumeStencils( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  //Testing NodeCenteredFiniteVolumeTransport.h
  //method testFiniteVolumeStencil

  cout << "Testing NodeCenteredFiniteVolumeTransport.h method testFiniteVolumeStencil " << endl;
  testFiniteVolumeStencil( sg.Database(), sg.Region( "Model" ), advector );

  test_NCFVT_methods( sg, advector, vtkOut );
  cout << "End of Testing!" << endl;
  
} //end TestNodeCenteredFiniteVolumeStencil






/**
  Calculates flux mismatch for predefined velocity field and normalizes it
  by finite volume. The finite volume and absolute value of the flux mismatch
  are reported to the variables "finite volume" and "nodal flux mismatch", respectively.
  do not use when surface elements are also present in model!
*/
static double  testNodeCenteredFiniteVolumeTransport_PrescribedVelocity( Model<3U>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  VectorVariable<3U>  velo(PLAIN,PLAIN,PLAIN, 3., 7., 1. );
  velo /= velo.Length(); // unit length
  velo.Out();
  sg.InputPropertyValue( "velocity", velo );

  cout <<"\n\tMeasuring the time required to build basic transport algorithm."<< endl;
  clock_t ticks = clock();
  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

  PropertyHandle<3U>  fv( sg,"finite volume",SCALAR,NODE);
  advector.FiniteVolume( "finite volume" );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Measuring the divergence of fluxes."<< endl;
  advector.Divergence( "velocity", "nodal flux mismatch" );

  // identifying the Dirichlet boundaries (since they will have in or outflow)
  csmp::Index  pf_key = sg.Database().StorageKey("fluid pressure");

  // zapping result values at model boundaries and normalizing divergence by finite volume
  csmp::Index          fv_key   = sg.Database().StorageKey("finite volume");
  csmp::Index          prop_key = sg.Database().StorageKey("nodal flux mismatch");
  ScalarVariable       sc;
  double             emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  // for all interior nodes we calculate the normalised flux balance
  for ( auto it=gref.NodesBegin(); it!=gref.PerimeterNodesBegin(); it++ ) {
       sc = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       (*it)->Store( prop_key, sc );
       emax = std::max( emax, fabs(sc()) );
    }

  // for all boundary nodes we set the balance to zero because we cannot evaluate it
  for ( auto it=gref.PerimeterNodesBegin(); it!=gref.NodesEnd(); it++ )
    (*it)->Store( prop_key, sc=0. );

  // finding the worst finite volume and analyzing it
  for ( auto it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - fabs((*it)->Read( prop_key ))) <= numeric_limits<double>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( auto i{0}; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

  return emax;

 } // end testNodeCenteredFiniteVolumeTransport_PrescribedVelocity










static void advectVariableExplicit( Model<3U>& sg, bool second_order )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  ExplicitNodeCenteredFiniteVolumeTransport<3U,ExplicitStencilProcessor>  explicit_advector(
                                                                            "Model", sg,
                                                                            "porosity",
                                                                            "concentration",
                                                                            "velocity",
                                                                            "nodal fluid volume source",
                                                                             second_order );

   cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
   if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
   else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
   cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement() << endl;
   cout <<"\nEnter advection time: ";
   double time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
   clock_t ticks = clock();
   explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
   ticks = clock() - ticks;
   cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit



static void advectVariableExplicit( Model<3U>& sg, const char* region, bool second_order )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  ExplicitNodeCenteredFiniteVolumeTransport<3U,ExplicitStencilProcessor>  explicit_advector(
                                                                            region, sg,
                                                                            "porosity",
                                                                            "concentration",
                                                                            "velocity",
                                                                            "nodal fluid volume source",
                                                                             second_order );

   cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
   if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
   else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
   cout <<"\nThe grid Courant number is "<< explicit_advector.AnisotropicCourantIncrement() << endl;
   cout <<"\nEnter advection time: ";
   double time_interval;
   cin >> time_interval;

   cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock();
  explicit_advector.AdvectVariable( time_interval, 0.1, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

} // end advectVariableExplicit








static void advectVariableFirstOrderImplicit( Model<3U>& sg, VTK_Interface<3U>& vtkOut )
 {
    // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
    //ostream &cout = *GetStream();

    NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source",false, false );

    cout <<"\n\nadvectVariableFirstOrderImplicit: Configuring TRANSPORT simulation: IMPIMS"<< endl;
    cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement();
    cout.flush();
    cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
    double time_interval, Courant_multiplier;
    cin >> time_interval >> Courant_multiplier;
    cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
    clock_t ticks = clock();
    for ( int i=0; i<5; ++i ) {
        advector.AdvectVariable( time_interval/5, Courant_multiplier );
        vtkOut.OutputDataToVTK( sg, "concentration", "concentration", i+1 );
     }
    ticks = clock() - ticks;
    cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableFirstOrderImplicit





// restricted to a group
static void advectVariableFirstOrderImplicit( Model<3U>& sg, const char* group )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( group, sg,
                                                  "porosity", "concentration", "velocity",
                                                  "nodal fluid volume source", false, false );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Configuring TRANSPORT simulation: IMPIMS for region'"<< group <<"'"<< endl;
  cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;

  cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
  double time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //         fluxbalancecorrection=true, updateporevols=false
  advector.AdvectVariable( time_interval, Courant_multiplier, false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableFirstOrderImplicit





static void advectVariableSecondOrderImplicit( Model<3U>& sg, bool bijective_mapping )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source", true, false );

  cout <<"\n\nadvectVariableSecondOrderImplicit: Configuring TRANSPORT simulation: ";
  if ( bijective_mapping ) cout <<" IMPIMS with BIJECTIVE MAPPING."<< endl;
  else                     cout <<" IMPIMS without BIJECTIVE MAPPING."<< endl;
  cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;
  cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
  double time_interval, Courant_multiplier;
  cin >> time_interval >> Courant_multiplier;

  cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
  clock_t ticks = clock(); //                  flux_balance_correction  update_pore_volumes
  advector.AdvectVariable( time_interval, Courant_multiplier, true, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

 } // end advectVariableSecondOrderImplicit




static  void advectVariableSecondOrderImplicitSecondOrderInTime( Model<3U>& sg, bool bijective_mapping )
   {
     // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
     //ostream &cout = *GetStream();

    NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                     "porosity", "concentration", "velocity",
                                                     "nodal fluid volume source", true, true );

    cout <<"\n\nadvectVariableSecondOrderImplicitSecondOrderInTime: Configuring TRANSPORT simulation: ";
    if ( bijective_mapping ) cout <<" IMPIMS with BIJECTIVE MAPPING."<< endl;
    else                     cout <<" IMPIMS without BIJECTIVE MAPPING."<< endl;
    cout <<"\nThe grid Courant number is "<< advector.AnisotropicCourantIncrement() << endl;
    cout <<"\nEnter advection time deduced from flow velocity and model-X extent (in seconds) and Courant multiplier: ";
    double time_interval, Courant_multiplier;
    cin >> time_interval >> Courant_multiplier;

    cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
    clock_t ticks = clock();
    advector.AdvectVariable( time_interval, Courant_multiplier, bijective_mapping );
    ticks = clock() - ticks;
    cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

   } // end advectVariableSecondOrderImplicit









static void testNodeCenteredFiniteVolumeTransport( Model<3U>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  cout <<"\n\tMeasuring the time required to build basic transport algorithm."<< endl;
  clock_t ticks = clock();
  NodeCenteredFiniteVolumeTransport<3U>  advector( "Model", sg,
                                                   "porosity", "concentration", "velocity",
                                                   "nodal fluid volume source",false, false );
  ticks = clock() - ticks;
  cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

  advector.FiniteVolume( "finite volume" );

  cout <<"\n\nadvectVariableFirstOrderImplicit: Measuring the divergence of fluxes."<< endl;
  advector.Divergence( "velocity", "nodal flux mismatch" );

  // identifying the Dirichlet boundaries (since they will have in or outflow)
  csmp::Index  pf_key = sg.Database().StorageKey("fluid pressure");

  // zapping result values at model boundaries and normalizing divergence by finite volume
  csmp::Index     fv_key   = sg.Database().StorageKey("finite volume");
  csmp::Index     prop_key = sg.Database().StorageKey("nodal flux mismatch");
  ScalarVariable  sc;
  double        emax(0.);
  Region<3>&  gref(sg.Region("Model"));

  for ( auto it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ ) {
       sc = fabs((*it)->Read( prop_key ) / (*it)->Read( fv_key ));
       if ( (*it)->AtBoundary() != NOT and (*it)->Status( pf_key ) == DIRICH )
         (*it)->Store( prop_key, sc=0. );
       else
         (*it)->Store( prop_key, sc );
       emax = std::max( emax, sc() );
    }

  // finding the worst finite volume and analyzing it
  for ( auto it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
    if ( fabs(emax - (*it)->Read( prop_key )) <= numeric_limits<double>::epsilon() ) {
         cout <<"\ntestNodeCenteredFiniteVolumeTransport: worst finite volume: "<< endl;
         (*it)->Out();
         cout <<"\ncomposed of the element types: "<< endl;
         for ( auto i{0}; i<(*it)->Parents(); i++ )
           cout << parseFiniteElementType( (*it)->Parent(i)->FE()->ElementType() ) << endl;
         cout << endl << endl;
      }

 } // end TestNodeCenteredFiniteVolumeTransport












// *************************************************************************************************
//
// definitions of main test function:  run()
//
// *************************************************************************************************


static void testSchemeAsComponent()
  {
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
      string  model_name("prism_test");
      //cout <<"\nmain: Enter name of 'ANSYS TETRA' input file (binary): ";
      //cin >> model_name;

      ANSYS_Model3D  model3D( model_name.c_str(), "example25.txt");
      printModelDimensions( model3D, true );


     // ------------------------------------------------------------
     // 2. configuring the model
     // ------------------------------------------------------------
      InputDataManager<3U>  model_configuration;
      ComputationalSettings  run_settings;
      model_configuration.ConfigureFromFile( model3D,
                                             model_name.c_str(),
                                             false, 
                                             true,   // 2) default prop.values
                                             true,   // 3) group prop.values
                                             true,   // 4) essential box-boundary conditions
                                             true,   // 5) essential flags
                                             true,   // 6) boundary conditions
                                             run_settings );
      Standard_IO_Handler  stdio;
      printRangeOfVariable( model3D, stdio, "permeability" );

      // optional visualization of the input permeability and boundary conditions
      VTK_Interface<3U>  vtk_output;
      vtk_output.OutputDataToVTK( model3D, "permeability", "permeability", 0 );
      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 0 );


     // -----------------------------------------------------------------------
     // 3. hydraulic conductivity computation
     // -----------------------------------------------------------------------
      const double  fluid_viscosity(1.0e-03);
      ConstantFactor<3U,divides>  conductivity( model3D.Database(),
                                               "conductivity", "permeability",
                                                fluid_viscosity );
      model3D.Apply( conductivity );
      printRangeOfVariable( model3D, "conductivity" );

      vtk_output.OutputDataToVTK( model3D, "conductivity", "conductivity", 0 );


     // -----------------------------------------------------------------------
     // 4. computing a steady-state fluid pressure distribution in the model
     // -----------------------------------------------------------------------
      SteadyStateDiffusor<3U,Element> steady_state_pressure( model3D,
                                                            "conductivity", "fluid pressure",
                                                            "fluid volume source" );
    
      // postprocessing of pressure gradients and flow velocities
      VelocityAndVolumeFlux<3U>  postpro0( model3D, "conductivity", "porosity", "fluid pressure" );
      steady_state_pressure.AddPostProcess( &postpro0 );

      // the calculation of fluid pressure
      steady_state_pressure.ComputeSteadyState( model3D.Region("Model") );

      // results: the pore velocity is the Darcy velocity divided by the porosity
      printRangeOfVariable( model3D, stdio, "fluid pressure" );
      printRangeOfVariable( model3D, stdio, "velocity" );
      printRangeOfVariable( model3D, stdio, "pore velocity" );

      vtk_output.OutputDataToVTK( model3D, "fluid-pressure", "fluid pressure", 1 );
      vtk_output.OutputDataToVTK( model3D, "velocity",       "velocity",       1 );


     // ---------------------------------------------------------------------------------
     // 5. Testing the generic transport scheme for tracer transport in single phase flow
     // ---------------------------------------------------------------------------------
     //   (here you can compare different schemes with one another and
     //    overstep CFL to see how this adds numerical diffusion to the solution)
     /*  
          Options
          -------
          A. prescribed divergence free velocity field
            1=explicit,
            2=explicit, O(2),
            3=implicit, 
            4=implicit O(2), 
            5=4+bijective mapping, 
            6,7=tests of volume integration,
            8=TestNodeCenteredFiniteVolumeStencil
      
          B. computed velocity field
     */
     // -----------------------------------------------------------------------

      int32_t     tmethod(1);
      double  max_error(1.0e-5);
    
      // 5.1 setting up the transport scheme
      /*
           We store the volume of the finite volumes and their pore volumes
           
           - finite volume
           - finite volume (effective) pore volume
           - sector volume (stored at sector integration point)
           - sector weight: sector pore volume / finite volume pore volume = weighting factor
      */
      Region<3U>&  flow_domain(model3D.Region("Model"));
      const bool initialize_flux(true);
      initializeFiniteVolumeProperties( model3D, flow_domain, initialize_flux );
    
// TESTING SECTOR INTEGRATION POINT STORAGE
      const csmp::Index swt_key(model3D.Database().StorageKey("node number"));
      // sector storage: writing global node numbers to sector IP's and reading them out
      for ( auto it=flow_domain.CellsBegin(); it!=flow_domain.CellsEnd(); ++it )
        for ( auto i{0}; i<(*it)->Sectors(); ++i )
          (*it)->Store( i, 0U, swt_key, makeScalar(PLAIN,(*it)->N(i)->Idx()) );
        
      // reading out node numbers and their double equivalents stored at the sector integration points
      for ( auto it=flow_domain.CellsBegin(); it!=flow_domain.CellsEnd(); ++it ) {
           cerr <<"\nelement: "<< (*it)->Idx() << endl;
           for ( auto i{0}; i<(*it)->Sectors(); ++i ) {
                cerr << (*it)->N(i)->Idx() <<":";
                cerr << (*it)->Read( i, 0U, swt_key ) <<" ";
             }
        }

//  TESTING FV VOLUME CALCULATIONS
      const csmp::Index fv_key(model3D.Database().StorageKey("finite element volume"));
      const csmp::Index sv_key(model3D.Database().StorageKey("fv sector volume"));
      const csmp::Index fvphi_key(model3D.Database().StorageKey("fv pore volume"));
      // accumulating matching sector volumes with finite element volumes
      double volume(0.);
      for ( auto it=flow_domain.NodesBegin(); it!=flow_domain.NodesEnd(); ++it )
        volume += (*it)->Read( fv_key );
      cerr <<"\nFV total volume: "<< volume;

      volume = 0.;
      for ( auto it=flow_domain.CellsBegin(); it!=flow_domain.CellsEnd(); ++it )
        for ( auto i{0}; i<(*it)->Sectors(); ++i )
          volume += (*it)->Read( i, 0U, sv_key );
      cerr <<"\nFV total volume: "<< volume;

      double pvolume(0.);
      for ( auto it=flow_domain.NodesBegin(); it!=flow_domain.NodesEnd(); ++it )
        pvolume += (*it)->Read( fvphi_key );
      cerr <<"\nFV total volume: "<< pvolume;



      // first order version only
      ExplicitTransport<3U>  explicit_advector( model3D, "Model" );
      printRangeOfVariable( model3D, "sector volume" );
      printRangeOfVariable( model3D, "sector pore volume" );
      printRangeOfVariable( model3D, "finite volume" );
      printRangeOfVariable( model3D, "FV pore volume" );

      cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
    //  if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
      cout <<" FIRST ORDER SCHEME."<< endl;
      cout <<"\nThe grid Courant number is "<< explicit_advector.TimeIncrement() << endl;
      cout <<"\nEnter advection time: ";
      double time_interval;
      cin >> time_interval;

      cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
      clock_t ticks = clock();
      explicit_advector.AdvectVariable( time_interval );
      ticks = clock() - ticks;
      cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;

      switch( tmethod ) {
           case 1:
              cout <<"\nmain: Would you like to restrict computation to group (yes=1, 0=no)? ";
              cin >> tmethod;
              if ( tmethod != 1 ) advectVariableExplicit( model3D, false /* second order=false */ );
              else {
                   string group_name;
                   cout <<"\nmain: Enter name of model region: ";
                   cin >> group_name;
                   advectVariableExplicit( model3D, group_name.c_str() );
                   vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "new-concentration", "new concentration", 1, true );
                }
             break;
           case 2:  advectVariableExplicit( model3D, true );
             break;
           case 3:
              cout <<"\nmain: Would you like to restrict computation to model region (yes=1, 0=no)? ";
              cin >> tmethod;
              if ( tmethod != 1 ) advectVariableFirstOrderImplicit( model3D, vtk_output );
              else {
                   string group_name;
                   cout <<"\nmain: Enter name of region: ";
                   cin >> group_name;
                   advectVariableFirstOrderImplicit( model3D, group_name.c_str() );
                   vtk_output.OutputDataToVTK( model3D, group_name.c_str(), "concentration", "concentration", 1, true );
                }
             break;
           case 4:  advectVariableSecondOrderImplicit( model3D, false );
             break;
           case 5:  advectVariableSecondOrderImplicit( model3D, true );
             break;
           case 6:
               cout <<"\nmain: Calculated flux mismatch: ";
               cout << testNodeCenteredFiniteVolumeTransport_PrescribedVelocity( model3D ) << endl;
               printRangeOfVariable( model3D, stdio, "finite volume" );
               vtk_output.OutputDataToVTK( model3D, "finite-volume", "finite volume", 1 );

               max_error = printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
               vtk_output.OutputDataToVTK( model3D, "nodal-flux-mismatch", "nodal flux mismatch", 1 );

               if ( max_error > 1.0e-7 ) {
                    model3D.FormRegionFrom( "corrupted-flux", "nodal flux mismatch", 1e-7, 100. );
                    vtk_output.OutputDataToVTK( model3D, "corrupted-flux", "nodal-flux-mismatch", "nodal flux mismatch", 1, true );
                }
             break;
           case 7:
              testNodeCenteredFiniteVolumeTransport( model3D );
              printRangeOfVariable( model3D, stdio, "finite volume" );
              vtk_output.OutputDataToVTK( model3D, "finite-volume", "finite volume", 1 );

              max_error = printRangeOfVariable( model3D, stdio, "nodal flux mismatch" );
              vtk_output.OutputDataToVTK( model3D, "nodal-flux-mismatch", "nodal flux mismatch", 1 );

              if ( max_error > 1.0e-7 ) {
                   model3D.FormRegionFrom( "corrupted-flux", "nodal flux mismatch", 1e-7, 100. );
                   vtk_output.OutputDataToVTK( model3D, "corrupted-flux", "nodal-flux-mismatch", "nodal flux mismatch", 1, true );
               }
             break;

           case 8:
             break;

           default:
               cout <<"\nmain: Transport method not recognized."<< endl;
             return;
        }

      printRangeOfVariable( model3D, stdio, "concentration" );
      vtk_output.OutputDataToVTK( model3D, "concentration", "concentration", 99 );

      cout <<"\nmain: That's it."<< endl;

  } // end run



} // end csmp
