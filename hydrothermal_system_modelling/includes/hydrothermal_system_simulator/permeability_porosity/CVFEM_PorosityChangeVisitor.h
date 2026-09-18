// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_POROSITY_CHANGE_VISITOR_H
#define CVFEM_POROSITY_CHANGE_VISITOR_H

/*   Changelog
     November 2014, Philipp Weis:
     - initial port to CSMP++ and strong simplification as compared to the csp5-version.
*/

#include "Visitor.h"
#include "Model.h"

#include "Exception.h"

namespace csmp
{

  /// Calculating the volume and pore volume of a node-centered control volume.

  template<size_t dim>
  class CVFEM_PorosityChangeVisitor : public Visitor<dim>
  {
    public:
      CVFEM_PorosityChangeVisitor( Model<dim>
                                   & ); // nodal variable to be used in further calculations and post-processing
      ~CVFEM_PorosityChangeVisitor();

      void ReadInputVariables( Model<dim> &model,
                              bool ext_with_gold, bool ext_with_quartz, bool ext_with_tracer,
                              bool ext_with_lithium, bool ext_with_magma);

      virtual void Visit(Node<dim> *n);     // Application level and target are set to NODE in the constructor.
      virtual void Visit(Model<dim> *m);    // Application level and target are set to NODE in the constructor.
      virtual void Visit(Region<dim> *r);   // Application level and target are set to NODE in the constructor.

    private:

      csmp::Index pore_volume_change_factor_key, pore_vol_key, phi_key, rr_key, cpr_key,
           cpf_key, rho_bulk_key, betaf_key, betar_key, sh_key;
      csmp::Index ncp_key, nQ_key, CT_key, mt_key, msp_key, Htp_key;
      csmp::Index mml_key, mmv_key, mmld_key, mmvd_key, eml_key, emv_key, emld_key, emvd_key,
           xml_key, xmv_key;
      csmp::Index ml_key, mv_key, mf_key, hCl_key, hCv_key, xCl_key, xCv_key;
      csmp::Index m_mml_key, m_mmv_key, m_xml_key, m_xmv_key,
           cml_key, cmv_key,
           znml_key, znmv_key,
           so2ml_key, so2mv_key,
           auhs0ml_key, auhs0mv_key,
           sio2ml_key, sio2mv_key,
           tracerml_key, tracermv_key,
           liml_key, limv_key;
      csmp::Index m_ml_key, m_mv_key, m_m_key, m_xCl_key, m_xCv_key, m_ms_key,
           cCl_key, cCv_key, cCf_key,
           znCl_key, znCv_key, znCf_key,
           so2Cl_key, so2Cv_key, so2Cf_key,
           auhs0Cl_key, auhs0Cv_key, auhs0Cf_key,
           sio2Cl_key, sio2Cv_key, sio2Cf_key,
           tracerCl_key, tracerCv_key, tracerCf_key,
           liCl_key, liCv_key, liCf_key;
      csmp::Index rl_transport_key, rv_transport_key;

      ScalarVariable pore_volume_change_factor, pore_vol, phi, rr, cpr, cpf, rho_bulk, betaf,
                     betar, sh;
      ScalarVariable ncp, nQ, CT, mt, msp, Htp;
      ScalarVariable mml, mmv, mmld, mmvd, eml, emv, emld, emvd, xml, xmv;
      ScalarVariable ml, mv, mf, hCl, hCv, xCl, xCv;
      ScalarVariable m_mml, m_mmv, m_xml, m_xmv,
                     cml, cmv,
                     znml, znmv,
                     h2ml, h2mv,
                     co2ml, co2mv,
                     so2ml, so2mv,
                     auhs0ml, auhs0mv,
                     sio2ml, sio2mv,
                     tracerml, tracermv,
                     liml, limv;
      ScalarVariable m_ml, m_mv, m_m, m_xCl, m_xCv, m_ms,
                     cCl, cCv, cCf,
                     znCl, znCv, znCf,
                     h2Cl, h2Cv, h2Cf,
                     co2Cl, co2Cv, co2Cf,
                     so2Cl, so2Cv, so2Cf,
                     auhs0Cl, auhs0Cv, auhs0Cf,
                     sio2Cl, sio2Cv, sio2Cf,
                     tracerCl, tracerCv, tracerCf,
                     liCl, liCv, liCf;
      ScalarVariable rl_transport, rv_transport;

      bool           with_gold, with_quartz, with_lithium, with_tracer, with_magma;


  };

  /**
       @class CVFEM_PorosityChangeVisitor CVFEM_PorosityChangeVisitor.h

       @author Philipp Weis, ETH Zuerich
       @section contact Contact
       philipp.weis@erdw.ethz.ch

       @changes changes Latest Changes

       @section motivation Motivation
        Intitially designed to avoid the usage of the FiniteVolumeManager in csp5.
        Output as nodal field variables allows to use them for post-processing (e.g. with ParaView or MatLab).

       @section usage Usage
        Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

       @code
       The visitor reads in a nodal variable for porosity, which has to be initialised before visitation
       It further requires existing nodal variables for volume and pore volume.
       The three variable names are provided as constructor arguments.

       The visitor calculates the node-centered volume and pore volume, assuming that all elements contribute equally to the control volumes of their respecetive nodes.

       The values for volume and pore volume are stored as nodal variables.

       @endcode

       @section dependencies Dependencies

       @section issues Known issues
       There may be new funcionality in CSMP++ that could make this visitor redundant.
       Only Visit(Node<dim>* n) is overwritten. However, application target and level are always set to NODE.
       How to limit to Region (former Group)?

       Does the Application Level have to be changed?
       How can it be restricted to a Region?
       Include the possibility to use porosity as an element variable?
       Change name to NodalCVFEM_PorosityChangeVisitor?

       @section testing Testing
       testing was done in the period before publication in 2014.

    */

} // csmp

#endif
