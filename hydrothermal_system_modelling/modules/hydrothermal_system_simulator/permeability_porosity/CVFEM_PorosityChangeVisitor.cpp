// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_PorosityChangeVisitor.h"

using namespace std;

namespace csmp
{

  /** custom constructor */
  template<size_t dim>
  CVFEM_PorosityChangeVisitor<dim>::CVFEM_PorosityChangeVisitor( Model<dim> &model )
  {

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    pore_volume_change_factor_key = model.Database().StorageKey("pore volume change factor"),//Scaling variable
    pore_vol_key = model.Database().StorageKey("pore volume");                //used to calculate nQ
    phi_key = model.Database().StorageKey("nodal porosity");             //used to calculate ncp and CT

    rr_key = model.Database().StorageKey("nodal density rock");         //used to calculate ncp
    ncp_key = model.Database().StorageKey("nodal heat capacity");        //CALCULATED LATER IN EQUILIBRATOR!
    cpr_key = model.Database().StorageKey("nodal heat capacity rock");   //used to calculate ncp
    cpf_key = model.Database().StorageKey("fluid heat capacity");        //used to calculate ncp
    nQ_key = model.Database().StorageKey("nodal fluid volume source");  //CALCULATED LATER IN EQUILIBRATOR!
    CT_key = model.Database().StorageKey("nodal total compressibility");//CALCULATED LATER IN EQUILIBRATOR!
    rho_bulk_key = model.Database().StorageKey("bulk fluid density");         //used to calculate nQ
    betaf_key = model.Database().StorageKey("compressibility");            //used to calculate CT
    betar_key = model.Database().StorageKey("nodal compressibility rock"); //used to calculate CT

    mt_key = model.Database().StorageKey("fluid density");              //"LHS"
    msp_key = model.Database().StorageKey("previous mass salt");         //"LHS"
    Htp_key = model.Database().StorageKey("previous total enthalpy");    //"LHS", not used at the moment
    sh_key  = model.Database().StorageKey("saturation halite");          //used to calculate CT

    // transport variable LHS, all ok
    ml_key = model.Database().StorageKey("fluid mass liquid");            //LHS
    mv_key = model.Database().StorageKey("fluid mass vapor");             //LHS
    mf_key = model.Database().StorageKey("fluid mass");                   //"LHS" = mt() - mh() in equilibrator;
    hCl_key = model.Database().StorageKey("enthalpy content liquid");      //LHS
    hCv_key = model.Database().StorageKey("enthalpy content vapor");       //LHS
    xCl_key = model.Database().StorageKey("salt content liquid");          //LHS
    xCv_key = model.Database().StorageKey("salt content vapor");           //LHS

    rl_transport_key = model.Database().StorageKey("density liquid transport");   //LHS ?? SHOULD IT BE SCALED?, no previous_V exist
    rv_transport_key = model.Database().StorageKey("density vapor transport");    //LHS ?? SHOULD IT BE SCALED?, no previous_V exist

    // variables for transport calculation
    mml_key =  model.Database().StorageKey("liquid mass mobility");
    mmv_key =  model.Database().StorageKey("vapor mass mobility");
    mmld_key = model.Database().StorageKey("liquid mass mobility density"); // = mml() * rl_transport() in equilibrator, rl_transport()=rl()*volfactRHS
    mmvd_key = model.Database().StorageKey("vapor mass mobility density");  // = mmv() * rv_transport() in equilibrator, rv_transport()=rv()*volfactRHS
    eml_key =  model.Database().StorageKey("liquid enthalpy mobility");     //
    emv_key =  model.Database().StorageKey("vapor enthalpy mobility");      //
    emld_key = model.Database().StorageKey("liquid enthalpy mobility density");// = eml() * rl_transport(); only used and calculated in equilibrator
    emvd_key = model.Database().StorageKey("vapor enthalpy mobility density"); // = emv() * rv_transport(); only used and calculated in equilibrator
    xml_key =  model.Database().StorageKey("liquid salt mobility");         //RHS
    xmv_key =  model.Database().StorageKey("vapor salt mobility");          //RHS


    if (pore_vol_key.type != SCALAR || pore_vol_key.place != NODE )
      throw csmp::Exception( FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                             "pore volume", " must be a nodal scalar property." );

    if ( phi_key.type != SCALAR || phi_key.place != NODE )
      throw csmp::Exception( FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                             "nodal porosity", " must be a scalar property." );

    if (rr_key.type != SCALAR || rr_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "nodal density rock", " must be a scalar property.");

    if (ncp_key.type != SCALAR || ncp_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "nodal heat capacity", " must be a scalar property.");

    if (cpr_key.type != SCALAR || cpr_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "nodal heat capacity rock", " must be a scalar property.");

    if (mt_key.type != SCALAR || mt_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "fluid density", " must be a scalar property.");

    if (msp_key.type != SCALAR || msp_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "previous mass salt", " must be a scalar property.");

    if (Htp_key.type != SCALAR || Htp_key.place != NODE)
      throw csmp::Exception(FATAL_ERROR, "CVFEM_PorosityChangeVisitor::(constructor)",
                            "previous total enthalpy", " must be a scalar property.");

  }

  /** deconstructor */
  template<size_t dim>
  CVFEM_PorosityChangeVisitor<dim>::~CVFEM_PorosityChangeVisitor()
  {}

  /** visit function for Node */
  template<size_t dim>
  void CVFEM_PorosityChangeVisitor<dim>::Visit(Node<dim> *n)
  {
    n->Read(pore_volume_change_factor_key, pore_volume_change_factor);
    n->Read(pore_vol_key, pore_vol);
    n->Read(phi_key, phi );
    n->Read(rr_key, rr);
    n->Read(cpr_key, cpr);
    n->Read(cpf_key, cpf);
    n->Read(rho_bulk_key, rho_bulk);
    n->Read(betaf_key, betaf);
    n->Read(betar_key, betar);
    n->Read(mt_key, mt);
    n->Read(msp_key, msp);
    n->Read(Htp_key, Htp);
    n->Read(sh_key, sh);

    // transport variable LHS
    n->Read(ml_key, ml);
    n->Read(mv_key, mv);
    n->Read(mf_key, mf);
    n->Read(hCl_key, hCl);
    n->Read(hCv_key, hCv);
    n->Read(xCl_key, xCl);
    n->Read(xCv_key, xCv);

    // transport magmatic and copper variables LHS
    if (with_magma)
      {
        n->Read(m_ml_key, m_ml);
        n->Read(m_mv_key, m_mv);
        n->Read(m_m_key, m_m);
        n->Read(m_xCl_key, m_xCl);
        n->Read(m_xCv_key, m_xCv);
        n->Read(m_ms_key, m_ms);
      }

    if (with_gold)
      {
        n->Read(so2Cl_key, so2Cl);
        n->Read(so2Cv_key, so2Cv);
        n->Read(so2Cf_key, so2Cf);
        n->Read(auhs0Cl_key, auhs0Cl);
        n->Read(auhs0Cv_key, auhs0Cv);
        n->Read(auhs0Cf_key, auhs0Cf);

        n->Read(so2ml_key, so2ml);
        n->Read(so2mv_key, so2mv);
        n->Read(auhs0ml_key, auhs0ml);
        n->Read(auhs0mv_key, auhs0mv);

        so2Cl() /= pore_volume_change_factor();
        so2Cv() /= pore_volume_change_factor();
        so2Cf() /= pore_volume_change_factor();
        auhs0Cl() /= pore_volume_change_factor();
        auhs0Cv() /= pore_volume_change_factor();
        auhs0Cf() /= pore_volume_change_factor();

        so2ml() /= pore_volume_change_factor();
        so2mv() /= pore_volume_change_factor();
        auhs0ml() /= pore_volume_change_factor();
        auhs0mv() /= pore_volume_change_factor();

        n->Store(so2Cl_key, so2Cl);
        n->Store(so2Cv_key, so2Cv);
        n->Store(so2Cf_key, so2Cf);
        n->Store(auhs0Cl_key, auhs0Cl);
        n->Store(auhs0Cv_key, auhs0Cv);
        n->Store(auhs0Cf_key, auhs0Cf);

        n->Store(so2ml_key, so2ml);
        n->Store(so2mv_key, so2mv);
        n->Store(auhs0ml_key, auhs0ml);
        n->Store(auhs0mv_key, auhs0mv);
      }

    if (with_gold || with_quartz)
      {
        n->Read(sio2Cl_key, sio2Cl);
        n->Read(sio2Cv_key, sio2Cv);
        n->Read(sio2Cf_key, sio2Cf);

        n->Read(sio2ml_key, sio2ml);
        n->Read(sio2mv_key, sio2mv);

        sio2Cl() /= pore_volume_change_factor();
        sio2Cv() /= pore_volume_change_factor();
        sio2Cf() /= pore_volume_change_factor();

        sio2ml() /= pore_volume_change_factor();
        sio2mv() /= pore_volume_change_factor();

        n->Store(sio2Cl_key, sio2Cl);
        n->Store(sio2Cv_key, sio2Cv);
        n->Store(sio2Cf_key, sio2Cf);

        n->Store(sio2ml_key, sio2ml);
        n->Store(sio2mv_key, sio2mv);
      }

    if (with_lithium)
      {
        n->Read(liCl_key, liCl);
        n->Read(liCv_key, liCv);
        n->Read(liCf_key, liCf);

        n->Read(liml_key, liml);
        n->Read(limv_key, limv);

        liCl() /= pore_volume_change_factor();
        liCv() /= pore_volume_change_factor();
        liCf() /= pore_volume_change_factor();

        liml() /= pore_volume_change_factor();
        limv() /= pore_volume_change_factor();

        n->Store(liCl_key, liCl);
        n->Store(liCv_key, liCv);
        n->Store(liCf_key, liCf);

        n->Store(liml_key, liml);
        n->Store(limv_key, limv);
      }

    if (with_tracer)
      {
        n->Read(tracerCl_key, tracerCl);
        n->Read(tracerCv_key, tracerCv);
        n->Read(tracerCf_key, tracerCf);

        n->Read(tracerml_key, tracerml);
        n->Read(tracermv_key, tracermv);

        tracerCl() /= pore_volume_change_factor();
        tracerCv() /= pore_volume_change_factor();
        tracerCf() /= pore_volume_change_factor();

        tracerml() /= pore_volume_change_factor();
        tracermv() /= pore_volume_change_factor();

        n->Store(tracerCl_key, tracerCl);
        n->Store(tracerCv_key, tracerCv);
        n->Store(tracerCf_key, tracerCf);

        n->Store(tracerml_key, tracerml);
        n->Store(tracermv_key, tracermv);
      }

    // density liquid/vapor transport (was not scaled previously)
    n->Read(rl_transport_key, rl_transport);
    n->Read(rv_transport_key, rv_transport);

    // transport variable
    n->Read(mml_key, mml);
    n->Read(mmv_key, mmv);
    n->Read(mmld_key, mmld);
    n->Read(mmvd_key, mmvd);
    n->Read(eml_key, eml);
    n->Read(emv_key, emv);
    n->Read(emld_key, emld);
    n->Read(emvd_key, emvd);
    n->Read(xml_key, xml);
    n->Read(xmv_key, xmv);

    // transport magmatic and copper variables
    if (with_magma)
    {
        n->Read(m_mml_key, m_mml);
        n->Read(m_mmv_key, m_mmv);
        n->Read(m_xml_key, m_xml);
        n->Read(m_xmv_key, m_xmv);
    }


    // Update porosity related variables
    ncp() = cpr() * rr() * (1. - phi());
    ncp() += cpf() * rho_bulk() * phi(); //BB add
    nQ() = (mt() - rho_bulk()) * pore_vol() * (1. - sh());

    // total compressibility, this is not used in poroelasticity
    CT()  = betaf() * phi() * (1. - sh());
    CT() += betar() * (1. - phi() + phi() * sh());
    //    CT() *= mt();
    CT() *= rho_bulk();
    //-----------------------------------
    //-----------------------------------

    // SCALING
    // main conservation variables
    //mt() /= pore_volume_change_factor();//Leave commented as this is already done before
    msp() /= pore_volume_change_factor();
    //Htp() -= d_phi()*rr()*cpr()*T(); //Use this if the change in porosity actually correspond to rock mass change (dissolution or precipitation)

    // transport variable LHS
    ml() /= pore_volume_change_factor();
    mv() /= pore_volume_change_factor();
    mf() /= pore_volume_change_factor();
    hCl() /= pore_volume_change_factor();
    hCv() /= pore_volume_change_factor();
    xCl() /= pore_volume_change_factor();
    xCv() /= pore_volume_change_factor();

    // transport magmatic and copper variables LHS
    if (with_magma)
      {
        m_ml() /= pore_volume_change_factor();
        m_mv() /= pore_volume_change_factor();
        m_m() /= pore_volume_change_factor();
        m_xCl() /= pore_volume_change_factor();
        m_xCv() /= pore_volume_change_factor();
        m_ms() /= pore_volume_change_factor();
      }

    // density liquid/vapor transport (was not scaled previously)
    //rl_transport/= pore_volume_change_factor();
    //rv_transport/= pore_volume_change_factor();

    // transport variable RHS
    mml() /= pore_volume_change_factor();
    mmv() /= pore_volume_change_factor();
    mmld() /= (pore_volume_change_factor() * pore_volume_change_factor());
    mmvd() /= (pore_volume_change_factor() * pore_volume_change_factor());
    eml() /= pore_volume_change_factor();
    emv() /= pore_volume_change_factor();
    emld() /= (pore_volume_change_factor() * pore_volume_change_factor());
    emvd() /= (pore_volume_change_factor() * pore_volume_change_factor());
    xml() /= pore_volume_change_factor();
    xmv() /= pore_volume_change_factor();

    // transport magmatic and copper variables RHS
    if (with_magma)
    {
        m_mml() /= pore_volume_change_factor();
        m_mmv() /= pore_volume_change_factor();
        m_xml() /= pore_volume_change_factor();
        m_xmv() /= pore_volume_change_factor();
    }
    //----------------------------------

    // Store
    // porosity and related variables
    n->Store(ncp_key, ncp);
    n->Store(nQ_key, nQ);
    n->Store(CT_key, CT);

    // main conservation variables
    n->Store(mt_key, mt);
    n->Store(msp_key, msp);
    n->Store(Htp_key, Htp);

    // transport variable LHS
    n->Store(ml_key, ml);
    n->Store(mv_key, mv);
    n->Store(mf_key, mf);
    n->Store(hCl_key, hCl);
    n->Store(hCv_key, hCv);
    n->Store(xCl_key, xCl);
    n->Store(xCv_key, xCv);

    // transport magmatic and copper variables LHS
    if (with_magma)
      {
        n->Store(m_ml_key, m_ml);
        n->Store(m_mv_key, m_mv);
        n->Store(m_m_key, m_m);
        n->Store(m_xCl_key, m_xCl);
        n->Store(m_xCv_key, m_xCv);
        n->Store(m_ms_key, m_ms);

        n->Store(m_mml_key, m_mml);
        n->Store(m_mmv_key, m_mmv);
        n->Store(m_xml_key, m_xml);
        n->Store(m_xmv_key, m_xmv);
      }

    // density liquid/vapor transport (was not scaled previously)
    n->Store(rl_transport_key, rl_transport);
    n->Store(rv_transport_key, rv_transport);

    // transport variable
    n->Store(mml_key, mml);
    n->Store(mmv_key, mmv);
    n->Store(mmld_key, mmld);
    n->Store(mmvd_key, mmvd);
    n->Store(eml_key, eml);
    n->Store(emv_key, emv);
    n->Store(emld_key, emld);
    n->Store(emvd_key, emvd);
    n->Store(xml_key, xml);
    n->Store(xmv_key, xmv);

  }

  /** visit function for Model */
  template<size_t dim>
  void CVFEM_PorosityChangeVisitor<dim>::Visit(Model<dim> *m)
  {
    // no calcution for the region
  }

  /** visit function for Region */
  template<size_t dim>
  void CVFEM_PorosityChangeVisitor<dim>::Visit(Region<dim> *r)
  {
    // no calcution for the region
  }

  template<size_t dim>
  void CVFEM_PorosityChangeVisitor<dim>::ReadInputVariables(
      Model<dim> &model,
      bool ext_with_gold,
      bool ext_with_quartz,
      bool ext_with_tracer,
      bool ext_with_lithium,
      bool ext_with_magma)
  {
    with_gold   = ext_with_gold;
    with_quartz = ext_with_quartz;
    with_lithium= ext_with_lithium;
    with_tracer = ext_with_tracer;
    with_magma  = ext_with_magma;

    //! create csmp indices for optional variables
    if (with_tracer)
    {
        tracerCl_key = model.Database().StorageKey("tracer content liquid");       //LHS,
        tracerCv_key = model.Database().StorageKey("tracer content vapor");        //LHS,
        tracerCf_key = model.Database().StorageKey("tracer content fluid");       //"LHS" = cCl+cCv, probably useless to scale

        tracerml_key = model.Database().StorageKey("liquid tracer mobility");
        tracermv_key = model.Database().StorageKey("vapor tracer mobility");
    }

    if (with_magma)
    {
        m_ml_key  = model.Database().StorageKey("magmatic fluid mass liquid");      //LHS
        m_mv_key  = model.Database().StorageKey("magmatic fluid mass vapor");       //LHS
        m_m_key   = model.Database().StorageKey("magmatic fluid mass");             //"LHS" = m_ml+m_mv   in ratio visitor, probably useless to scale
        m_xCl_key = model.Database().StorageKey("magmatic salt content liquid");    //LHS
        m_xCv_key = model.Database().StorageKey("magmatic salt content vapor");     //LHS
        m_ms_key  = model.Database().StorageKey("magmatic mass salt");              //"LHS" = m_xCl+m_xCv in ratio visitor, probably useless to scale

        m_mml_key = model.Database().StorageKey("magmatic liquid mass mobility");   //RHS
        m_mmv_key = model.Database().StorageKey("magmatic vapor mass mobility");    //RHS
        m_xml_key = model.Database().StorageKey("magmatic liquid salt mobility");   //RHS
        m_xmv_key = model.Database().StorageKey("magmatic vapor salt mobility");    //RHS
    }

    if (with_gold)
    {
        so2Cl_key   = model.Database().StorageKey("so2 content liquid");            //LHS,
        so2Cv_key   = model.Database().StorageKey("so2 content vapor");             //LHS,
        so2Cf_key   = model.Database().StorageKey("so2 content fluid");             //"LHS" = cCl+cCv, probably useless to scale
        auhs0Cl_key = model.Database().StorageKey("auhs0 content liquid");          //LHS,
        auhs0Cv_key = model.Database().StorageKey("auhs0 content vapor");           //LHS,
        auhs0Cf_key = model.Database().StorageKey("auhs0 content fluid");           //"LHS" = cCl+cCv, probably useless to scale

        so2ml_key   = model.Database().StorageKey("liquid so2 mobility");
        so2mv_key   = model.Database().StorageKey("vapor so2 mobility");
        auhs0ml_key = model.Database().StorageKey("liquid auhs0 mobility");
        auhs0mv_key = model.Database().StorageKey("vapor auhs0 mobility");
    }

    if (with_quartz)
    {
        sio2Cl_key = model.Database().StorageKey("sio2 content liquid");            //LHS,
        sio2Cv_key = model.Database().StorageKey("sio2 content vapor");             //LHS,
        sio2Cf_key = model.Database().StorageKey("sio2 content fluid");             //"LHS" = cCl+cCv, probably useless to scale

        sio2ml_key = model.Database().StorageKey("liquid sio2 mobility");
        sio2mv_key = model.Database().StorageKey("vapor sio2 mobility");
    }

    if (with_lithium)
    {
        liCl_key = model.Database().StorageKey("li content liquid");               //LHS,
        liCv_key = model.Database().StorageKey("li content vapor");                //LHS,
        liCf_key = model.Database().StorageKey("li content fluid");               //"LHS" = cCl+cCv, probably useless to scale

        liml_key = model.Database().StorageKey("liquid li mobility");
        limv_key = model.Database().StorageKey("vapor li mobility");
    }

  }

  template class CVFEM_PorosityChangeVisitor<1U>;
  template class CVFEM_PorosityChangeVisitor<2U>;
  template class CVFEM_PorosityChangeVisitor<3U>;

} // csmp
