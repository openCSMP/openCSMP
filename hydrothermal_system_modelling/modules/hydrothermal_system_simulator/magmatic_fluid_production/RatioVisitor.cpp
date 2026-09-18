// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "RatioVisitor.h"

using namespace std;

namespace csmp
{

  template<size_t dim>
  RatioVisitor<dim>::RatioVisitor( Model<dim> &sg )
    : open_boundaries(false),
      dt(0.)
  {
    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);

    // mass

    mr_key = sg.Database().StorageKey("magmatic ratio");
    m_wr_key = sg.Database().StorageKey("magmatic water ratio");
    mt_key = sg.Database().StorageKey("fluid density");
    ml_key    = sg.Database().StorageKey("fluid mass liquid");
    mv_key    = sg.Database().StorageKey("fluid mass vapor");
    mml_key   = sg.Database().StorageKey("liquid mass mobility");
    mmv_key   = sg.Database().StorageKey("vapor mass mobility");
    m_m_key   = sg.Database().StorageKey("magmatic fluid mass");
    m_ml_key  = sg.Database().StorageKey("magmatic fluid mass liquid");
    m_mv_key  = sg.Database().StorageKey("magmatic fluid mass vapor");
    m_mml_key = sg.Database().StorageKey("magmatic liquid mass mobility");
    m_mmv_key = sg.Database().StorageKey("magmatic vapor mass mobility");
    m_mlp_key = sg.Database().StorageKey("previous magmatic fluid mass liquid");
    m_mvp_key = sg.Database().StorageKey("previous magmatic fluid mass vapor");
    src_rate_key = sg.Database().StorageKey("fluid source rate");

    // salt
    m_sr_key   = sg.Database().StorageKey("magmatic salt ratio");
    msp_key    = sg.Database().StorageKey("previous mass salt");
    xCl_key    = sg.Database().StorageKey("salt content liquid");
    xCv_key    = sg.Database().StorageKey("salt content vapor");
    xCh_key    = sg.Database().StorageKey("salt content halite");
    xml_key    = sg.Database().StorageKey("liquid salt mobility");
    xmv_key    = sg.Database().StorageKey("vapor salt mobility");
    m_ms_key   = sg.Database().StorageKey("magmatic mass salt");
    m_xCl_key  = sg.Database().StorageKey("magmatic salt content liquid");
    m_xCv_key  = sg.Database().StorageKey("magmatic salt content vapor");
    m_xCh_key  = sg.Database().StorageKey("magmatic salt content halite");
    m_xml_key  = sg.Database().StorageKey("magmatic liquid salt mobility");
    m_xmv_key  = sg.Database().StorageKey("magmatic vapor salt mobility");
    m_xClp_key = sg.Database().StorageKey("previous magmatic salt content liquid");
    m_xCvp_key = sg.Database().StorageKey("previous magmatic salt content vapor");
    src_wt_key = sg.Database().StorageKey("fluid source wt");

    // copper
    cff_key = sg.Database().StorageKey("copper fraction fluid");
    cCf_key   = sg.Database().StorageKey("copper content fluid");
    cCl_key  = sg.Database().StorageKey("copper content liquid");
    cCv_key  = sg.Database().StorageKey("copper content vapor");
    cml_key = sg.Database().StorageKey("liquid copper mobility");
    cmv_key = sg.Database().StorageKey("vapor copper mobility");
    cClp_key = sg.Database().StorageKey("previous copper content liquid");
    cCvp_key = sg.Database().StorageKey("previous copper content vapor");

    // boundary variables
    bfm_key    = sg.Database().StorageKey("boundary flow mass");
    bfs_key    = sg.Database().StorageKey("boundary flow salt");

    pv_key      = sg.Database().StorageKey("pore volume");

    if ( mr_key.type != SCALAR || mr_key.place != NODE )
      throw Exception(FATAL_ERROR, "RatioVisitor::(constructor)",
                      "magmatic ratio", " must be a nodal scalar property." );

    if ( mt_key.type != SCALAR || mt_key.place != NODE )
      throw Exception(FATAL_ERROR, "RatioVisitor::(constructor)",
                      "fluid density", " must be a nodal scalar property." );

    // TO DO:
    // add more checks

  }


  template<size_t dim>
  RatioVisitor<dim>::~RatioVisitor()
  {}

  /** visit function for Region */
  template<size_t dim>
  void RatioVisitor<dim>::Visit(Model<dim> *n)
  {
    // no calcution for the region
  }

  template<size_t dim>
  void RatioVisitor<dim>::Visit(Node<dim> *n)
  {

    if (n->Status( ml_key ) != DIRICH || open_boundaries)
      {

        ReadVariables( n );

        PerformMagmaticFluidAdvection();
        AddBoundaryTerms();
        AddSourceTerms();
        CalculateMagmaticRatio();
        CalculateMagmaticSaltRatio( );
        CalculateCopperMassFraction();
        SubtractBoundaryTerms();
        CalculateMagmaticTransportVariables();
        StoreVariables( n );

      }

  }


  template<size_t dim>
  void RatioVisitor<dim>::ReadVariables(Node<dim> *n)
  {

    // ratios
    n->Read(mr_key, mr );
    n->Read(m_wr_key, m_wr);
    n->Read(m_sr_key, m_sr);

    // mass advection variables
    n->Read(mt_key, mt );
    n->Read(ml_key, ml );
    n->Read(mv_key, mv );
    n->Read(mml_key, mml );
    n->Read(mmv_key, mmv );

    // magmatic mass advection variables
    n->Read(m_m_key, m_m );
    n->Read(m_ml_key, m_ml );
    n->Read(m_mv_key, m_mv );
    n->Read(m_mlp_key, m_mlp );
    n->Read(m_mvp_key, m_mvp );
    n->Read(m_mml_key, m_mml );
    n->Read(m_mmv_key, m_mmv );

    // salt advection variables
    n->Read(msp_key, msp );
    n->Read(xCl_key, xCl );
    n->Read(xCv_key, xCv );
    n->Read(xCh_key, xCh );
    n->Read(xml_key, xml );
    n->Read(xmv_key, xmv );

    // copper advection variables
    n->Read(cCf_key, cCf );
    n->Read(cCl_key, cCl );
    n->Read(cCv_key, cCv );
    n->Read(cClp_key, cClp );
    n->Read(cCvp_key, cCvp );
    n->Read(cml_key, cml );
    n->Read(cmv_key, cmv );

    // magmatic salt advection variables
    n->Read(m_ms_key, m_ms );
    n->Read(m_xCl_key, m_xCl );
    n->Read(m_xCv_key, m_xCv );
    n->Read(m_xCh_key, m_xCh );
    n->Read(m_xml_key, m_xml );
    n->Read(m_xmv_key, m_xmv );
    n->Read(m_xClp_key, m_xClp );
    n->Read(m_xCvp_key, m_xCvp );

    // boundary flow
    n->Read(bfm_key, bfm );
    n->Read(bfs_key, bfs );

    // constant source terms
    n->Read(src_rate_key, src_rate );
    n->Read(src_wt_key, src_wt );
    n->Read(pv_key, pv );

  }

  template<size_t dim>
  void RatioVisitor<dim>::PerformMagmaticFluidAdvection()
  {
    m_m()  += (m_ml() - m_mlp());
    m_m()  += (m_mv() - m_mvp());

    m_ms() += (m_xCl() - m_xClp());
    m_ms() += (m_xCv() - m_xCvp());

    cCf()  += (cCl() - cClp());
    cCf()  += (cCv() - cCvp());

  }

  template<size_t dim>
  void RatioVisitor<dim>::CalculateMagmaticSaltRatio( )
  {

    extern ErrorHandler   skm_err;

    if (essentiallyEqual(m_ms(), 0., numeric_limits<double>::epsilon()))
      {
        m_sr() = 0.0;
        m_ms() = 0.0;
      }
    else if (essentiallyEqual(m_ms(), msp(), numeric_limits<double>::epsilon()))
      {
        m_sr() = 1.0;
        m_ms() = msp();
      }
    else if (definitelyGreaterThan( m_ms(), msp(), numeric_limits<double>::epsilon() ))
      {
        m_ms() = msp();
        //      throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticSaltRatio",
        //                     "magmatic mass salt is greater than total mass salt." );
        //      int stop;
        //      cout << "Node at " << n->x() << ", " << n->y() << endl;
        //      cout << "m_ms: " << m_ms << ", msp: " << msp << endl;
        //      cin >> stop;
      }
    else if (definitelyLessThan( m_ms(), 0.0, numeric_limits<double>::epsilon() ))
      {
        m_ms() = 0.0;
        //      throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticSaltRatio",
        //                     "magmatic mass salt is less than zero." );
      }
    else if (definitelyGreaterThan( msp(), 0.0, numeric_limits<double>::epsilon() ))
      {
        m_sr() = m_ms() / msp();
      }
    else if (definitelyLessThan( msp(), 0.0, numeric_limits<double>::epsilon() ))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticSaltRatio",
                        "mass salt is less than zero." );
      }
    else
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticSaltRatio",
                        "magmatic salt ratio could not be calculated." );
      }

  }

  template<size_t dim>
  void RatioVisitor<dim>::CalculateCopperMassFraction( )
  {

    extern ErrorHandler   skm_err;

    if (essentiallyEqual(cCf(), 0., numeric_limits<double>::epsilon()))
      {
        cff() = 0.0;
      }
    else if (essentiallyEqual(mt(), 0., numeric_limits<double>::epsilon()))
      {
        cff() = 0.0;
      }
    else if (essentiallyEqual(cCf(), mt(), numeric_limits<double>::epsilon()))
      {
        cff() = 1.0;
      }
    else if (definitelyGreaterThan( cCf(), mt(), numeric_limits<double>::epsilon() ))
      {
        cff() = 1.0;
      }
    else if (definitelyLessThan( cCf(), 0.0, numeric_limits<double>::epsilon() ))
      {
        cff() = 0.0;
      }
    else if (definitelyGreaterThan( cCf(), 0.0, numeric_limits<double>::epsilon() ))
      {
        cff() = cCf() / mt();
      }
    else
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateCopperContent",
                        "copper fraction fluid could not be calculated." );
      }

  } // end CalculateCopperMassFraction

  template<size_t dim>
  void RatioVisitor<dim>::CalculateMagmaticRatio()
  {

    extern ErrorHandler   skm_err;

    if (essentiallyEqual(m_m(), 0., numeric_limits<double>::epsilon()))
      {
        mr() = 0.0;
        m_wr() = 0.0;
        m_m() = 0.0;
      }
    else if (essentiallyEqual(m_m(), mt(), numeric_limits<double>::epsilon()))
      {
        mr() = 1.0;
        m_wr() = 1.0;
        m_m() = mt();
      }
    else if (definitelyGreaterThan( m_m(), mt(), numeric_limits<double>::epsilon() ))
      {
        m_m() = mt();
      }
    else if (definitelyLessThan( m_m(), 0.0, numeric_limits<double>::epsilon() ))
      {
        m_m() = 0.0;
      }
    else if (definitelyGreaterThan( mt(), 0.0, numeric_limits<double>::epsilon() ))
      {
        // ratio of bulk fluid
        mr() = m_m() / mt();

        // ratio of H2O component
        if (definitelyGreaterThan(mt(), msp(), numeric_limits<double>::epsilon()))
          m_wr() = (m_m() - m_ms()) / (mt() - msp());
        else
          m_wr() = 0.0;
      }
    else if (essentiallyEqual(mt(), 0.0, numeric_limits<double>::epsilon()))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticRatio",
                        "fluid density is zero." );
      }
    else if (definitelyLessThan( mt(), 0.0, numeric_limits<double>::epsilon() ))
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticRatio",
                        "fluid density is less than zero." );
      }
    else
      {
        throw Exception(FATAL_ERROR, "RatioVisitor::CalculateMagmaticRatio",
                        "magmatic ratio could not be calculated." );
      }

  }

  template<size_t dim>
  void RatioVisitor<dim>::CalculateMagmaticTransportVariables()
  {

    // salt
    m_xCl() = xCl() * m_sr();
    m_xCv() = xCv() * m_sr();
    m_xCh() = xCh() * m_sr();
    m_xml() = xml() * m_sr();
    m_xmv() = xmv() * m_sr();

    // mass
    //      m_ml()  = ml()  * mr();
    //      m_mv()  = mv()  * mr();
    //      m_mml() = mml() * mr();
    //      m_mmv() = mmv() * mr();
    m_ml() = ( ml() - xCl() ) * m_wr() + m_xCl();
    m_mv() = ( mv() - xCv() ) * m_wr() + m_xCv();
    m_mml() = mml() * m_ml() / ml();
    m_mmv() = mmv() * m_mv() / mv();

    // copper - no preferential partitioning
    cCl()  = ml()  * cff();
    cCv()  = mv()  * cff();
    cml() = mml() * cff();
    cmv() = mmv() * cff();

    // in include some consistency checks?

  }

  template<size_t dim>
  void RatioVisitor<dim>::StoreVariables(Node<dim> *n)
  {

    // ratios
    n->Store(mr_key, mr );
    n->Store(m_wr_key, m_wr);
    n->Store(m_sr_key, m_sr);
    n->Store(cff_key,   cff );

    // magmatic mass advection variables
    n->Store(m_m_key,   m_m );
    n->Store(m_ml_key,  m_ml );
    n->Store(m_mv_key,  m_mv );
    n->Store(m_mlp_key, m_mlp );
    n->Store(m_mvp_key, m_mvp );
    n->Store(m_mml_key, m_mml );
    n->Store(m_mmv_key, m_mmv );

    // magmatic salt advection variables
    n->Store(m_ms_key,   m_ms );
    n->Store(m_xCl_key,  m_xCl );
    n->Store(m_xCv_key,  m_xCv );
    n->Store(m_xCh_key,  m_xCh );
    n->Store(m_xml_key,  m_xml );
    n->Store(m_xmv_key,  m_xmv );
    n->Store(m_xClp_key, m_xClp );
    n->Store(m_xCvp_key, m_xCvp );

    // copper advection variables
    n->Store(cCf_key,   cCf );
    n->Store(cCl_key,  cCl );
    n->Store(cCv_key,  cCv );
    n->Store(cClp_key, cClp );
    n->Store(cCvp_key, cCvp );
    n->Store(cml_key, cml );
    n->Store(cmv_key, cmv );

  }

  template<size_t dim>
  void RatioVisitor<dim>::AddBoundaryTerms()
  {
    mt() += std::max(0.0, bfm());
    msp() += std::max(0.0, bfs());
  }

  template<size_t dim>
  void RatioVisitor<dim>::AddSourceTerms()
  {
    new_mass = src_rate() * dt / pv();
    m_m() += new_mass;
    m_ms() += new_mass * src_wt() * 0.01;
    // copper source term?
  }

  template<size_t dim>
  void RatioVisitor<dim>::SubtractBoundaryTerms()
  {
    m_m()  -= std::max(0.0, bfm() * mr());
    m_ms() -= std::max(0.0, bfs() * m_sr());
    cCf()  -= std::max(0.0, bfm() * cff());
  }

  template<size_t dim>
  void RatioVisitor<dim>::SetOpenTopTo( bool open)
  {
    open_boundaries = open;
  }

  template<size_t dim>
  void RatioVisitor<dim>::SetTimeIncrement( double time_increment )
  {
    dt = time_increment;
  }

  template class RatioVisitor<1U>;
  template class RatioVisitor<2U>;
  template class RatioVisitor<3U>;

} // csmp
