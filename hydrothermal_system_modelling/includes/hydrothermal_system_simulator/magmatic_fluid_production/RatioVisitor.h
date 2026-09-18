// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RatioVisitor_h
#define RatioVisitor_h

#include "Visitor.h"
#include "Model.h"
#include "compareFloats.h"
#include "ErrorHandler.h"

namespace csmp
{

  template<size_t dim>
  class RatioVisitor : public Visitor<dim>
  {
    public:
      RatioVisitor( Model<dim> &sg);
      ~RatioVisitor();

      virtual void Visit(Node<dim> *n);
      virtual void Visit(Model<dim> *n);
      void SetOpenTopTo(bool open);
      void SetTimeIncrement( double time_increment );

    private:

      void ReadVariables( Node<dim> *n );
      void PerformMagmaticFluidAdvection();
      void CalculateMagmaticRatio();
      void CalculateMagmaticSaltRatio( );
      void CalculateCopperMassFraction( );
      void CalculateMagmaticTransportVariables();
      void StoreVariables( Node<dim> *n );
      void AddSourceTerms( );
      void AddBoundaryTerms( );
      void SubtractBoundaryTerms( );

      bool open_boundaries;
      double dt, new_mass;

      // mass
      csmp::Index         mr_key, m_wr_key, mt_key, ml_key, mv_key, mml_key, mmv_key,
           m_m_key, m_ml_key, m_mv_key, m_mml_key, m_mmv_key,
           m_mlp_key, m_mvp_key, src_rate_key;
      ScalarVariable  mt, ml, mv, mml, mmv, m_m, m_ml, m_mv, m_mml,
                      mr, m_wr, m_mmv, m_mlp, m_mvp, src_rate;

      // salt
      csmp::Index          msp_key, xCl_key, xCv_key, xCh_key, xml_key, xmv_key,
           m_ms_key, m_xCl_key, m_xCv_key, m_xCh_key, m_xml_key, m_xmv_key,
           m_sr_key, m_xClp_key, m_xCvp_key, src_wt_key;
      ScalarVariable  msp, xCl, xCv, xCh, xml, xmv,
                      m_ms, m_xCl, m_xCv, m_xCh, m_xml, m_xmv,
                      m_sr, m_xClp, m_xCvp, src_wt;
      //copper
      csmp::Index          cCf_key, cCl_key, cCv_key, cClp_key, cCvp_key,
           cff_key, cml_key, cmv_key;
      ScalarVariable  cCf, cCl, cCv, cClp, cCvp,
                      cff, cml, cmv;

      // open top
      csmp::Index         bfm_key, bfs_key, pv_key;
      ScalarVariable bfm, bfs, pv;

  };

} // csmp

#endif
