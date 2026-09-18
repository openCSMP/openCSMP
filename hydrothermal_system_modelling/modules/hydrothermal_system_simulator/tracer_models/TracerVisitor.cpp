// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TracerVisitor.h"

using namespace std;

namespace csmp
{

  template <uint32_t dim>
  TracerVisitor<dim>::TracerVisitor( Model<dim> &model)
    : Visitor<dim>( MODEL, NODE ),
      model_ref_ ( model )

  {
    c_tr_fl_key = model.Database().StorageKey("tracer content fluid");
    c_tr_l_key = model.Database().StorageKey("tracer content liquid");
    c_tr_v_key = model.Database().StorageKey("tracer content vapor");

    massfrac_tr_l_key  = model.Database().StorageKey("tracer fraction liquid");
    massfrac_tr_v_key = model.Database().StorageKey("tracer fraction vapor");
    massfrac_tr_fl_key = model.Database().StorageKey("tracer fraction fluid");

    mob_tr_l_key  = model.Database().StorageKey("liquid tracer mobility");
    mob_tr_v_key = model.Database().StorageKey("vapor tracer mobility");

    tr_mass_out_total_key = model.Database().StorageKey("tracer mass out total");

    bfm_key = model.Database().StorageKey("boundary flow mass");
    pore_volume_key = model.Database().StorageKey("pore volume");

    mml_key = model.Database().StorageKey("liquid mass mobility"); //
    mmv_key = model.Database().StorageKey("vapor mass mobility"); //

    ml_key = model.Database().StorageKey("fluid mass liquid"); //
    mv_key = model.Database().StorageKey("fluid mass vapor"); //
    mt_key = model.Database().StorageKey("fluid density"); //

    temperature_key = model.Database().StorageKey("temperature");
    fluid_pressure_key = model.Database().StorageKey("fluid pressure");
    lithostatic_pressure_key = model.Database().StorageKey("lithostatic pressure");

    coord_x_key = model.Database().StorageKey("coordinate x");
    coord_y_key = model.Database().StorageKey("coordinate y");

    fluid_enthalpy_key = model.Database().StorageKey("fluid enthalpy");
  }

  template<uint32_t dim>
  TracerVisitor<dim>::~TracerVisitor()
  {
  }


  template <uint32_t dim>
  void TracerVisitor<dim>::Visit(Node<dim> *n)
  {
  }


  template<uint32_t dim>
  void TracerVisitor<dim>::CalculateMassFractions(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(c_tr_l_key, c_tr_l);
        (*nit)->Read(c_tr_v_key, c_tr_v);

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        if (ml() != 0.)
          {
            massfrac_tr_l()  =   c_tr_l() / ml();
          }
        else
          {
            massfrac_tr_l()  =   0.;
          }

        if (mv() != 0.)
          {
            massfrac_tr_v()  =   c_tr_v() / mv();
          }
        else
          {
            massfrac_tr_v()  =   0.;
          }

        massfrac_tr_fl() = (c_tr_l() + c_tr_v()) / (ml() + mv());
        massfrac_tr_fl_rescaled() = (c_tr_l() + c_tr_v()) / (ml() + mv());

        (*nit)->Store(massfrac_tr_l_key, massfrac_tr_l);
        (*nit)->Store(massfrac_tr_v_key, massfrac_tr_v);
        (*nit)->Store(massfrac_tr_fl_key, massfrac_tr_fl);
      }
  }


  template<uint32_t dim>
  void TracerVisitor<dim>::CalculateMobilities(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(mml_key, mml );
        (*nit)->Read(mmv_key, mmv );

        (*nit)->Read(massfrac_tr_l_key, massfrac_tr_l);
        (*nit)->Read(massfrac_tr_v_key, massfrac_tr_v);

        mob_tr_l() = mml() * massfrac_tr_l();
        mob_tr_v() = mmv() * massfrac_tr_v();

        (*nit)->Store(mob_tr_l_key, mob_tr_l);
        (*nit)->Store(mob_tr_v_key, mob_tr_v);
      }
  }

  template<uint32_t dim>
  void TracerVisitor<dim>::SetTracerContent(Model<dim> *model, double timestep)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(temperature_key, temperature );
        (*nit)->Read(fluid_pressure_key, fluid_pressure );

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(c_tr_l_key, c_tr_l);
        (*nit)->Read(c_tr_v_key, c_tr_v);

        (*nit)->Read(fluid_enthalpy_key, fluid_enthalpy);


        //        double rate=1.E-6/31540000.; // 1 ppm tracer per year

        if (temperature() > 374. && fluid_enthalpy() > 20860000.) // example for super critical regions
          {
            //            if (ml()!=0.) c_tr_l() += rate*timestep*ml(); else c_tr_l() = 0.;
            //            if (mv()!=0.) c_tr_v() += rate*timestep*mv(); else c_tr_v() = 0.;
            if (ml() != 0.)
              c_tr_l() = 1.E-6 * ml();
            else
              c_tr_l() = 0.;

            if (mv() != 0.)
              c_tr_v() = 1.E-6 * mv();
            else
              c_tr_v() = 0.;
          }

        (*nit)->Store(c_tr_l_key, c_tr_l);
        (*nit)->Store(c_tr_v_key, c_tr_v);
      }
  }

  template<uint32_t dim>
  void TracerVisitor<dim>::InitializeTracerVisitor(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        c_tr_l() = 0.;
        c_tr_v() = 0.;

        (*nit)->Store(c_tr_l_key, c_tr_l);
        (*nit)->Store(c_tr_v_key, c_tr_v);
      }
  }

  template<uint32_t dim>
  void TracerVisitor<dim>::DistributeTracerBetweenPhases(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(c_tr_l_key, c_tr_l);
        (*nit)->Read(c_tr_v_key, c_tr_v);

        double massfractionwaterliquid = ml() / (ml() + mv());

        c_tr_fl() = c_tr_l() + c_tr_v();

        //Assuming equal partitioning between vapor and liquid
        c_tr_l() = massfractionwaterliquid * c_tr_fl();
        c_tr_v() = (1. - massfractionwaterliquid) * c_tr_fl();

        (*nit)->Store(c_tr_l_key, c_tr_l);
        (*nit)->Store(c_tr_v_key, c_tr_v);
        (*nit)->Store(c_tr_fl_key, c_tr_fl);
      }
  }

  template<uint32_t dim>
  void TracerVisitor<dim>::ApplyBoundaryFlowCorrections(Model<dim> *model)
  {
    model_ref_.Read( tr_mass_out_total_key, tr_mass_out_total );
    double tottrout(tr_mass_out_total());

    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(bfm_key, bfm);

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(c_tr_l_key, c_tr_l);
        (*nit)->Read(c_tr_v_key, c_tr_v);

        (*nit)->Read(mml_key, mml);
        (*nit)->Read(mmv_key, mmv);

        double scalingfactor(0.);
        double tr_frac_fluid_inflow(0.); //must be adapted if recharge water has tracer content

        double tr_mass_before = pore_volume() * (c_tr_l() + c_tr_v());

        if (bfm() > 0.) //outflow
          {
            scalingfactor = 1. / (1. + ((mml() / (mmv() + mml())) * bfm() / (mv() + ml())));

            //NK may 2024, : the relative mass proportions of vapor and liquid contributing to the total boundary flow should be approximately proportional
            //to the relative proportions of mobilities. Since the previous mobilities are not accessible, we take the current mobilities
            //as an approximation. Therefore the contribution of the liquid to the total boundary flow is (mmv/(mmv+mml))*bfm,
            //analogous for vapor

            c_tr_l() *= scalingfactor;

            scalingfactor = 1. / (1. + ((mmv() / (mmv() + mml())) * bfm() / (mv() + ml())));
            c_tr_v() *= scalingfactor;
          }
        else //inflow
          {
            c_tr_l() -= bfm() * tr_frac_fluid_inflow;
            c_tr_v() -= bfm() * tr_frac_fluid_inflow;
          }

        c_tr_fl() = c_tr_l() + c_tr_v();

        tottrout += tr_mass_before - pore_volume() * c_tr_fl();

        (*nit)->Store(c_tr_l_key, c_tr_l);
        (*nit)->Store(c_tr_v_key, c_tr_v);
        (*nit)->Store(c_tr_fl_key, c_tr_fl);
      }

    tr_mass_out_total() += tottrout;
    model_ref_.Store( tr_mass_out_total_key, tr_mass_out_total );
  }


  template<uint32_t dim>
  void TracerVisitor<dim>::Visit( Region<dim> *r )
  {}

  template<uint32_t dim>
  void TracerVisitor<dim>::Visit( Model<dim> *m )
  {}

  template<uint32_t dim>
  void TracerVisitor<dim>::Visit( Element<dim> *e )
  {}




  template class TracerVisitor<1>;
  template class TracerVisitor<2>;
  template class TracerVisitor<3>;

}// csmp
