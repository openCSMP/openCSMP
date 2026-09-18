// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "LithiumVisitor.h"


using namespace std;

namespace csmp
{

  template <uint32_t dim>
  LithiumVisitor<dim>::LithiumVisitor( Model<dim> &model)
    : Visitor<dim>( MODEL, NODE ),
      model_ref_ ( model )
  {

    li_mass_liquid_key    = model.Database().StorageKey("li mass liquid");
    li_mass_vapor_key     = model.Database().StorageKey("li mass vapor");
    li_mass_fluid_key     =  model.Database().StorageKey("li mass fluid");
    li_mass_melt_key      = model.Database().StorageKey("li mass melt");
    li_mass_crystals_key  = model.Database().StorageKey("li mass crystals");
    li_mass_total_key     =  model.Database().StorageKey("li mass total output");

    li_content_fluid_key = model.Database().StorageKey("li content fluid");
    li_content_liquid_key = model.Database().StorageKey("li content liquid");
    li_content_vapor_key = model.Database().StorageKey("li content vapor");
    li_content_crystals_key = model.Database().StorageKey("li content crystals");
    li_content_melt_key = model.Database().StorageKey("li content melt");

    li_concentration_fluid_key = model.Database().StorageKey("li concentration fluid");
    li_concentration_liquid_key = model.Database().StorageKey("li concentration liquid");
    li_concentration_vapor_key = model.Database().StorageKey("li concentration vapor");
    li_concentration_melt_key = model.Database().StorageKey("li concentration melt");
    li_concentration_crystals_key = model.Database().StorageKey("li concentration crystals");

    bfm_mass_out_total_key = model.Database().StorageKey("total bfm mass out");
    li_mass_out_total_key = model.Database().StorageKey("li mass out total");

    massfrac_li_liquid_key  = model.Database().StorageKey("li fraction liquid");
    massfrac_li_vapor_key = model.Database().StorageKey("li fraction vapor");
    massfrac_li_fluid_key = model.Database().StorageKey("li fraction fluid");
    massfrac_li_fluid_rescaled_key =
      model.Database().StorageKey("li fraction fluid rescaled");

    mob_li_liquid_key  = model.Database().StorageKey("liquid li mobility");
    mob_li_vapor_key = model.Database().StorageKey("vapor li mobility");

    rho_l_key = model.Database().StorageKey("density liquid");
    rho_v_key = model.Database().StorageKey("density vapor");

    bfm_key = model.Database().StorageKey("boundary flow mass");
    pore_volume_key = model.Database().StorageKey("pore volume");
    previous_pore_volume_key = model.Database().StorageKey("previous pore volume");
    bulk_volume_key = model.Database().StorageKey("bulk volume");

    mml_key = model.Database().StorageKey("liquid mass mobility"); //
    mmv_key = model.Database().StorageKey("vapor mass mobility"); //

    ml_key = model.Database().StorageKey("fluid mass liquid"); //
    mv_key = model.Database().StorageKey("fluid mass vapor"); //
    mt_key = model.Database().StorageKey("fluid density"); //

    temperature_key = model.Database().StorageKey("temperature");
    fluid_state_key = model.Database().StorageKey("fluid state");

    fluid_pressure_key = model.Database().StorageKey("fluid pressure");
    lithostatic_pressure_key = model.Database().StorageKey("lithostatic pressure");


    meltmass_key = model.Database().StorageKey("melt mass");
    crystalmass_key = model.Database().StorageKey("crystal mass");

    D_li_vl_key =
      model.Database().StorageKey("partition coefficient lithium vapor liquid node");

  }

  template<uint32_t dim>
  LithiumVisitor<dim>::~LithiumVisitor()
  {
  }

  template <uint32_t dim>
  void LithiumVisitor<dim>::Visit(Node<dim> *n)
  {
  }

  // test
  template<uint32_t dim>
  void LithiumVisitor<dim>::CalculateMassFractionsAndMobilities(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(li_content_liquid_key, li_content_liquid);
        (*nit)->Read(li_content_vapor_key, li_content_vapor);
        (*nit)->Read(pore_volume_key, pore_volume);
        (*nit)->Read(previous_pore_volume_key, previous_pore_volume);


        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(mml_key, mml );
        (*nit)->Read(mmv_key, mmv );

        if (ml() != 0.)
          {
            massfrac_li_liquid()  =   li_content_liquid() * /*previous_*/pore_volume() /
                                      (ml() * pore_volume());
          }
        else
          {
            massfrac_li_liquid()  =   0.;
          }

        if (mv() != 0.)
          {
            massfrac_li_vapor()  =   li_content_vapor() * /*previous_*/pore_volume() /
                                     (mv() * pore_volume());
          }
        else
          {
            massfrac_li_vapor()  =   0.;
          }

        massfrac_li_fluid() = (li_content_liquid() + li_content_vapor())
                              * /*previous_*/pore_volume() / ((ml() + mv()) * pore_volume());

        mob_li_liquid() = mml() * massfrac_li_liquid();
        mob_li_vapor() = mmv() * massfrac_li_vapor();

        (*nit)->Store(massfrac_li_liquid_key, massfrac_li_liquid);
        (*nit)->Store(massfrac_li_vapor_key, massfrac_li_vapor);
        (*nit)->Store(massfrac_li_fluid_key, massfrac_li_fluid);
        (*nit)->Store(mob_li_liquid_key, mob_li_liquid);
        (*nit)->Store(mob_li_vapor_key, mob_li_vapor);

      }
  }

  template<uint32_t dim>
  void LithiumVisitor<dim>::InitializeLithiumVisitor(Model<dim> *model)
  {
    const Region<dim>   &ref1 = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref1.NodesBegin();
          nit != ref1.NodesEnd();
          ++nit )
      {
        (*nit)->Read(pore_volume_key, pore_volume);

        li_content_liquid() = 0.;
        li_content_vapor() = 0.;
        li_content_melt() = 0;
        li_content_crystals() = 0.;
        li_content_fluid() = 0.;

        li_concentration_liquid() = 0.;
        li_concentration_vapor() = 0.;
        li_concentration_melt() = 0;
        li_concentration_crystals() = 0.;
        li_concentration_fluid() = 0.;

        // to check mass balance
        li_mass_liquid() = 0.;
        li_mass_vapor() = 0.;
        li_mass_melt() = 0.;
        li_mass_crystals() = 0.;
        li_mass_total_output() = 0.;

        D_li_vl() = 0.; // TS

        (*nit)->Store(li_content_liquid_key, li_content_liquid);
        (*nit)->Store(li_content_vapor_key, li_content_vapor);
        (*nit)->Store(li_content_fluid_key, li_content_fluid);
        (*nit)->Store(li_content_melt_key, li_content_melt);
        (*nit)->Store(li_content_crystals_key, li_content_crystals);

        (*nit)->Store(li_concentration_liquid_key, li_concentration_liquid);
        (*nit)->Store(li_concentration_vapor_key, li_concentration_vapor);
        (*nit)->Store(li_concentration_fluid_key, li_concentration_fluid);
        (*nit)->Store(li_concentration_melt_key, li_concentration_melt);
        (*nit)->Store(li_concentration_crystals_key, li_concentration_crystals);

        (*nit)->Store(previous_pore_volume_key, pore_volume);
        (*nit)->Store(D_li_vl_key, D_li_vl); // TS
      }

    const Region<dim>  &ref2 = model_ref_.Region("CHAMBER");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref2.NodesBegin();
          nit != ref2.NodesEnd();
          ++nit )
      {
        (*nit)->Read(meltmass_key, meltmass);
        (*nit)->Read(pore_volume_key, pore_volume);
        (*nit)->Read(bulk_volume_key, bulk_volume);

        li_content_melt() = li_concentration_melt_initialize * meltmass() / bulk_volume();
        li_concentration_melt() = 1.e6 * li_concentration_melt_initialize;

        (*nit)->Store(li_content_melt_key, li_content_melt);
        (*nit)->Store(li_concentration_melt_key, li_concentration_melt);


      }

    li_mass_out_total() = 0.;
    bfm_mass_out_total() = 0.;

    model_ref_.Store( li_mass_out_total_key, li_mass_out_total );
    model_ref_.Store( bfm_mass_out_total_key, bfm_mass_out_total );
  }


  template<uint32_t dim> // reads the constant partioning coefficients for single phase fluid-melt partitioning and xtel-melt partitioning
  void LithiumVisitor<dim>::ReadInputArguments(double
                                               ext_partition_coefficient_lithium_fluid_melt,
                                               double ext_partition_coefficient_lithium_xtel_melt,
                                               double ext_partition_coefficient_lithium_vl_slope,
                                               double ext_li_concentration_melt_initialize)

  {
    D_li_fm                          = ext_partition_coefficient_lithium_fluid_melt;
    D_li_cm                          = ext_partition_coefficient_lithium_xtel_melt;
    D_li_vl_slope                    = ext_partition_coefficient_lithium_vl_slope;
    li_concentration_melt_initialize = ext_li_concentration_melt_initialize;
  }


  template<uint32_t dim>
  double LithiumVisitor<dim>::ReadPartitionCoefficient(double lv_partition_slope,
                                                       double rho_liq, double rho_vap) // e.g. slope = 2
  {
    //input arguments: log(K v/l) at log(rhov/rhol)=-1, liquid density, vapor density
    double log_ratio_rhov_rhol(0.);

    if (rho_liq != 0.)
      {
        log_ratio_rhov_rhol = log10(rho_vap / rho_liq);
      }

    return pow(10., log_ratio_rhov_rhol * lv_partition_slope); // returns D_li_vl
  }


  template<uint32_t dim>
  void LithiumVisitor<dim>::Partitioning(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(li_content_liquid_key, li_content_liquid);
        (*nit)->Read(li_content_vapor_key, li_content_vapor);
        (*nit)->Read(li_content_fluid_key, li_content_fluid);
        (*nit)->Read(li_content_melt_key, li_content_melt);
        (*nit)->Read(li_content_crystals_key, li_content_crystals);

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(rho_l_key, rho_l);
        (*nit)->Read(rho_v_key, rho_v);

        (*nit)->Read(meltmass_key, meltmass );
        (*nit)->Read(crystalmass_key, crystalmass );
        (*nit)->Read(pore_volume_key, pore_volume );
        (*nit)->Read(previous_pore_volume_key, previous_pore_volume );
        (*nit)->Read(bulk_volume_key, bulk_volume );
        (*nit)->Read(fluid_state_key, fluid_state );

        //initialize Li masses in phases
        double li_mass_tot(0.);

        li_mass_liquid() = li_content_liquid()*/*previous_*/pore_volume();
        li_mass_vapor() = li_content_vapor()*/*previous_*/pore_volume();
        li_mass_melt() = li_content_melt() * bulk_volume();
        li_mass_crystals() = li_content_crystals() * bulk_volume();
        li_mass_fluid() = li_mass_liquid() + li_mass_vapor();

        //if there is melt present, Li mass in crystals and in melt are taken into account
        if (!essentiallyEqual(meltmass(), 0., numeric_limits<double>::epsilon()))
          {
            li_mass_tot = li_mass_fluid() + li_mass_melt() + li_mass_crystals(); //kg/m^3
          }
        //if there is no melt present, only Li mass in fluid is taken into account
        else
          {
            li_mass_tot = li_mass_fluid(); //kg/m^3
          }

        double massfractionwaterliquid = ml() / (ml() + mv());
        double rho_vap = rho_v();
        double rho_liq = rho_l();





        //if node is in magma chamber (meltmass non-zero)
        if (!essentiallyEqual(meltmass(), 0., numeric_limits<double>::epsilon()))
          {

            li_mass_melt() = li_mass_tot / (1.
                                            + D_li_fm * (ml() + mv()) / (meltmass() / pore_volume())
                                            + D_li_cm * crystalmass() / meltmass());
            li_mass_fluid() = li_mass_melt() * D_li_fm * (ml() + mv()) / (meltmass() / pore_volume());
            li_mass_crystals() = li_mass_melt() * D_li_cm * crystalmass() / meltmass();

          }
        //if there is no melt present
        else
          {
            li_mass_melt()  = 0.;
          }

        //partitioning between the aqueous phases
        if (mv() == 0. || ml() == 0.)
          {
            li_mass_liquid() = massfractionwaterliquid * li_mass_fluid();
            li_mass_vapor() = li_mass_fluid() - li_mass_liquid();
          }
        else
          {
            D_li_vl = ReadPartitionCoefficient(D_li_vl_slope, rho_liq,
                                               rho_vap); //input arguments: log(K v/l) at log(rhov/rhol)=-1, liquid density, vapor density
            li_mass_liquid() = li_mass_fluid() / (1 + D_li_vl() * mv() / ml());
            li_mass_vapor() = li_mass_fluid() - li_mass_liquid();
          }

        li_mass_fluid() = li_mass_liquid() + li_mass_vapor();

        li_content_liquid() = li_mass_liquid() / pore_volume();
        li_content_vapor() = li_mass_vapor() / pore_volume();
        li_content_fluid() = li_content_liquid() + li_content_vapor();
        li_content_melt() = li_mass_melt() / bulk_volume();
        li_content_crystals() = li_mass_crystals() / bulk_volume();

        //calculate concentrations
        //all concentrations in ppm (mg Li/kg fluid or melt)
        /////
        if (ml() != 0.)
          {
            li_concentration_liquid() = 1.e6 * li_content_liquid() / ml();
          }
        else
          {
            li_concentration_liquid() = 0.;
          }

        /////
        if (mv() != 0.)
          {
            li_concentration_vapor() = 1.e6 * li_content_vapor() / mv();
          }
        else
          {
            li_concentration_vapor() = 0.;
          }

        /////
        if (!essentiallyEqual(meltmass(), 0., numeric_limits<double>::epsilon()))
          {
            li_concentration_melt() = 1.e6 * li_mass_melt() / meltmass();
          }
        else
          {
            li_concentration_melt() = 0.;
          }

        /////
        if (!essentiallyEqual(crystalmass(), 0., numeric_limits<double>::epsilon()))
          {
            li_concentration_crystals() = 1.e6 * li_content_crystals() * bulk_volume() /
                                          crystalmass();

          }
        else
          {
            li_concentration_crystals() = 0.;
          }

        /////
        li_concentration_fluid() = 1.e6 * (li_content_liquid() + li_content_vapor()) /
                                   (mv() + ml());

        li_mass_total_output() = li_mass_fluid() + li_mass_melt() + li_mass_crystals();


        (*nit)->Store(li_concentration_liquid_key, li_concentration_liquid);
        (*nit)->Store(li_concentration_vapor_key, li_concentration_vapor);
        (*nit)->Store(li_concentration_melt_key, li_concentration_melt);
        (*nit)->Store(li_concentration_fluid_key, li_concentration_fluid);
        (*nit)->Store(li_concentration_crystals_key, li_concentration_crystals);

        (*nit)->Store(li_content_melt_key, li_content_melt);
        (*nit)->Store(li_content_liquid_key, li_content_liquid);
        (*nit)->Store(li_content_vapor_key, li_content_vapor);
        (*nit)->Store(li_content_fluid_key, li_content_fluid);
        (*nit)->Store(li_content_crystals_key, li_content_crystals);
        (*nit)->Store(D_li_vl_key, D_li_vl);

        (*nit)->Store(li_mass_liquid_key, li_mass_liquid);
        (*nit)->Store(li_mass_vapor_key, li_mass_vapor);
        (*nit)->Store(li_mass_fluid_key, li_mass_fluid);
        (*nit)->Store(li_mass_melt_key, li_mass_melt);
        (*nit)->Store(li_mass_crystals_key, li_mass_crystals);
        (*nit)->Store(li_mass_total_key, li_mass_total_output);

      }
  }




  template<uint32_t dim>
  void LithiumVisitor<dim>::ApplyBoundaryFlowCorrections(Model<dim> *model)
  {
    const Region<dim>   &ref = model_ref_.Region("Model");

    model_ref_.Read( li_mass_out_total_key, li_mass_out_total );
    model_ref_.Read( bfm_mass_out_total_key, bfm_mass_out_total );

    double totmassout(bfm_mass_out_total()),
           totliout(li_mass_out_total());


    for ( typename std::vector<Node<dim>*>::const_iterator
          nit = ref.NodesBegin();
          nit != ref.NodesEnd();
          ++nit )
      {
        (*nit)->Read(bfm_key, bfm);

        (*nit)->Read(ml_key, ml);
        (*nit)->Read(mv_key, mv);

        (*nit)->Read(li_content_liquid_key, li_content_liquid);
        (*nit)->Read(li_content_vapor_key, li_content_vapor);

        (*nit)->Read(mml_key, mml);
        (*nit)->Read(mmv_key, mmv);

        double scalingfactor(0.);
        double li_frac_fluid_inflow(0.); //must be adapted if recharge water has lithium content

        double li_mass_before = pore_volume() * (li_content_liquid() + li_content_vapor());


        if (bfm() > 0.) //outflow
          {
            //            scalingfactor = 1./(1.+(bfm()/(mv()+ml())));

            scalingfactor = 1. / (1. + ((mml() / (mmv() + mml())) * bfm() / (mv() + ml())));

            //NK may 2024, : the relative mass proportions of vapor and liquid contributing to the total boundary flow should be approximately proportional
            //to the relative proportions of mobilities. Since the previous mobilities are not accessible, we take the current mobilities
            //as an approximation. Therefore the contribution of the liquid to the total boundary flow is (mmv/(mmv+mml))*bfm,
            //analogous for vapor

            li_content_liquid() *= scalingfactor;

            scalingfactor = 1. / (1. + ((mmv() / (mmv() + mml())) * bfm() / (mv() + ml())));
            li_content_vapor() *= scalingfactor;
          }

        else   //inflow
          {
            li_content_liquid() -= bfm() * li_frac_fluid_inflow;
            li_content_vapor() -= bfm() * li_frac_fluid_inflow;
          }

        li_content_fluid() = li_content_liquid() + li_content_vapor();

        totliout += li_mass_before - pore_volume() * li_content_fluid();
        totmassout += bfm() * pore_volume();

        (*nit)->Store(li_content_liquid_key, li_content_liquid);
        (*nit)->Store(li_content_vapor_key, li_content_vapor);
        (*nit)->Store(li_content_fluid_key, li_content_fluid);
      }

    li_mass_out_total() += totliout;
    bfm_mass_out_total() += totmassout;

    model_ref_.Store( li_mass_out_total_key, li_mass_out_total );
    model_ref_.Store( bfm_mass_out_total_key, bfm_mass_out_total );

  }

  template<uint32_t dim>
  void LithiumVisitor<dim>::Visit( Region<dim> *r )
  {}

  template<uint32_t dim>
  void LithiumVisitor<dim>::Visit( Model<dim> *m )
  {}

  template<uint32_t dim>
  void LithiumVisitor<dim>::Visit( Element<dim> *e )
  {}




  template class LithiumVisitor<1>;
  template class LithiumVisitor<2>;
  template class LithiumVisitor<3>;

}// csmp
