// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SplitRegionDiffusiveHeating.h"
#include "Exception.h"
#include "Model.h"

using namespace std;
namespace csmp {

template<uint32_t dim>
SplitRegionDiffusiveHeating<dim>::SplitRegionDiffusiveHeating( Model<dim>& model, const std::string& sb_name)
    :
    // 1D leakage computation
    // Only works with split boundaries

    sb_name     (sb_name),

    model_ref_                    ( model ),
    prop_ref_                     ( model.Database() ),
    sb                            ( model.SplitBoundary( sb_name ) ),

    q_key( model.Database().StorageKey("split region diffused heat source") ),
    t_key( model.Database().StorageKey("temperature") ),
    nA_key( model.Database().StorageKey("split region nodal area") ),

    kth_key( model.Database().StorageKey("thermal conductivity") ),
    thick_key( model.Database().StorageKey("thickness") )

{
    ComputeNodalArea();
} // end SplitRegionDiffusiveHeating

template<uint32_t dim>
SplitRegionDiffusiveHeating<dim>::~SplitRegionDiffusiveHeating()
{} // end ~SplitRegionDiffusiveHeating

// The main task:
template<uint32_t dim>
void SplitRegionDiffusiveHeating<dim>::Apply()
{
    for (auto sb_c : sb.CellVector() )
    {

        kth = sb_c->InterveningElement()->Read(kth_key);
        thick = sb_c->InterveningElement()->Read(thick_key);

        //test
        //thick()=0.;

        uint32_t n_nodes = sb_c->FE()->Nodes();

        for ( uint32_t n{0U}; n<n_nodes;++n)
        {
            in_node  = sb_c->MatchingN(n,INSIDE);
            out_node  = sb_c->MatchingN(n,OUTSIDE);
            mid_node  = sb_c->MatchingN(n,MIDDLE);

            Point<dim> in_mid = in_node->Coordinate() - mid_node->Coordinate();
            Point<dim> out_mid = out_node->Coordinate() - mid_node->Coordinate();

            //Point<dim> unrml( sb_c->UnitNormal(MIDDLE) );
            double distance_in_mid = in_mid.Length();
            double distance_out_mid = out_mid.Length();

            // cerr<<endl<<"distance_in_mid: "<<distance_in_mid;
            // cerr<<endl<<"distance_out_mid: "<<distance_out_mid;
            // cerr<<endl<<distance_out_mid + thick()/2.;

            //on perimeter, we set an arbitrary distance of 0.5 meter
            //if( !(distance_in_mid >0.) ) distance_in_mid = 0.5;
            //if( !(distance_out_mid >0.) ) distance_out_mid = 0.5;

            t_in = in_node->Read(t_key);
            t_out = out_node->Read(t_key);
            t_mid = mid_node->Read(t_key);

            nA = mid_node->Read(nA_key);

            //flux calculation
            q_in() = -kth() * (t_in()-t_mid())/( distance_in_mid + thick()/2. )*nA();
            q_mid()= -kth() * (t_mid()-t_in())/( distance_in_mid + thick()/2. )*nA()
                      -kth() * (t_mid()-t_out())/( distance_out_mid + thick()/2. )*nA();
            q_out()= -kth() * (t_out()-t_mid())/( distance_out_mid + thick()/2. )*nA();

            in_node->Store(q_key, q_in);
            out_node->Store(q_key, q_out);
            mid_node->Store(q_key, q_mid);
        }
    }
}

template<uint32_t dim>
void SplitRegionDiffusiveHeating<dim>::ComputeNodalArea()
{

    for (auto sb_c : sb.CellVector() )
    {

        uint32_t n_nodes = sb_c->FE()->Nodes();

        double area3 = sb_c->InterveningElement()->Volume() / static_cast<double>( sb_c->InterveningElement()->Nodes() );

        cerr<<endl<<endl<<"area: "<<sb_c->InterveningElement()->Volume();
        cerr<<", num nodes: "<<sb_c->InterveningElement()->Nodes();

        for ( uint32_t n{0U}; n<n_nodes;++n)
        {
            in_node  = sb_c->MatchingN(n,INSIDE);
            out_node  = sb_c->MatchingN(n,OUTSIDE);
            mid_node  = sb_c->MatchingN(n,MIDDLE);

            nA = mid_node->Read(nA_key);

            nA()+=area3;

            in_node->Store(nA_key, nA);
            out_node->Store(nA_key, nA);
            mid_node->Store(nA_key, nA);

        }
    }
}

template class SplitRegionDiffusiveHeating<1U>;
template class SplitRegionDiffusiveHeating<2U>;
template class SplitRegionDiffusiveHeating<3U>;

} // csmp
