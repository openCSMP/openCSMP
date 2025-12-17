#include "AuxillaryFunctions.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Exception.h"

using namespace std;

namespace csmp
{

template<uint32_t dim>
void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<dim>& sg,
                                                              const char* element_property,
                                                              const char* node_property,
                                                              const char* region,
                                                              double multiplication_factor )
{

    const csmp::Index  e_key = sg.Database().StorageKey(element_property);
    if ( e_key.place != ELEMENT || e_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "volumeWeightedDistributionOfElementPropertyToNodeAndAdd",
                               element_property, "variable must be a scalar placed on the element." );

    const csmp::Index  n_key = sg.Database().StorageKey(node_property);
    if ( n_key.place != NODE || n_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "volumeWeightedDistributionOfElementPropertyToNodeAndAdd",
                               node_property, "variable must be a scalar placed on the node." );

    for ( auto eit=sg.Region(region).CellsBegin();
          eit!=sg.Region(region).CellsEnd(); eit++ )
    {
        uint32_t number_of_nodes((*eit)->Nodes());
        double e_property_value((*eit)->Read(e_key));
        e_property_value *= (*eit)->Volume();
        e_property_value *= multiplication_factor;

        for ( uint32_t i = 0U; i < number_of_nodes; i++ )
        {
            double n_property_value( ((*eit)->N(i))->Read(n_key) );
            n_property_value += ( e_property_value / number_of_nodes );
            ((*eit)->N(i))->Store( n_key, ScalarVariable( PLAIN, n_property_value )  );
        }
    }

} // end volumeWeightedDistributionOfElementPropertyToNodeAndAdd

template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<1U>&, const char*, const char*, const char*, double );
template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<2U>&, const char*, const char*, const char*, double );
template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<3U>&, const char*, const char*, const char*, double );




template<uint32_t dim>
void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<dim>& sg,
                                                              const char* element_property,
                                                              const char* node_property,
                                                              const char* volume_modifier,
                                                              const char* region,
                                                              double multiplication_factor )
{

    const csmp::Index  e_key = sg.Database().StorageKey(element_property);
    if ( e_key.place != ELEMENT || e_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "volumeWeightedDistributionOfElementPropertyToNodeAndAdd",
                               element_property, "variable must be a scalar placed on the element." );

    const csmp::Index  n_key = sg.Database().StorageKey(node_property);
    if ( n_key.place != NODE || n_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "volumeWeightedDistributionOfElementPropertyToNodeAndAdd",
                               node_property, "variable must be a scalar placed on the node." );

    const csmp::Index  vol_mod_key = sg.Database().StorageKey(volume_modifier);
    if ( vol_mod_key.place != ELEMENT || vol_mod_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "volumeWeightedDistributionOfElementPropertyToNodeAndAdd",
                               volume_modifier, "variable must be a scalar placed on the element." );

    for ( auto eit=sg.Region(region).CellsBegin();
          eit!=sg.Region(region).CellsEnd(); eit++ )
    {
        uint32_t number_of_nodes((*eit)->Nodes());
        double e_property_value((*eit)->Read(e_key));
        double vol_mod((*eit)->Read(vol_mod_key));
        e_property_value *= (*eit)->Volume();
        e_property_value *= vol_mod;
        e_property_value *= multiplication_factor;

        for ( uint32_t i{0U}; i < number_of_nodes; i++ )
        {
            double n_property_value( ((*eit)->N(i))->Read(n_key) );
            n_property_value += ( e_property_value / number_of_nodes );
            ((*eit)->N(i))->Store( n_key, ScalarVariable( PLAIN, n_property_value )  );
        }
    }

} // end volumeWeightedDistributionOfElementPropertyToNodeAndAdd

template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<1U>&, const char*, const char*, const char*, const char*, double );
template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<2U>&, const char*, const char*, const char*, const char*, double );
template void volumeWeightedDistributionOfElementPropertyToNodeAndAdd( Model<3U>&, const char*, const char*, const char*, const char*, double );




template<uint32_t dim>
void areaWeightedDistributionOfNodePropertyToNodeAndAdd(Model<dim>& sg,
                                                        const char* nodal_flux,
                                                        const char* nodal_value,
                                                        const char* region,
                                                        double multiplication_factor)
{
    const csmp::Index  nf_key = sg.Database().StorageKey(nodal_flux);
    if ( nf_key.place != NODE || nf_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "areaWeightedDistributionOfNodePropertyToNodeAndAdd",
                               nodal_flux, "variable must be a scalar placed on the node." );

    const csmp::Index  nv_key = sg.Database().StorageKey(nodal_value);
    if ( nv_key.place != NODE || nv_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "areaWeightedDistributionOfNodePropertyToNodeAndAdd",
                               nodal_value, "variable must be a scalar placed on the node." );

    for ( auto eit=sg.Region(region).PerimeterCellsBegin();
          eit!=sg.Region(region).CellsEnd(); eit++ )
    {
        for ( auto n{0U}; n < (*eit)->Nodes(); n++)
        {
            if ((*eit)->N(n)->AtBoundary() != NOT)
            {
                const Point<dim> p1((*eit)->N(n)->Coordinate());
                const double nf1((*eit)->N(n)->Read(nf_key));
                const double nv1((*eit)->N(n)->Read(nv_key));
                const uint32_t right_node(n == (*eit)->Nodes() - 1U ? 0U : n + 1U);
                const uint32_t left_node(n == 0U ? (*eit)->Nodes() - 1U : n - 1U);                
                double new_nv(0.);
                
                if ((*eit)->N(right_node)->AtBoundary() != NOT)
                {
                    const Point<dim> p2((*eit)->N(right_node)->Coordinate());
                    const double nf2((*eit)->N(right_node)->Read(nf_key));
                    new_nv += 0.5 * (nf1 + nf2) * 0.5 * p2.DistanceTo(p1);

                }

                if ((*eit)->N(left_node)->AtBoundary() != NOT)
                {
                    const Point<dim> p2((*eit)->N(left_node)->Coordinate());
                    const double nf2((*eit)->N(left_node)->Read(nf_key));
                    new_nv += 0.5 * (nf1 + nf2) * 0.5 * p2.DistanceTo(p1);

                }
                new_nv *= multiplication_factor;
                new_nv += nv1;
                (*eit)->N(n)->Store(nv_key, ScalarVariable(PLAIN, new_nv));            
            }
        }

    }

}

template void areaWeightedDistributionOfNodePropertyToNodeAndAdd(Model<1U>&, const char* ,const char* ,const char* ,double);
template void areaWeightedDistributionOfNodePropertyToNodeAndAdd(Model<2U>&, const char* ,const char* ,const char* ,double);
template void areaWeightedDistributionOfNodePropertyToNodeAndAdd(Model<3U>&, const char* ,const char* ,const char* ,double);


template<uint32_t dim>
void setPropertyToZero( Model<dim>& sg, const char* property )
{

    const csmp::Index  p_key = sg.Database().StorageKey(property);
    if ( ( p_key.place != ELEMENT &&  p_key.place != NODE ) || p_key.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "setPropertyToZero",
                               property, "variable must be a scalar placed on the element or node." );

    sg.InputPropertyValue( property, ScalarVariable( PLAIN, 0. ) );

}



template void setPropertyToZero( Model<1U>&, const char* );
template void setPropertyToZero( Model<2U>&, const char* );
template void setPropertyToZero( Model<3U>&, const char* );

} // end csmp
