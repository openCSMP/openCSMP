#include "GravityProjectionVisitor.h"
#include "Region.h"

using namespace std;

namespace csmp{

template<uint32_t dim>
GravityProjectionVisitor<dim>::GravityProjectionVisitor(
                                                        Model<dim>& m,
                                                        const Index& prop_idx,
                                                        const Index& result_idx,
                                                        VectorVariable<dim> vec )
    : Visitor<dim>( MODEL, ELEMENT ),
      prop_idx_   ( prop_idx),
      result_idx_ ( result_idx),
      const_vec_  ( ( prop_idx==Index() ? true: false) )
{
    Index nan;

    if( result_idx == nan )
    {
        std::string result_name;

        if( prop_idx!=nan )
            result_name = m.Database().Name( prop_idx );
        else
            result_name = "vector field";

        result_name +=" low dimensional projection";

        if ( !m.Database().IsDefined( result_name.c_str() ) )
            m.CreateProperty( result_name.c_str(),  "SI", VECTOR, ELEMENT );

        result_idx_  = m.Database().StorageKey( result_name.c_str() );
    }

    debug_ = false;

    if( const_vec_ )
        vec_ = vec;
}

template<uint32_t dim>
GravityProjectionVisitor<dim>::~GravityProjectionVisitor()
{

}

template<uint32_t dim>
Index GravityProjectionVisitor<dim>::Get_ResultIndex( )
{
    return result_idx_;
}

template<uint32_t dim>
Index GravityProjectionVisitor<dim>::Get_PropertyIndex( )
{
    return prop_idx_;
}

template<uint32_t dim>
void GravityProjectionVisitor<dim>::Get_Result( Element<dim>* eptr, VectorVariable<dim>& result )
{
    eptr->Read( result_idx_, result );
}

template<uint32_t dim>
void GravityProjectionVisitor<dim>::Visit( Model<dim>* m ){
    csmp::Region<dim>&  mref(m->Region("Model"));
    for ( auto e_it=mref.CellsBegin(); e_it!=mref.CellsEnd(); e_it++ )
        (*e_it)->Accept( *this );
}

template<uint32_t dim>
void GravityProjectionVisitor<dim>::Visit(Region<dim>* region ){
    for ( auto e_it=region->CellsBegin(); e_it!=region->CellsEnd(); e_it++ )
        (*e_it)->Accept( *this );
}

template<uint32_t dim>
void GravityProjectionVisitor<dim>::Visit( Element<dim>* eptr )
{
    if( dim == 1U )
    {
        if(!const_vec_ )
            eptr->Read(prop_idx_,vec_);

        //! projection of a vector field in 1D is simply the original fields

        proj_(0) = vec_[0];

        if( debug_ )
            cout<<"\n Element["<<eptr->Idx()<<"]: vec_proj[0]="<<proj_[0]<<endl;

        eptr->Store( result_idx_, proj_);

    }
    else if( dim== 2U )
    {
        if(!const_vec_ )
            eptr->Read(prop_idx_,vec_);

        if( eptr->FE()->IsLineElement() )
        {
            //! projection of a vector field onto a unit norm vector r(r[0],r[1]),|r|=1
            //! VecProj = ( vec, r ) r

            Point<dim> line_vector( eptr->N(0)->Coordinate() - eptr->N(1)->Coordinate() );

            line_vector.NormalizeLengthTo(1.);

            line_vector *= ( line_vector[0]*vec_[0] + line_vector[1]*vec_[1] );

            //! ( VecProj, vec ) > 0 Note that the direction of the result vector should be align with the original one
            if( line_vector[0]*vec_[0] + line_vector[1]*vec_[1] < 0.0 )
                line_vector *= -1.;

            proj_(0) = line_vector[0];
            proj_(1) = line_vector[1];
        }
        else
        {
            proj_(0) = vec_(0);
            proj_(1) = vec_(1);
        }


        if( debug_ )
            cout<<"\n Element["<<eptr->Idx()<<"]: vec_proj[0]="<<proj_[0]<<" ; vec_proj[1]="<<proj_[1]<<endl;

        eptr->Store( result_idx_, proj_);

    }else{


        if(!const_vec_ )
            eptr->Read(prop_idx_,vec_);

        if( eptr->FE()->IsLineElement() )
        {
            //! projection of a vector field onto a unit norm vector r(r[0],r[1],r[2]),|r|=1
            //! VecProj = ( vec, r ) r

            Point<dim> line_vector( eptr->N(0)->Coordinate() - eptr->N(1)->Coordinate() );

            line_vector.NormalizeLengthTo(1.);

            line_vector *= ( line_vector[0]*vec_[0] + line_vector[1]*vec_[1] + line_vector[2]*vec_[2] );

            //! (n,vec) >0 Note that the direction of the result vector should be align with the original one
            if( line_vector[0]*vec_[0] + line_vector[1]*vec_[1] + line_vector[2]*vec_[2] < 0.0 )
                line_vector *= -1.;

            proj_(0) = line_vector[0];
            proj_(1) = line_vector[1];
            proj_(2) = line_vector[2];
        }
        else if( eptr->FE()->IsSurfaceElement())
        {

            //! projection of a vector field onto a surface with a unit normal vector n(n[0], n[1], n[2]), |n|=1
            //! VecProj = n x ( n x vec ) = ( vec, n ) n - vec

            Point<dim> line_vector1( eptr->N(1)->Coordinate() - eptr->N(0)->Coordinate() );
            Point<dim> line_vector2( eptr->N(2)->Coordinate() - eptr->N(0)->Coordinate() );

            Point<dim> surface_normal( crossProduct( line_vector1,line_vector2) );

            surface_normal.NormalizeLengthTo (1.);

            surface_normal *= ( surface_normal[0]*vec_[0] + surface_normal[1]*vec_[1] + surface_normal[2]*vec_[2] );

            surface_normal[0] -= vec_[0];
            surface_normal[1] -= vec_[1];
            surface_normal[2] -= vec_[2];


            //! ( VecProj, vec ) > 0 Note that the direction of the result vector should be align with the original one
            if( surface_normal[0]*vec_[0] + surface_normal[1]*vec_[1] + surface_normal[2]*vec_[2] < 0.0 )
                surface_normal *= -1.;

            proj_(0) = surface_normal[0];
            proj_(1) = surface_normal[1];
            proj_(2) = surface_normal[2];
        }
        else
        {
            proj_(0) = vec_[0];
            proj_(1) = vec_[1];
            proj_(2) = vec_[2];
        }

        if( debug_ )
            cout<<"\n Element["<<eptr->Idx()<<"]: vec_proj[0]="<<proj_[0]<<" ; vec_proj[1]="<<proj_[1]<<" ; vec_proj[2]="<<proj_[2]<<endl;

        eptr->Store( result_idx_, proj_);


    }

} // end Visit

template class GravityProjectionVisitor<1U>;
template class GravityProjectionVisitor<2U>;
template class GravityProjectionVisitor<3U>;

} //csmp
