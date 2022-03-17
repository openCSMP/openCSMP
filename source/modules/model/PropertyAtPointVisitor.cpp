#include "PropertyAtPointVisitor.h"
#include "Element.h"

using namespace std;

namespace csmp {

/**

    @author Andrey Mezentsev (10/2004)
*/
template<uint32_t dim>
PropertyAtPointVisitor<dim>::PropertyAtPointVisitor( const Model<dim>& m,
                                                     const map <size_t, vector<double> > & inXYZ,
                                                     const char *propertyName )
    : pref_(m.Database()),
      maxElementsInTheMesh_(m.Region("Model").Elements()),  // Number of elements in the mesh
      TargetElement_        ( NULL ),                       // Found element
      CurrentElement_       ( NULL ),                       // Intermediate element
      NumberOfElementNodes_ ( 8 ),                          // Max Number of nodes
      numberOfPointsFound_  ( 0 ),                          // Number of currently found points
      numberOfPoints_       ( 0 ),                          // Number of points needed to be found
      xyz_                  ( inXYZ ),                      // Points of interest
      NI_                   ( NumberOfElementNodes_ ),      // Shape function values
      NPS_                  ( NumberOfElementNodes_ ),      // Scalar nodal properties
      NPV_                  ( NumberOfElementNodes_ ),      // Vector nodal properties
      NPT_                  ( NumberOfElementNodes_ ),      // Tensor nodal properties
      NPA_                  ( NumberOfElementNodes_ ),      // Array nodal properties
      NPFA_                 ( NumberOfElementNodes_ ),      // FlaggedArray nodal properties
      maxIterations_        ( 5 ),                          // Neighbour search
      currIteration_        ( 0 ),                          // Neighbour search
      precision_            ( 1.0e-8),                      // Precision for shape function values
      continueSearch_       ( true ),                       // Status of search algorithm
      bruteforce_           ( false ),                      // Set Brute Force methods by default
      debug_                ( false )                       // Debug mode
{

    Visitor<dim>::ApplicationLevel(MODEL);

    Visitor<dim>::ApplicationTarget(ELEMENT);

    SetPropertyKeys( propertyName );

    // Set targets to -1 index

    for( typename map <size_t,vector<double> >::iterator itX = xyz_.begin(); itX!=xyz_.end(); itX++)
    {
        numberOfPoints_++;
        elements_ids_found_[itX->first] = (size_t)-1;
        //elements_found_[itX->first] = NULL;
    }

    if( numberOfPointsFound_ == numberOfPoints_ )
    {
        if(debug_)
            cout<<"  PropertyAtPointVisitor<> All points found..."<<endl;

        SetContinueSearchOff();
    }

}


template<uint32_t dim>
PropertyAtPointVisitor<dim>::~PropertyAtPointVisitor()
{

}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::SetDebugOn()
{
    debug_ = true;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::SetContinueSearchOff()
{
    continueSearch_ = false;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::SetBruteForceOff()
{
    bruteforce_ = false;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::SetBruteForceOn()
{
    bruteforce_ = true;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::SetPropertyKeys( const char* prop )
{
    csmp::ErrorHandler& csmp_error(csmp::ErrorHandler::Instance());

    prop_idx_ = pref_.StorageKey(prop);

    if (prop_idx_.type != SCALAR && prop_idx_.type != VECTOR && prop_idx_.type != TENSOR )
    {
        csmp_error.notice( FATAL_ERROR, "PropertyAtPointVisitor<dim>::SetNodeAndElementPropertyKeys",
                           prop, "Element property must be a scalar, vector or tensor" );
    }
}

template<uint32_t dim>
bool PropertyAtPointVisitor<dim>::CheckResults( )
{
    bool success( ( numberOfPointsFound_ == numberOfPoints_ ) ? true: false );

    if( success )
        cout<<"  PropertyAtPointVisitor<dim> All points were found..."<<endl;
    else
      {
         cout<<" WARNING: PropertyAtPointVisitor<dim>: Unfortunately not all of the points were found!!!"<<endl;
         cout<<" There were "<<numberOfPointsFound_<<"  Points found from "<<numberOfPoints_<<endl;
      }
    return success;
}

/**
   returns NULL pointer if the element cannot be found
*/
template<uint32_t dim>
Element<dim>* PropertyAtPointVisitor<dim>::ElementThatContains( size_t point ) const
{
    typename map<size_t,Element<dim>*>::const_iterator it=elements_found_.find(point);
    if ( it == elements_found_.end() ) return NULL;
    return (*it).second;
}

template<uint32_t dim>
map<size_t,size_t>&  PropertyAtPointVisitor<dim>::TargetFound()
{
    return elements_ids_found_;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::PropertyValueAt( size_t pt, ScalarVariable& sc )

{
    typename std::map<size_t,ScalarVariable >::const_iterator it = propS_.find(pt);

    if ( it == propS_.end() )
        throw std::range_error("PropertyAtPointVisitor<dim>::PropertyValueAt(scalar)");

    sc = (*it).second;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::PropertyValueAt( size_t pt,
                                                          VectorVariable<dim>& vc )
{
    typename std::map<size_t,VectorVariable<dim> >::const_iterator it = propV_.find(pt);

    if ( it == propV_.end() )
        throw std::range_error("PropertyAtPointVisitor<dim>::PropertyValueAt(vector)");

    vc = (*it).second;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::PropertyValueAt( size_t pt,
                                                          TensorVariable<dim>& ts )
{
    typename std::map<size_t,TensorVariable<dim> >::const_iterator it = propT_.find(pt);

    if ( it == propT_.end() )
        throw std::range_error("PropertyAtPointVisitor<dim>::PropertyValueAt(tensor)");

    ts = (*it).second;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::PropertyValueAt( size_t pt,
                                                          ArrayVariable& av )

{
    typename std::map<size_t,ArrayVariable >::const_iterator it = propA_.find(pt);

    if ( it == propA_.end() )
        throw std::range_error("PropertyAtPointVisitor<dim>::PropertyValueAt(array)");

    av = (*it).second;
}

template<uint32_t dim>
void PropertyAtPointVisitor<dim>::PropertyValueAt( size_t pt,
                                                          FlaggedArrayVariable& fav )

{
    typename std::map<size_t,FlaggedArrayVariable >::const_iterator it = propFA_.find(pt);

    if ( it == propFA_.end() )
        throw std::range_error("PropertyAtPointVisitor<dim>::PropertyValueAt(flagged array)");

    fav = (*it).second;
}






template<uint32_t dim>
bool PropertyAtPointVisitor<dim>::isCloseToBarycenter( const vector<double> currXyz,
                                                       Element<dim>* e, VectorVariable<dim> bc )
{

    // Check whether point is lying inside the neghborghood of the element BaryCenter
	e->CoordinateMatrix( );

    double distance(0.0);
    double max_distance(0.0);
    for ( size_t i=0; i<e->Nodes(); i++ )
    {
        distance = 0.0;
        for ( auto j=0; j<dim; j++ )
            distance += (e->FE()->XYZ(i,j)- bc[j])*(e->FE()->XYZ(i,j)- bc[j]);
        distance = sqrt (distance );
        max_distance = std::max( max_distance, distance);
    }

    distance = 0.0;
    for ( auto j=0; j<dim; j++ )
        distance += (currXyz[j]- bc[j])*(currXyz[j]- bc[j]);
    distance = sqrt (distance );

    double distance_factor( 1.1 );

    if( distance <= distance_factor*max_distance )
        return true;

    return false;
}


template<uint32_t dim>
bool PropertyAtPointVisitor<dim>::FindPoint_BruteForceSearch( const vector<double> currXyz,
                                                              Element<dim>* e,
                                                              VectorVariable<dim> bc )
{
    if( isCloseToBarycenter(currXyz,e,bc) )
    {
        //cout<<"Visiting Element["<<e->Idx()<<"] of type: "<< parseFiniteElementType( e->FE_Type() ) << endl;

        e->CoordinateMatrix( );
        e->FE()->N(NI_,currXyz);

        typename vector<double>::const_iterator Nmin = min_element( NI_.begin(), NI_.end() );
        typename vector<double>::const_iterator Nmax = max_element( NI_.begin(), NI_.end() );

        if (debug_)
        {
            if( dim == 2 )
                std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<"; z="<<currXyz[2]<<std::endl;
            else
                std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<std::endl;

            std::cout<<"Nmax="<<*Nmax<<"; Nmin="<<*Nmin<<std::endl;
        }

        // found element, containing given point
        if((*Nmax<=(1.0+precision_)) && (*Nmin>=(0.0-precision_)))
        {
            TargetElement_  = e;

            if(debug_)
            {
                if( dim == 2 )
                    std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<"; z="<<currXyz[2]<<std::endl;
                else
                    std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<std::endl;

                vector<double> xyz(dim,0.0);
                for ( auto i=0; i<e->Nodes(); i++ ) {
                    for ( auto j=0; j<dim; j++ )
                        xyz[j] += NI_[i] * e->FE()->XYZ(i,j);
                }

                cout<<" Found point:("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")  vs.("
                    <<currXyz[0]<<","<<currXyz[1]<<","<<currXyz[2]<<")"<<endl;
            }

            return true;
        }

    }

    //cout<<"Finish Visiting Element["<<e->Idx()<<"] of type: "<< parseFiniteElementType( e->FE_Type() ) << endl;

    return false;

}






template<uint32_t dim>
bool PropertyAtPointVisitor<dim>::FindPoint_NeighborSearch( const vector<double> currXyz,Element<dim>* e, VectorVariable<dim> bc )
{

    if( isCloseToBarycenter(currXyz,e,bc) )
    {
        csmp::ErrorHandler& csmp_error(csmp::ErrorHandler::Instance());
        //
        //go through the neighbouring elements of element e unless target is found
        //
        set<uint32_t> listOfCheckedElements;  // set of previously checked elements for a given point
        typename set<uint32_t>::iterator itChd;
        VectorVariable<dim> current_bc;
        typename vector<double>::const_iterator Nmin;
        typename vector<double>::const_iterator Nmax;
        double minimumNi;

        bool stopFlag   = false;
        bool Find       = false;
        bool checked    = false;

        currIteration_ = 0;
        CurrentElement_ = e;
        TargetElement_  = e;

        do
        {
            if( TargetElement_ != NULL )
                CurrentElement_ = TargetElement_;

            CurrentElement_->CoordinateMatrix();
            CurrentElement_->FE()->N( NI_, currXyz );

            Nmin = min_element( NI_.begin(), NI_.end() );
            Nmax = max_element( NI_.begin(), NI_.end() );

            if (debug_)
            {
                if( dim == 2 )
                    std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<"; z="<<currXyz[2]<<std::endl;
                else
                    std::cout<<"Point: x="<<currXyz[0]<<"; y="<<currXyz[1]<<std::endl;

                std::cout<<"Nmax="<<*Nmax<<"; Nmin="<<*Nmin<<std::endl;
                cout<<"\nPropertyAtPointVisitor<dim>::FindPoint: do-cycle: Element id="<<CurrentElement_->Idx()<<", minNi="<<*Nmin<<", maxNI="<<*Nmax<<endl;
            }

            minimumNi = *Nmax;

            // found element, containing given point
            if( ( *Nmax <= (1.0+precision_) ) && (*Nmin >= (0.0-precision_)) )
            {
                stopFlag = true;
                Find     = true;
            }
            else
            {
                stopFlag        = true;
                TargetElement_  = NULL;
                minimumNi       = *Nmax;

                // Visitation of the element neigbours in order to find the element with the min value of shape function
                for( auto i=0; i<CurrentElement_->Neighbors(); i++ )
                {
                    // if not a boundary element
                    if( CurrentElement_->Neighbor(i) != NULL)
                    {
                        checked = false;

                        itChd = listOfCheckedElements.find( CurrentElement_->Neighbor(i)->Idx() );
                        if(itChd != listOfCheckedElements.end() )
                            checked = true;

                        // Element not in the list of checked
                        if( !checked )
                        {

                            current_bc = CurrentElement_->Neighbor(i)->BaryCenter();

                            if( isCloseToBarycenter(currXyz,CurrentElement_->Neighbor(i),current_bc) )
                            {
                                // Shape function in neib
                                CurrentElement_->Neighbor(i)->CoordinateMatrix();
                                CurrentElement_->Neighbor(i)->FE()->N( NI_, currXyz );

                                for(size_t j=0; j<CurrentElement_->Neighbor(i)->Nodes(); j++)
                                    NI_[j] = fabs(NI_[j]);

                                Nmin        = min_element( NI_.begin(), NI_.end() );
                                minimumNi   = std::min(*Nmin,minimumNi);

                                if(debug_)
                                    cout<<" do-cycle: Neigbour Element="<< CurrentElement_->Neighbor(i)->Idx()
                                        <<" minNi="<<*Nmin<<", minimumNi="<<minimumNi<<endl;

                                // if one element stays in the list it stops
                                if( minimumNi == *Nmin )
                                {
                                    TargetElement_ = CurrentElement_->Neighbor(i);
                                    stopFlag = false;

                                    if(debug_)
                                        cout<<" do-cycle: setting NEXT Element="<< CurrentElement_->Neighbor(i)->Idx()<<endl;

                                }

                            }

                            listOfCheckedElements.insert( CurrentElement_->Neighbor(i)->Idx() );

                            if( listOfCheckedElements.size() >= maxElementsInTheMesh_-1 )
                                stopFlag = true;

                        } // !checked

                    }//Neighbor(i)  not a NULL pointer

                }// all neighbours

            } //else - shape f-s not within range - element not found yet

            if(debug_)
            {
                if(TargetElement_ !=NULL)
                    cout<<endl<<" Found Next Neighbor id="<< TargetElement_->Idx()<<endl;
                else
                    cout<<" NextElement not found "<<endl;
            }

            if( ++currIteration_ > maxIterations_ )
            {
                stopFlag     = true;
                TargetElement_ = NULL;
            }

        }// end of do until target found
        while ( !stopFlag );

        //
        // Process special cases: when search is stuck or point out of domain
        //

        // Process is stuck
        if( TargetElement_== NULL )
            //Find = FindPoint_BruteForceSearch( currXyz, e, bc );
            Find = false;

        // warning if the element not found at all
        if( listOfCheckedElements.size() >= maxElementsInTheMesh_-1 )
        {
            csmp_error.notice( ERROR," PropertyAtPointVisitor<>::Visit:Element: Current point not found",
                               " Probably point is out of the meshed domain");
            TargetElement_ = NULL;
        }

        return Find;
    }

    return false;
}






template<uint32_t dim>
void PropertyAtPointVisitor<dim>::Visit( Element<dim>* e)
{
    if( continueSearch_ )
    {

        csmp::ErrorHandler& csmp_error(csmp::ErrorHandler::Instance());

        if(debug_)
            cout<<endl<<" Visit(Element): Initially Visiting Element="<<e->Idx()<<endl;

        // Cycle over all points in the input map of vector (X,Y,Z) (numerous points in the cloud)
        // Cycle over the elements not implemented here
        bool resultOfSearch( false );
        VectorVariable<dim> bc;

        //e->CoordinateMatrix();
        bc = e->BaryCenter();

        for( typename map <size_t,vector<double> >::iterator itX = xyz_.begin(); itX != xyz_.end(); itX++ )
        {
            if(debug_)
                cout<<" Visit(Element): point N="<< itX->first<<" targetFound="<< elements_ids_found_[itX->first]<<endl;

            if( elements_ids_found_[itX->first] == (size_t)-1 )
            {
                if( bruteforce_ )
                    resultOfSearch = FindPoint_BruteForceSearch ( xyz_[itX->first], e, bc );
                 else
                    resultOfSearch = FindPoint_NeighborSearch   ( xyz_[itX->first], e, bc );
            }
            else if(debug_ && itX==xyz_.end() )
                cout<<" PropertyAtPointVisitor<dim>::Visit(Element) Point "<<itX->first<<" not found in domain"<<endl;

            // Point found
            if( resultOfSearch )
            {
                numberOfPointsFound_++;

                elements_found_[itX->first]     = TargetElement_;
                elements_ids_found_[itX->first] = TargetElement_->Idx();

                /////////////////////////////////////////////
                /// Processing of the properties for output
                /////////////////////////////////////////////

                if(debug_)
                    cout<<" For point="<<itX->first<<" found element N="<<TargetElement_->Idx()<<endl;

                if( prop_idx_.type == SCALAR )
                {
                    if (prop_idx_.place == NODE)
                    {
                        TargetElement_->NodePropertyVector( prop_idx_, NPS_ );
                        // Interpolate properties
                        propS_[itX->first]() = 0.0;
                        for ( size_t i=0; i<TargetElement_->Nodes(); i++ )
                            propS_[itX->first]() += NI_[i] * NPS_[i]();

                    }else if (prop_idx_.place == ELEMENT)
                        propS_[itX->first]() = TargetElement_->Read(prop_idx_);

                }
                else if(prop_idx_.type == VECTOR )
                {
                    if (prop_idx_.place == NODE)
                    {
                        TargetElement_->NodePropertyVector( prop_idx_, NPV_ );

                        for(auto i=0;i<dim;i++)
                        {
                            propV_[itX->first](i)=0.0;
                            for ( size_t j=0; j<TargetElement_->Nodes(); j++ )
                                propV_[itX->first](i)+=NI_[j] * NPV_[j](i);
                        }

                    }else if (prop_idx_.place == ELEMENT)
                        TargetElement_->Read(prop_idx_,propV_[itX->first]);

                }
                else if(prop_idx_.type == TENSOR )
                {
                    if (prop_idx_.place == NODE)
                    {
                        TargetElement_->NodePropertyVector( prop_idx_, NPT_ );
                        for(auto i=0;i<dim;i++)
                        {
                            for(auto k=0;k<dim;k++)
                            {
                                propT_[itX->first](i,k)=0.0;
                                for ( size_t j=0; j<TargetElement_->Nodes(); j++ )
                                    propT_[itX->first](i,k)+=NI_[j] * NPT_[j](i,k);
                            }
                        }
                    }else if (prop_idx_.place == ELEMENT)
                        TargetElement_->Read(prop_idx_,propT_[itX->first]);

                }
                else if(prop_idx_.type == ARRAY )
                {
                    if (prop_idx_.place == NODE)
                    {
                        TargetElement_->NodePropertyVector( prop_idx_, NPA_ );
                        for(size_t i=0;i<dim;i++)
                        {
                            propA_[itX->first](i)=0.0;
                            for ( size_t j=0; j<TargetElement_->Nodes(); j++ )
                                propA_[itX->first](i)+=NI_[j] * NPA_[j](i);
                        }
                    }else if (prop_idx_.place == ELEMENT)
                        TargetElement_->Read(prop_idx_,propA_[itX->first]);
                }
                else if(prop_idx_.type == FLAGGEDARRAY )
                {
                    if (prop_idx_.place == NODE)
                    {
                        TargetElement_->NodePropertyVector( prop_idx_, NPFA_ );
                        for(size_t i=0;i<dim;i++)
                        {
                            propA_[itX->first](i)=0.0;
                            for ( size_t j=0; j<TargetElement_->Nodes(); j++ )
                                propFA_[itX->first](i)+=NI_[j] * NPFA_[j](i);
                        }
                    }else if (prop_idx_.place == ELEMENT)
                        TargetElement_->Read(prop_idx_,propFA_[itX->first]);
                }
                else
                {
                    csmp_error.notice( ERROR," PropertyAtPointVisitor<>::Visit:Element() Property type not supported"," Probably wrong type of property");
                }

            }// if  target element found - resulOfSearch, we are checking the properties

            resultOfSearch = false;

        }// cycle over all points in the cloud

        if(debug_)
            cout<<" numberOfPoints = "<<numberOfPoints_<<", numberOfPointsFound = "<<numberOfPointsFound_<<endl;

        if( numberOfPointsFound_ == numberOfPoints_ )
        {
            if(debug_)
                cout<<"  PropertyAtPointVisitor<> All points found..."<<endl;

            SetContinueSearchOff();

        }

    }

}// end Visit(Element<dim>*)





template<uint32_t dim>
void PropertyAtPointVisitor<dim>::Visit(Model<dim>* m)
{

    csmp::Region<dim>&  res(m->Region("Model"));
    for ( auto el_it=res.CellVector().begin(); el_it!=res.CellVector().end(); el_it++ )
       (*el_it)->Accept( *this );
}


template class PropertyAtPointVisitor<1U>;
template class PropertyAtPointVisitor<2U>;
template class PropertyAtPointVisitor<3U>;
}
