#include "ANSYS_Model_Quadratic.h"

using namespace std;

namespace csmp {


/**
@ Description:
Default constructor of Model is called.
The model is built from the 'icem_file_set' '*.asc' and '*.dat', variables file,
and the regions file prefix is used to read the regions file.
*/

  template<size_t dim>
  ANSYS_Model_Quadratic<dim>
  ::ANSYS_Model_Quadratic( const std::string& mesh_file_set,
                       const std::string& regions_file_prefix,
                       const std::string& variable_file )
   : Model<dim>( variable_file.c_str(), false ),
     meshFileName_      ( mesh_file_set),
     topologyFileName_  ( regions_file_prefix ),
     interpolationOrder_( LINEAR )
   {
      Initialize( );
   }



  template<size_t dim>
  ANSYS_Model_Quadratic<dim>::~ANSYS_Model_Quadratic()
    {
      EraseElements();
      EraseNodes();
    }


  /// Builds Model after it was constructed with the default constructor.

  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::Initialize()
    {
      // IO
      AnnounceStart();
      // read in regions and boundaries to consider
      ReadInSubDomains();
      // IO
      AnnounceSubDomains();
      // asc file
      ReadInFamilies();
      // dat file
      ReadInData();
      // initialize line element neighbor connectivity (icem does not do that)
      NeighborConnectivity();
      // delete low dim elements in higher dim regions
      DeleteRedundantElements();
      // form regions and initialize containers
      CreateRegions();
      // boundaries from regions
      CreateBoundaries();
      // split boundaries from boundaries
      CreateSplitBoundaries();
      // IO
      AnnounceEnd();
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::ReadInSubDomains()
    {
      string fileName( (string)(topologyFileName_ + "-subdomains.txt") );
      cout << endl << "Reading subdomain file " << fileName << "..." << endl;
      ifstream file( fileName.data(), ios::in );

      if( !file.is_open() )
        {
          string errorMessage( (string)("Unable to open subdomain input file " + fileName ) );
          throw csmp::Exception( ERROR, "ANSYS_Model_Quadratic<dim>::ReadInSubDomains()", errorMessage );
        }        

      string line; 

      regions_.clear(); 
      boundaries_.clear();
      splitBoundaries_.clear();
      subdomains_.clear();

      // skipping to regions
      while( IsNot( file, line, "regions" ) );
      // reading regions
      while( IsNot( file, line, "boundaries" ) ){regions_.push_back(line);}
      // reading boundaries
      while( IsNot( file, line, "split" ) ){boundaries_.push_back(line);}
      // reading split boundaries
      while( IsNot( file, line, "eof" ) ){splitBoundaries_.push_back(line);}

      subdomains_.insert( regions_.begin(), regions_.end() );
      subdomains_.insert( boundaries_.begin(), boundaries_.end() );
      subdomains_.insert( splitBoundaries_.begin(), splitBoundaries_.end() );
      regions_.clear();
      regions_.insert( regions_.begin(), subdomains_.begin(), subdomains_.end() );
      boundaries_.insert( boundaries_.begin(), splitBoundaries_.begin(), splitBoundaries_.end() );

      file.close();
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::ReadInFamilies()
    {
      string fileName( (string)(meshFileName_ + ".asc") );
      ifstream file( fileName.data(), ios::in );

      if( !file.is_open() )
        throw csmp::Exception( ERROR, "ANSYS_Model_Quadratic<dim>::ReadInRegions()", "Unable to open ascii input file!" );

      string line; 

      // region overview
      while( IsNot( file, line, "families" ) );
      const int regionCount( ExtractEntityCount(line) );

      cout << endl << "Reading " << fileName << " Icem model topology with " << regionCount << " families..." << endl;

      regionElements_.clear();
      while( IsNot( file, line, "Objectname" ) );
      for( size_t r(0); r < regionCount; ++r )
        ExtractRegionInfo( file, line );

      CheckSubDomainRequest();

      elementsInUse_.clear();
      for( size_t r(0); r < regionElements_.size(); ++r )
        {
          if( subdomains_.find( regionElements_[r].first ) == subdomains_.end() )
            {
              cout << endl << "  Family " <<  regionElements_[r].first << " discarded..." << endl;
              continue;
            }
          while( IsNot( file, line, regionElements_[r].first ) );
          for( size_t e(0); e != regionElements_[r].second.size(); ++e )
            {
              file >> regionElements_[r].second[e];
              elementsInUse_.insert( regionElements_[r].second[e] );
            }
        }

      file.close();
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::ExtractRegionInfo( ifstream& file, string& line )
    {      
      getline( file, line );
      size_t runningPos( line.find_first_of( " " ) );
      string regionName( line.substr( 0, runningPos ) );
      families_.insert(regionName);
      if( RegionInfoExists(regionName) )
        return;
      runningPos = line.find_last_of( " " );
      string regionElementNumber( line.substr( runningPos) );
      size_t regionElementCount = std::stol(regionElementNumber);
      cout << endl << "  Family " << regionName << " found with " << regionElementCount << " eligible element(s)." << endl;
      regionElements_.push_back( make_pair( regionName, vector<size_t>(regionElementCount) ) );
    }




  template<size_t dim>
  bool ANSYS_Model_Quadratic<dim>::ReadInData()
    {
      string fileName( (string)(meshFileName_ + ".dat") );
      ifstream file( fileName.data(), ios::in );

      if( !file.is_open() )
        throw csmp::Exception( ERROR, "ANSYS_Model_Quadratic<dim>::ReadInData()", "Unable to open data input file!" );

      cout << endl << "Reading " << fileName << " Icem model data file..." << endl;

      string line;     

      // PX PY PZ
      while( IsNot( file, line, "PX" ) );
      const int nodeCount( ExtractEntityCount( line ) );

      cout << endl << "  Reading px, pz, py with a total of " << nodeCount << " coordinates..." << endl;

      vector<Point<dim> > points(nodeCount);
      for( size_t d(0); d < dim; ++d )
        for( size_t i(0); i < nodeCount; ++i )
          {
            double c(0.);
            file >> c;
            points[i][d] = c;
          }

      cout << endl << "  Creating a total of " << nodeCount << " csmp::Nodes..." << endl;

      // we create the nodes here
      EraseNodes();
      nodes_.resize( nodeCount, NULL );
      const LocalVariables nvars( this->Database().LocalVariablesAt(NODE) );
      for( size_t n(0); n < nodeCount; ++n )
        {
          Node<dim>* newNode = new Node<dim>;
          newNode->Idx(n);
          newNode->Coordinate( points[n] );
          newNode->ResizePropertyStorage(nvars);
          nodes_[n] = newNode;
        }

      // SKIP UNTIL PELEMENTS
      while( IsNot( file, line, "PELEMENTS" ) );

      // PELEMENTS
      const int elementCount( ExtractEntityCount( line ) );
      
      cout << endl << "  Reading types of " << elementCount  << " Elements..." << endl;

      int32 ansysFemType(999);
      vector<CSMP_FEM_TYPE> femTypes;
      typedef ANSYS_ElementSpecifications ansysSpecs;
      for( size_t e(0); e < elementCount; ++e )
        {
          file >> ansysFemType;
          femTypes.push_back( ansysSpecs::CSMP_TypeFrom_ANSYS_Type( ansysFemType, true, dim ) );
        }
      // we assume either all linear, all quadratic or all cubic here
      interpolationOrder_ = ansysSpecs::InterpolationOrder(ansysFemType);
      cout << endl << "  Elements are ";
      if( interpolationOrder_ == LINEAR ) cout << "linear." << endl;
      else if( interpolationOrder_ == QUADRATIC ) cout << "quadratic." << endl;
      else cout << "cubic." << endl;

      // SKIP UNTIL PLIST
      while( IsNot( file, line, "PLIST" ) );

      // PLIST
      const int plistCount( ExtractEntityCount( line ) );

      cout << endl << "  Reading " << plistCount  << " element node entries and creating " << elementCount<< " csmp::Elements..." << endl;

      // we track how many parents there are per element
      multiset<size_t> parentCount;

      // we create the Elements here
      this->FE_Manager().InitializeElements( dim, interpolationOrder_, true );
      EraseElements();
      const LocalVariables evs( this->Database().LocalVariablesAt(ELEMENT) );
      const IntegrationPointVariables ivs( this->Database().IntegrationPointVariablesAt(ELEMENT) );
      elements_.resize( elementCount, NULL );
      for( size_t e(0); e < elementCount; ++e )
        {
          Element<dim>* newElement = new Element<dim>( this->FE_Manager().E( femTypes[e] ), NULL, evs, ivs );
          newElement->Idx(e);
          for( size_t en(0); en < newElement->Nodes(); ++en )
            {
              size_t enIdx(0);
              file >> enIdx;
              parentCount.insert(enIdx);
              newElement->Assign( en, nodes_.at(enIdx) );              
            }
          for( size_t enb(0); enb < newElement->Neighbors(); ++enb )
            newElement->Assign( enb, static_cast<Element<dim>*>(NULL) );
          elements_[e] = newElement;
        }

      // SKIP UNTIL PFVERTS
      while( IsNot( file, line, "PFVERTS" ) );

      // PFVERTS
      const int pfvertsCount( ExtractEntityCount( line ) );

      cout << endl << "  Reading " << pfvertsCount  << " element neighbor connection entries..." << endl;

      size_t pfverts(0);
      for( size_t e(0); e < elementCount; ++e )
        {
          Element<dim>* const currentElement( elements_[e] );
          for( size_t enb(0); enb < currentElement->Neighbors(); ++enb )
            {
              size_t enbIdx(0);
              file >> enbIdx;
              currentElement->Assign( enb, elements_.at(enbIdx) );
              ++pfverts;
            }
          if( pfverts == pfvertsCount )
            break;
        }

      cout << endl << "  Done reading raw Icem mesh data." << endl;

      file.close();

      cout << endl << "  Connecting nodes to parent elements..." << endl;

      // Initialize node parents
      for( size_t n(0); n < nodeCount; ++n )
        nodes_[n]->ResizeParentStorage( parentCount.count(n) );

      // Assign node parents
      for( size_t e(0); e < elementCount; ++e )
        {
          Element<dim>* const currentElement( elements_[e] );
          for( size_t en(0); en < currentElement->Nodes(); ++en )
            currentElement->N(en)->Assign( en, currentElement );
        }

      return true;
    }


  template<size_t dim>
  bool ANSYS_Model_Quadratic<dim>::IsNot( std::ifstream& file, std::string& line, const std::string& isNot ) const
    {
      getline( file, line );
      if( line.find( isNot ) == string::npos )
        return true;
      return false;
    }


  template<size_t dim>
  int ANSYS_Model_Quadratic<dim>::ExtractEntityCount( const std::string& line )
    {
      const size_t endOfCount( line.find_first_of(" ") );
      const string countString( line.substr(0,endOfCount) );
      return std::stoi(countString);
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::EraseNodes()
    {
      for( size_t n(0); n < nodes_.size(); ++n )
          if( nodes_[n] )
            delete nodes_[n];
      nodes_.clear();
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::EraseElements()
    {
      for( size_t n(0); n < elements_.size(); ++n )
        if( elements_[n] )
          delete elements_[n];
      elements_.clear();
    }


  /// For sake of speed, this assumes line element only region and does not test for that.
  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::InitializeLineElementNeighbors( size_t regionIndex )
    {
      vector<size_t>& regionElements( regionElements_[regionIndex].second );
      const size_t regionElementCount( regionElements.size() );
      for( size_t e(0); e < regionElementCount; ++e )
        {
          Element<dim>* const currentElement( elements_.at( regionElements[e] ) );
          const size_t node1idx( currentElement->N(0)->Idx() );
          const size_t node2idx( currentElement->N(1)->Idx() );
          for( size_t en(0); en < regionElementCount; ++en )
            {
              if( en == e )
                continue;
              Element<dim>* possibleNeighbor( elements_.at( regionElements[en] ) );
              if( possibleNeighbor->N(0)->Idx() == node1idx || possibleNeighbor->N(1)->Idx() == node1idx )
                currentElement->Assign( 0, possibleNeighbor );
              if( possibleNeighbor->N(0)->Idx() == node2idx || possibleNeighbor->N(1)->Idx() == node2idx )
                currentElement->Assign( 1, possibleNeighbor );
            }
        }     
    }


  template<size_t dim>
  bool ANSYS_Model_Quadratic<dim>::RegionInfoExists( const string& regionName )
    {
      for( size_t r(0); r < regionElements_.size(); ++r )
        if( regionElements_[r].first == regionName )
          return true;
      return false;
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::NeighborConnectivity()
    {
      cout << endl << "Connecting line elements..." << endl;
      for( size_t r(0); r < regionElements_.size(); ++r )
        if( elements_.at( regionElements_[r].second.at(0) )->IsLineElement() )
          InitializeLineElementNeighbors(r);
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::CreateRegions()
    {
      cout << endl << "Forming csmp::Regions..." << endl;
      // a. 'Model'
      CreateRegion( "Model", elements_.begin(), elements_.end() );

      // b. all .asc regions which are in subdomains file
      const size_t elementCount( elements_.size() );
      for( size_t r(0); r < regionElements_.size(); ++r )
        {
          if( subdomains_.find( regionElements_[r].first ) == subdomains_.end() )
            continue;
          const size_t regionElementCount( regionElements_[r].second.size() );
          vector<Element<dim>*> regionElements( regionElementCount, NULL );
          for( size_t re(0); re != regionElementCount; ++re )
            {
              for( size_t e(0); e != elementCount; ++e )
                if( elements_[e]->Idx() == regionElements_[r].second[re] )
                  regionElements[re] = elements_[e];
            }
          this->CreateRegion( regionElements_[r].first, regionElements.begin(), regionElements.end() );
        }
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::CreateBoundaries()
    {
      cout << endl << "Forming csmp::Boundaries...";
      if( boundaries_.empty() )
        cout << "...none specified." << endl;
      else
        cout << endl;

      for( size_t b(0); b < boundaries_.size(); ++b )
        this->InsertBoundary( boundaries_[b].data() );
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::CreateSplitBoundaries()
    {
      cout << endl << "Forming csmp::SplitBoundaries...";
      if( splitBoundaries_.empty() )
        cout << "...none specified." << endl;
      else
        cout << endl;

      for( size_t sb(0); sb < splitBoundaries_.size(); ++sb )
          this->InsertSplitBoundary( splitBoundaries_[sb], false/* do not delete region */ );
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::CreateRegion( const std::string& regionName, typename std::vector<Element<dim>*>::const_iterator elementsBegin,
                                       typename std::vector<Element<dim>*>::const_iterator elementsEnd )
    {
      pair<typename map<string,csmp::Region<dim> >::iterator,bool>
        it = this->uniqueGroupMap_.insert( make_pair( regionName, csmp::Region<dim>( regionName, this->Database() ) ) );
      if ( it.second )
        {
          (*it.first).second.SimplexVector().assign( elementsBegin, elementsEnd );
          (*it.first).second.CreateNodePointerVector();
          (*it.first).second.IdentifyPerimeter( );
        }
      cout << endl << "  csmp::Region " << regionName << " created." << endl;
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::DeleteRedundantElements()
    {
      size_t redundantElements(0);
      cout << endl << "Deleting lower dimensional elements from higher dimensional regions...";
      const set<size_t>::const_iterator elementsInUseEnd( elementsInUse_.end() );
      size_t e(0);
      while( e < elements_.size() )
        {
          if( elementsInUse_.find( elements_[e]->Idx() ) == elementsInUseEnd )
            {
              // Detach Elements from its Nodes and Neighbors
              elements_[e]->UnassignNodes();
              elements_[e]->UnassignNeighbors();

              delete elements_[e];
              elements_.erase( elements_.begin()+e );

              ++redundantElements;
            }
          else
            ++e;
        }
      cout << "..." << redundantElements << " element(s) removed." << endl;
    }



  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::AnnounceSubDomains() const
    {
      AnnounceSubDomains( regions_, "Regions" );
      AnnounceSubDomains( boundaries_, "Boundaries" );
      AnnounceSubDomains( splitBoundaries_, "SplitBoundaries" );
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::AnnounceSubDomains( const std::vector<std::string>& subdomains, const std::string& subdomain ) const
    {
      cout << endl << "  " << subdomain << " to form:";
      if( subdomains.empty() )
        { cout << endl << "    NONE" << endl; return; }
      for( size_t i(0); i < subdomains.size(); ++i )
        cout << endl << "    " << subdomains[i];
      cout << endl;
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::AnnounceStart() const
    {
      cout << endl << "===================================" << endl;
      cout << "ANSYS_Model<" << dim << "> for Icem CSP FE mesh\n\n";
      cout << "Model file: " << meshFileName_<< endl;
      cout << "SubDomain file: " << topologyFileName_ << endl;
      cout << "===================================" << endl;
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::AnnounceEnd() const
    {
      cout << endl << "===================================" << endl;
      cout << "ANSYS_Model<" << dim << "> for Icem CSP FE mesh\n\n";
      cout <<  meshFileName_ << " created!" << endl;
      cout <<  "Elements: " << elements_.size() << endl;
      cout <<  "Nodes: " << nodes_.size() << endl;
      cout <<  "Mesh type: " << ParseInterpolationOrder() << endl;
      cout <<  "Regions: " << this->Regions()  << endl;
      cout <<  "Boundaries: " << this->Boundaries()  << endl;
      cout <<  "Split Boundaries: " << this->SplitBoundaries()  << endl;
      cout << "===================================" << endl;
    }


  template<size_t dim>
  std::string ANSYS_Model_Quadratic<dim>::ParseInterpolationOrder() const
    {
      if( interpolationOrder_ == LINEAR ) return "linear";
      if( interpolationOrder_ == QUADRATIC ) return "quadratic";
      if( interpolationOrder_ == CUBIC ) return "cubic";

	  string errorMessage( "Undefined order of approximation!"  );
      throw csmp::Exception( ERROR, "ANSYS_Model_Quadratic<dim>::ParseInterpolationOrder", errorMessage );

	  return errorMessage;
    }


  template<size_t dim>
  void ANSYS_Model_Quadratic<dim>::CheckSubDomainRequest() const
    {
      for( set<string>::const_iterator it( subdomains_.begin() ); it != subdomains_.end(); ++it )
        if( families_.find(*it) == families_.end() ) 
          {
            string errorMessage(  string( "SubDomain " + *it + " specified in -subdomain.txt file not available in mesh!" ) );
            throw csmp::Exception( ERROR, "ANSYS_Model_Quadratic<dim>::CheckSubDomainRequest", errorMessage );
          }
    }


  template class ANSYS_Model_Quadratic<2>;
  template class ANSYS_Model_Quadratic<3>;

} // end csmp
