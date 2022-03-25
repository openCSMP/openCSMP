#include "VSet.h"
#include "ANSYS_ElementSpecifications.h"
#include "Exception.h"
#include "binaryReadWrite.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
VSet<dim>::VSet()
{
}

/**

Allocates storage for mono-element mesh from finite element blue-print
and the number of nodes and elements.

@section implementation Implementation

Element types gets a size of 1 in this case which is just enough to
store the type of the single element of which the mesh consists.
*/
template<uint32_t dim>
VSet<dim>::VSet( uint32_t nodes_per_element,
                 uint32_t nbors_per_element,
                 int8_t etype,
                 size_t nodes, size_t elmts )
: VData(nodes_per_element, nbors_per_element, nodes, elmts),
  pmtrl_( elmts, UNSPECIFIED )
{
	SingleElementType(etype);
}


/// creates empty VSet of the desired dimensions
template<uint32_t dim>
VSet<dim>::VSet(const deque<uint32_t>& npes,
                const deque<uint32_t>& epes,
                size_t nodes)
: VData(npes, epes, nodes),
  pmtrl_( epes.size(), UNSPECIFIED )
{
	SingleElementType(0);
}



/// copy constructor
template<uint32_t dim>
VSet<dim>::VSet( const VSet<dim>& a )
 : VData(a), pmtrl_(a.pmtrl_), property_map_(a.property_map_)
 {
 }


/** Assignment operator
*/
template<uint32_t dim>
VSet<dim>& VSet<dim>::operator=( const VSet<dim>& a )
  {
    if (&a != this) {
         VData::operator=( a ); 
         pmtrl_        = a.pmtrl_; 
         property_map_ = a.property_map_;
      }
    return *this;
  }


/**

Resizes the VSet internal containers to hold a mono-element type mesh
of the specified size.
*/
template<uint32_t dim>
void VSet<dim>::Resize( uint32_t nodes_per_element,
                        uint32_t nbors_per_element,
                        int8_t csmp_etype,
                        size_t nodes,
                        size_t elmts )
{
  // ascertaining that the information matches the other entries in the VSet
  if ( HybridElementTypeMesh() )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "this method is only for single-element type meshes.");

  if ( Elements() > 0 && ElementType(0) != csmp_etype )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "element type does not match the ones stored in VData.");

	VData::Resize(nodes_per_element, nbors_per_element, csmp_etype, nodes, elmts);
  pmtrl_.resize( elmts );
  
  if ( !property_map_.empty() )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "resizing of property map not handled yet.");
}


/**

Just resizes all storage according to specifications without actually
assigning any values to the subdeques.
*/
template<uint32_t dim>
void VSet<dim>::Resize( const deque<int8_t>& etypes,
                        const deque<uint32_t>& npes,
                        const deque<uint32_t>& epes,
                        size_t nodes, size_t faces, size_t interfaces )
{
  if ( etypes.size() == 1 )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "use this method only for hybride meshes with multiple element types" );

  if ( etypes.size() != npes.size() || npes.size() != epes.size() )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "supplied deques must all have the same size equal to the number of finite elements that VSet shall have");

  // resizing the polygonal data
	VData::Resize(etypes, npes, epes, nodes, faces, interfaces);

  if ( etypes.size() != TotalNumberOfCells() )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "element types does not match number of finite elements in 'plist'");

  pmtrl_.resize( epes.size() - faces - interfaces );
  
  if ( !property_map_.empty() )
    throw csmp::Exception( ERROR, "VSet<dim>::Resize", "resizing of property map not handled yet.");
}




template<uint32_t dim>
VSet<dim>::~VSet()
{
}




/**
Adds node coordinates to VSet
*/
template<uint32_t dim>
void VSet<dim>::AddXYZ( const deque<double>& x,
                        const deque<double>& y,
                        const deque<double>& z )
{
	if (x.size() != Vertices())
    {
      if (Vertices() > 0) {
          cout << "\nVSet<"<< dim <<">::AddXYZ: Warning: Changing node-coordinate " << endl;
          cout << "array size from " << Vertices() << " to " << x.size() << endl;
        }
      VData::ResizeNodes(x.size());
    }

	for (auto i = 0; i<x.size(); i++) Px(i, x[i]);

	if constexpr ( dim > 1 )
    {
      assert(x.size() == y.size());
      for (auto i = 0; i<y.size(); i++) Py(i, y[i]);
    }

	if constexpr ( dim > 2 )
    {
      assert(y.size() == z.size());
      for (auto i = 0; i<x.size(); i++) Pz(i, z[i]);
    }
    
   ResizeBFlags();
}



/**
Assigns record of nodes per element to VSet.
Since the entries are serialised, the record 'pelmt' is needed in order
to recuperate element types and the nodes per element information when
reading the VSet.
*/
template<uint32_t dim>
void VSet<dim>::AddPlist( typename map<size_t, vector<int64_t> >::const_iterator first,
						              typename map<size_t, vector<int64_t> >::const_iterator last)
{
	typename deque<vector<int64_t> >::iterator it = PlistBegin();

	while (first != last && it != PlistEnd())
    {
      (*it) = (*first).second;
      first++;
      it++;
    }
    
} // end


	/**
	Assigns record of nodes per element to VSet.
	Since the entries are serialised, the record 'pelmt' is needed in order
	to recuperate element types and the nodes per element information when
	reading the VSet.
	*/
template<uint32_t dim>
void VSet<dim>::AddPlist( typename deque<vector<int64_t> >::const_iterator first,
                          typename deque<vector<int64_t> >::const_iterator last )
{
	typename deque<vector<int64_t> >::iterator it = PlistBegin();

	while (first != last && it != PlistEnd())
    {
      (*it) = (*first);
      first++;
      it++;
    }
    
} // end


	/**
	Adds neighborhood information: volumetric element neighbors for volume elements,
	surface element neighbors for surface elements, and line element neighbors for
	line elements.

	@attention lower-dimensional elements may have multiple neighbors per face.
	These manifolds cannot be captured by this classical neighbor record.
	*/
template<uint32_t dim>
void VSet<dim>::AddPfverts( typename map<size_t, vector<int64_t> >::const_iterator first,
						                typename map<size_t, vector<int64_t> >::const_iterator last )
{
	typename deque<vector<int64_t> >::iterator it = PfvertsBegin();

	while (first != last && it != PfvertsEnd())
    {
      (*it) = (*first).second;
      first++;
      it++;
    }
    
} // end



	/**
	Adds neighborhood information: volumetric element neighbors for volume elements,
	surface element neighbors for surface elements, and line element neighbors for
	line elements.

	@attention lower-dimensional elements may have multiple neighbors per face.
	These manifolds cannot be captured by this classical neighbor record.
	*/
template<uint32_t dim>
void VSet<dim>::AddPfverts( typename deque<vector<int64_t> >::const_iterator first,
						                typename deque<vector<int64_t> >::const_iterator last )
{
	typename deque<vector<int64_t> >::iterator it = PfvertsBegin();

	while (first != last && it != PfvertsEnd())
    {
      (*it) = (*first);
      first++;
      it++;
    }
    
} // end


	/**
	Adds a map that stores the indices of the boundary nodes as keys and
	BOX_BOUNDARY flags as values.
	*/
template<uint32_t dim>
void VSet<dim>::AddBFlags( typename vector<std::int8_t>::const_iterator first,
                           typename vector<std::int8_t>::const_iterator last )
{
   const size_t bflags_size(distance(first,last));
   if ( bflags_size != Vertices() ) {
        cerr <<"\nVSet<dim>::AddBFlags: supplied BOX_BOUNDARY flag range does not match the number of nodes: ";
        cerr << bflags_size <<" vs. "<< Vertices() << endl;
        cerr <<"No assignments were made.\n";
        return;
     }
     
   ResizeBFlags();
   auto bit=BFlagsBegin();
	 while ( first != last ) {
       (*bit) = (*first);
       first++;
       bit++;
     }
    
} // end




/**
       Material ID identifiers need to be provided for all elements, boundaries and split boundaries.
*/
template<uint32_t dim>
void VSet<dim>::AddPmtrl( typename std::vector<int32_t>::const_iterator first,
                          typename std::vector<int32_t>::const_iterator last )
 {
    pmtrl_.clear();
    pmtrl_.assign( first, last );
    if ( pmtrl_.size() != Elements() )
      cerr <<"\nVSet<"<< dim <<">::AddPmtrl: mismatch between the number of elements and material IDs in vset.\n";
 }


template<uint32_t dim>
std::vector<int32_t>::const_iterator VSet<dim>::PmtrlBegin() const
 { return pmtrl_.begin(); }

template<uint32_t dim>
std::vector<int32_t>::const_iterator VSet<dim>::PmtrlEnd() const
 { return pmtrl_.end(); }




/**
Reports whether the VSet contains distributed variable values
stored in PropertyData object.s
*/
template<uint32_t dim>
bool csmp::VSet<dim>::DataEmpty() const
{
	return property_map_.empty();
}




/**
Adds variable dataset to 'property_map_' data member of VSet
New much more memory efficient variable storage framework.

@author SKM
@data 4/4/2016

@note NEW!
*/
template<uint32_t dim>
bool VSet<dim>::AddData( const char* s, const PropertyData& data )
{
	auto result = property_map_.insert(make_pair(s, data));
	return result.second;

} // end AddData





/**
   To delete a dataset
*/
template<uint32_t dim>
void VSet<dim>::RemoveData( const char* s )
{
   property_map_.erase(s);

} // end RemoveData





/**
Retrieves property data from VSet (if any)

@author SKM
@data 4/4/2016

@note NEW!
*/
template<uint32_t dim>
PropertyData  VSet<dim>::Data(const char* s) const
{
	auto prop_it = property_map_.find(s);
	if (prop_it == property_map_.end())
		throw csmp::Exception(ERROR, "VSet<dim>::Data:", s, "property data was not found in VSet.");

	return (*prop_it).second;

} // end Data


/**
Iterator to map of PropertyData records
*/
template<uint32_t dim>
std::map<std::string, PropertyData>::const_iterator VSet<dim>::PropertyValuesBegin() const {
	return property_map_.begin();
}


/**
Iterator to map of PropertyData records
*/
template<uint32_t dim>
std::map<std::string, PropertyData>::const_iterator VSet<dim>::PropertyValuesEnd() const {
	return property_map_.end();
}


/// if the property records contain FV data
template<uint32_t dim>
bool  VSet<dim>::ContainsFiniteVolumeIntegrationPointData() const
{
	for (auto it = property_map_.begin(); it != property_map_.end(); ++it)
		if ( (*it).second.Placement() == SECTOR_INTEGRATION_POINT or
			   (*it).second.Placement() == FACET_INTEGRATION_POINT or
			   (*it).second.Placement() == FACE_SECTOR_INTEGRATION_POINT or
			   (*it).second.Placement() == FACE_FACET_INTEGRATION_POINT or
			   (*it).second.Placement() == INTER_FACE_SECTOR_INTEGRATION_POINT or
			   (*it).second.Placement() == INTER_FACE_FACET_INTEGRATION_POINT) return true;
	return false;
}



/**

Writes the content of the VSet to a binary data file (including property
data. This output is ordered in the following way:

1. File header: char* s.
2. Mesh connectivity: px, py, pz, plist, pfverts, bflags, bvals.
3. Property data: scalars, vectors, tensors.

@section implementation Implementation

The binary writing is done with the templatized set of functions declared
in 'binaryReadWrite.h'. These can read and write all CSMP type of datasets.
*/
template<uint32_t dim>
bool  VSet<dim>::OutputTo( const char* bin_file, double time ) const
{
	char file_name[200], num[20];
	sprintf(num, "%lf", time);
	strcpy(file_name, bin_file);
	size_t  records(0);

	// 1. opening the file
	fstream fp(file_name, ios::out | ios::binary);
	if (!fp.is_open()) {
		cout << "\nVSet<dim>::OutputTo: File: " << file_name;
		cout << " could not be opened" << endl;
		return false;
	}

	// 2. writing the file header
	{
		BinaryFileSectionWrite sect(fp, "VSETHEDR");
		char heading[200];
		strcpy(heading, "VSet<dim>::OutputTo: Binary version of VSet: ");
		strcat(heading, bin_file);
		strcat(heading, " saved at time: ");
		strcat(heading, num);
		binaryFileWrite(fp, heading);
	}

	// 3. Writing the mesh connectivity to file
	{
		BinaryFileSectionWrite sect(fp, "VSETCONN");
		OutBinary(fp);
	}

	// 4. Writing material / rocktype identifiers ('pmtrl' keys from the elements)
	{
		BinaryFileSectionWrite sect(fp, "VSETMTRL");
    // number of property records
    records = pmtrl_.size();
    fp.write( reinterpret_cast<const char*>(&records), sizeof(int32_t));
    // individual records (all together)
    fp.write( reinterpret_cast<const char*>(&pmtrl_[0]), sizeof(int32_t) * records );
	}

	// 5. Writing the property data records to file
	{
		BinaryFileSectionWrite sect(fp, "VSETPROP");
		if (!property_map_.empty()) {
        // number of property records
        records = property_map_.size();
        fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t));
        // individual records
        for (auto it = property_map_.begin(); it != property_map_.end(); ++it) {
          // writing the property name
          binaryFileWrite(fp, (*it).first.c_str());
          // writing the dataset
          (*it).second.OutBinary(fp);
        }
     }
		else { // no data record is registered for later reading
			records = 0U;
			fp.write( reinterpret_cast<const char*>(&records), sizeof(size_t));
		}
	}


	// 6. cleaning up
  BinaryFileSectionWrite sect(fp, "VSETFOTR");

	fp.close();
	cout << "\nVSet<" << dim << ">::OutputTo: VSet has been successfully written to: ";
	cout << file_name << endl;

	return true;

} // end OutputTo



/**
Key method for recovery of a model from binary file. 
*/
template<uint32_t dim>
bool  VSet<dim>::InputFrom( const char* bin_file, double& time )
{
   const set<string> empty_subset;
   return InputFrom( bin_file, time, empty_subset );
}


/**
Key method for recovery of a model from binary file. It can load only a subset of variables if neccesary.
*/
template<uint32_t dim>
bool  VSet<dim>::InputFrom( const char* bin_file, double& time, const set<string>& subset_variables )
{
	char file_name[NAME_STRING];
	strcpy(file_name, bin_file);
	size_t  records(0);
	string  dname;

	// 1. opening the file
	fstream fp(file_name, ios::in | ios::binary);
	if (!fp.is_open()) {
		cout << "\nVSet<dim>::InputFrom: File: " << file_name;
		cout << " could not be opened" << endl;
		return false;
	}
	// 2. reading the file header and extracting time
	{
		BinaryFileSectionRead sect(fp, "VSETHEDR");
		char heading[INFO_STRING];

		binaryFileRead(fp, heading);
		cout << "\nVSet<dim>::InputFrom: Reading: " << heading << endl;
		strtok(heading, ":");
		strtok(NULL, ":");
		strtok(NULL, ":");
		strtok(NULL, ":");
		time = atof(strtok(NULL, ":"));
	}

	// 3. reading the mesh connectivity to file (VData)
	{
		BinaryFileSectionRead sect(fp, "VSETCONN");
		cout << "\nVSet<dim>::InputFrom: reading finite element mesh..." << endl;
		InBinary(fp);
    // consistency check
    if ( SpatialDimension() != dim ) {
        throw csmp::Exception( ERROR, "VSet<dim>::InputFrom",
                              "mismatch between spatial dimension of VData (mesh) and VSet<(spatial) dim>.");
      }
	}

	// 4. reading material / rocktype identifiers ('pmtrl' keys from the elements)
	{
		BinaryFileSectionRead sect(fp, "VSETMTRL");
    // number of property records
    fp.read( reinterpret_cast<char*>(&records), sizeof(int32_t));
    pmtrl_.resize( records );
    // individual records (all together)
    fp.read( reinterpret_cast<char*>(&pmtrl_[0]), records * sizeof(int32_t) );
	}

	// Reading the property data records from file (PropertyData)
	{
		BinaryFileSectionRead sect(fp, "VSETPROP");

		fp.read( reinterpret_cast<char*>(&records), sizeof(size_t));
		if (records > 0)
			// reading the datasets sequentially
			for (auto i = 0; i<records; ++i)
			{
				// reading the property name
				char heading[INFO_STRING];
				binaryFileRead(fp, heading);
				dname = heading;

        auto prop = make_pair( dname, inBinaryPropertyData( fp ) );
        if ( !subset_variables.empty() ) {
            if ( subset_variables.find( dname ) != subset_variables.end() )
              property_map_.insert( prop );
          }
        else  property_map_.insert( prop );
			}
		else cout << "\nVSet<dim>::InputFrom: no PropertyData objects detected." << endl;
	}

	// 5. cleaning up
  BinaryFileSectionRead sect(fp, "VSETFOTR");

	fp.close();

	cout << "\nVSet<" << dim << ">::InputFrom: VSet has been successfully read from: '";
	cout << file_name << "'." << endl;

	return true;

} // end InputFrom


/**

Writes the content of the VSet to a binary data file (including property
data. This output is ordered in the following way:

1. File header: char* s.
2. Mesh connectivity: px, py, pz, plist, pfverts, bflags, bvals.
3. Property data: scalars, vectors, tensors.

@section implementation Implementation

The binary writing is done with the templatized set of functions declared
in 'binaryReadWrite.h'. These can read and write all CSMP type of datasets.

TODO: @todo Refactor to work with PropertyData based variable storage

*/
template<uint32_t dim>
bool  VSet<dim>::ParallelOutputTo(const char* bin_file, double time, size_t first_outerhalo) const
{
	char file_name[200], num[20];
	sprintf(num, "%lf", time);
	strcpy(file_name, bin_file);
	char heading[200];
	strcpy(heading, "VSet<dim>::OutputTo: Binary version of VSet: ");
	strcat(heading, bin_file);
	strcat(heading, " saved at time: ");
	strcat(heading, num);
	sprintf(num, "%lu", first_outerhalo );
	strcat(heading, ", first outerhalo: ");
	strcat(heading, num);

	// 1. opening the file
	fstream fp(file_name, ios::out | ios::binary);
	if (!fp.is_open()) {
		cout << "\nVSet<dim>::ParallelOutputTo: File: " << file_name;
		cout << " could not be opened" << endl;
		return false;
	}
	// 2. writing the file header   
	binaryFileWrite(fp, heading);

	// 3. Writing the mesh connectivity to file
	OutBinary(fp);

	cerr << "\nVSet<dim>::ParallelOutputTo: variable output has not been implemented yet.\n";

	// 5. cleaning up
	fp.close();
	cout << "\nVSet<" << dim << ">::ParallelOutputTo: VSet has been successfully written to: ";
	cout << file_name << endl;

	return true;
} // end ParallelOutputTo


/**
TODO: @todo Refactor to work with PropertyData based variable storage
*/
template<uint32_t dim>
bool  VSet<dim>::ParallelInputFrom(const char* bin_file, double& time, size_t& first_outerhalo)
{
	cerr << "\nVSet<dim>::ParallelInputFrom: variable output has not been implemented yet.\n";

	char file_name[200], heading[300];
	strcpy(file_name, bin_file);
	string dname;

	// 1. opening the file
	fstream fp(file_name, ios::in | ios::binary);
	if (!fp.is_open()) {
		cout << "\nVSet<dim>::ParallelInputFrom: File: " << file_name;
		cout << " could not be opened" << endl;
		return false;
	}
	// 2. reading the file header and extracting time
	binaryFileRead(fp, heading);
	cout << "\nVSet<dim>::ParallelInputFrom: Reading: " << heading << endl;
	strtok(heading, ":");
	strtok(NULL, ":");
	strtok(NULL, ":");
	strtok(NULL, ":");
	time = atof(strtok(NULL, ":"));
	first_outerhalo = atoi(strtok(NULL, ":"));

	// 3. reading the mesh connectivity to file
	cout << "\nVSet<dim>::ParallelInputFrom: reading finite element mesh..." << endl;
	InBinary(fp);

	// 5. cleaning up
	fp.close();

	cout << "\nVSet<" << dim << ">::ParallelInputFrom: VSet has been successfully read from: ";
	cout << file_name << endl;

	return true;

} // end ParallelInputFrom




/** reads -in polygonal dataset to internal VData container
 
 @todo read 'pmtrl' and property data as well
 
 */
template<uint32_t dim>
bool  VSet<dim>::InputFromTextFile( const char* text_file )
{
	string file(text_file);
	file += ".txt";

	ifstream  ifs(file.c_str());

	if (!ifs.is_open()) return false;

	// reading the first line of the file and echoing it to screen
	char  input_line[200];
	ifs.getline(input_line, 200, '\n');
	cout << "\nVSet<dim>::InputFromTextFile: first line of input file: '";
	cout << file << "'" << endl;
	cout << "\n" << input_line << endl;

	InText(ifs);

	CheckFix();

	cout << "\n\tFile read successfully." << endl;

	return true;

} // end InputFrom(textfile)




/**
       Writes VSet to text file.
       
        @todo write non-scalar property data as well

*/
template<uint32_t dim>
void VSet<dim>::Out( bool print_data_as_well ) const
{
	VData::Out();

	if ( print_data_as_well )
    {
      if ( distance(PmtrlBegin(),PmtrlEnd()) > 0 )
        cout << "\nVSet<"<< dim <<">::Out: material identifiers (rocktypes) for each element stored in VSet:\n";
      // material records
      for ( auto& it : pmtrl_ ) 
        cout << it <<" ";
      // scalar type data
      if ( distance(PropertyValuesBegin(),PropertyValuesEnd()) > 0 )
        cout << "\nproperty records stored in VSet:\n";
      for ( map<string, PropertyData>::const_iterator
            it = property_map_.begin(); it != property_map_.end(); it++ )
        {
          cout << "\n" << (*it).first << endl;
          (*it).second.Out();
        }
    }

} // end Out







    /// writes C++17  code that reproduces a hardwired version of the current VSet
template<uint32_t dim>
void VSet<dim>::OutCPP17( const char* cpp_file ) const
 {
    ofstream  ofs( string(cpp_file) +".cpp" );
    
    //----------------------------FILE HEADER
    ofs <<"// '"<< cpp_file <<"' - cplusplus source code file for the generatioon of a VSet.\n";
    
    ofs <<"\n\nVSet<dim>  vset;";
    
    //----------------------------WRITING THE VDATA
    VData::OutCPP17( ofs );
    
    //----------------------------MATERIAL IDENFIFIERS FOR EACH ELEMENT
    // vector<int32_t> pmtrl(45,1); // matrix
    ofs <<"\nconst int32_t  material_identifier{1};";
    ofs <<"\nvector<int32_t> pmtrl( "<< Elements() <<", material_identifier );";
    //fill( next(pmtrl.begin(),25), next(pmtrl.begin(),31), 2 ); // fine because wrong values will be overwritten next
    ofs <<"\n\nvset.AddPmtrl( pmtrl.begin(), pmtrl.end() );";
    
    //----------------------------NODE & ELEMENT NUMBERS AS PROPERTIES
    // adding node and element numbers for comparisons
    ofs <<"\n\nPropertyData elmt_nums( ELEMENT, SCALAR, 2U );";
    ofs <<"\nelmt_nums.Reserve( vset.Elements() );";
    ofs <<"\n\nfor ( auto i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, static_cast<double>(i) ) );";
    ofs <<"\nvset.AddData( \"element number\" , elmt_nums );";
    // node numbers
    ofs <<"\n\nPropertyData node_nums( NODE, SCALAR, 2U );";
    ofs <<"\nnode_nums.Reserve( vset.Vertices() );";
    ofs <<"\n\nfor ( auto i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, static_cast<double>(i) ) );";
    ofs <<"\nvset.AddData( \"node number\", node_nums );";
 
 } // end OutCPP17






/**
   updates pmtrl and property storage to size changes in VData
*/
template<uint32_t dim>
void VSet<dim>::UpdatePropertyStorage()
 {
    const int32_t  pmtrl_default_value{1};
     if ( pmtrl_.size() != Elements() )
       pmtrl_.resize( Elements(), pmtrl_default_value );
       
    // storage of discretised variables
    if ( !property_map_.empty() ) {
        for ( auto& it : property_map_ ) {
             switch ( it.second.Placement() ) {
               case NODE:
                   if ( Vertices() != it.second.Size() )
                     it.second.Resize( Vertices() );
                 break;
               case ELEMENT:
                   if ( Elements() != it.second.Size() )
                     it.second.Resize( Elements() );
                 break;
               default:
                 cerr <<"\nplacement of '"<< it.first <<"', not handled yet; no resizing done.\n";
             }
          }
    }
    
 } // end UpdatePropertyStorage




/**

Reduces the number of elements in the 'vset' to those identified by their
ID number (1...n) in the supplied set. Then the correspond element
material property data containers are also reduced.

@parameter o_n_elmt_ids is a mapping that allows to retrieve the new element IDs by their old ones.

@note revised by SKM 4/12/20 - to include properties placed on the node as well and warning user about errors for properties with any other placement.

*/
template<uint32_t dim> //             old    new
void VSet<dim>::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids )
{
	if (o_n_elmt_ids.empty())
		throw csmp::Exception(ERROR, "VSet<dim>::ReduceTo:", "new element ID set is empty.");

  std::map<size_t, size_t> o_n_node_ids;
	if ( Elements() > 0 )
	// adjusting 'plist' and 'pfverts' in base class, retrieving which nodes are kept
   VData::ReduceTo( o_n_elmt_ids, o_n_node_ids );
   
  // if there are no properties, all is done already
	if ( property_map_.empty() && pmtrl_.empty() ) return;

  // condensing the pmtrl vector and property maps
  // ---------------------------------------------
  // required vector of elements that are kept
	std::vector<size_t> n_o_elmt_ids;
	n_o_elmt_ids.resize(o_n_elmt_ids.size());
	for ( const auto& o_n : o_n_elmt_ids )
    n_o_elmt_ids[o_n.second] = o_n.first;
  // node vector
	std::vector<size_t> n_o_node_ids;
	n_o_node_ids.resize(o_n_node_ids.size());
	for ( const auto& o_n : o_n_node_ids )
    n_o_node_ids[o_n.second] = o_n.first;
  
  // material indentifiers
  vector<int32_t> new_pmtrl;
  new_pmtrl.reserve( n_o_elmt_ids.size() );
  for ( auto id : n_o_elmt_ids )
    new_pmtrl.push_back( pmtrl_[id] );
  pmtrl_ = new_pmtrl;
  
  // distributed properties
	for ( auto& property : property_map_ ) 
    {
      auto& oldprop = property.second;
      // handling properties with different placements
      if ( oldprop.Placement() == ELEMENT ) {
          PropertyData new_data( oldprop.Placement(), oldprop.Type(), dim );
          new_data.Reserve(n_o_elmt_ids.size());
          for (auto& id : n_o_elmt_ids)
            new_data.PushBackFrom(oldprop, id);
          // reassigning reduced set
          property.second = std::move(new_data);
       }
      else if ( oldprop.Placement() == NODE ) {
          PropertyData new_data( oldprop.Placement(), oldprop.Type(), dim );
          new_data.Reserve(n_o_node_ids.size());
          for (auto& id : n_o_node_ids)
            new_data.PushBackFrom(oldprop, id);
          // reassigning reduced set
          property.second = std::move(new_data);
       }
      else {
          cerr <<"\nERROR: VSet<"<< dim <<">::ReduceTo: placement of property '"<< property.first <<"' not handled. yet.\n";
          cerr <<"\n\tProperties with the placement "<< parsePlacement(oldprop.Placement()) <<" were not transferred correctly\n";
          throw logic_error("VSet<dim>::ReduceTo");
       }
	  }
    
} // end ReduceTo





/** Frees up all the storage in the VSet<dim>; all contained data are deleted.

@section application Application

Use method to economize on memory usage in computations.
*/
template<uint32_t dim>
void VSet<dim>::Erase()
{
	// mesh connnectivity
	VData::Erase();
	// property data
	property_map_.clear();

} // end Erase



/**
   uses the node coordinates to infer the model dimension: if Z-range=zero, dim=2, if Y-range=2, dim=1, else dim=3
*/
template<uint32_t dim>
int32_t  VSet<dim>::MeshDimension( bool check_coordinates ) const
 {
    if ( !check_coordinates ) return dim;
    
    // if they are not empty but all coordinate values are zero, the model has no extent in these directions
    pair<double,double> z_range = Z_Range(), y_range = Y_Range();
    bool z_zero = ( fabs( z_range.second - z_range.first ) <= numeric_limits<double>::epsilon() ) ? true : false; 
    bool y_zero = ( fabs( y_range.second = y_range.first ) <= numeric_limits<double>::epsilon() ) ? true : false; 
    
    if ( z_zero and  y_zero ) return 1;
    if ( z_zero and !y_zero ) return 2;
    
    // if there is nothing that allows a diagnosis
    return UNSPECIFIED;
    
 } // MeshDimension






template class VSet<1U>;
template class VSet<2U>;
template class VSet<3U>;

} // end namespace csmp
