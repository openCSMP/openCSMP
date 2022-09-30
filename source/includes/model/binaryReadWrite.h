#ifndef CSMP_READ_WRITE_H
#define CSMP_READ_WRITE_H

#include "CSMP_definitions.h"
#include "PropertyDatabase.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/**
@file binaryReadWrite.h
*/

constexpr int CSMP_BINARY_FILE_HDR_SIZE(8); // JCK header for data record in file

/**
   Helper class for reading fixed-size strings from binary files
   
      defines strings that are used to read the possible datablocks that are abbreviated as
        VSETCORD
        VSETPELT
        VSETPLST
        VSETPFVT
        VSETBFLG
        VSETHEDR
        VSETCONN
        VSETFOTR
        VSETMTRL - new for pmtrl record
        
     BoundaryInterface
        BNDFHEDR - file header
        BOUNDARY - object description
        ONE_BDRY
        BOUNDVAR - variable datablock
        BNDFFOTR - file footer
 
    SplitBoundaryInterface
        SBDFHEDR
        SPLITBDRY
        ONE_BDRY
        SBDRYVAR
        SBDFFOTR

    RegionInterface
        REGGHEDR
        UNIQREGN
        ONE_REGN
        NONUREGN
        MODELVARS
        REGFFOTR
        
    NodeManifolds
        NDMFHEDR
        MANIFLDS
*/
class BinaryFileSectionRead
{
public:
	BinaryFileSectionRead( std::fstream& fp, const char* header );
	~BinaryFileSectionRead();

private:
	char hdr_[CSMP_BINARY_FILE_HDR_SIZE + 1]; // header JCK - string + null terminator
	std::fstream& fp_;  
};


/// Helper class to write sections to a binary file
class BinaryFileSectionWrite {
public:
	BinaryFileSectionWrite( std::fstream& fp, const char* header );
	~BinaryFileSectionWrite();

private:
	std::fstream& fp_;
};


/**
@addtogroup CSMPglobalFunctions
*/

/// checks whether the integer value fits within the range of a size_t  and returns it
size_t checkContainerSize( std::fstream& fp );

template<class T>
bool binaryFileWrite(std::fstream& fp, const std::vector<T>& stl_ctner);

template<class T>
bool binaryFileRead(std::fstream& fp, std::vector<T>& stl_ctner);

template<typename T>
bool binaryFileWrite(std::fstream& fp, const std::deque<T>& stl_ctner);

template<typename T>
bool binaryFileRead(std::fstream& fp, std::deque<T>& stl_ctner);

template<class T>
bool binaryFileWrite(std::fstream& fp, const std::deque<std::vector<T> >& stl_ctner);

template<class T>
bool binaryFileRead(std::fstream& fp, std::deque<std::vector<T> >& stl_ctner);

// maps

template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::map<M, T>& stl_ctner);

template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::unordered_map<M, T>& stl_ctner);

template<class M, class T>
bool binaryFileRead(std::fstream& fp, std::map<M, T>& stl_ctner);

// maps of vectors

template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::map<M, std::vector<T> >& stl_ctner);

template<class M, class T>
bool binaryFileRead(std::fstream& fp, std::map<M, std::vector<T> >& stl_ctner);


// character strings

bool binaryFileWrite( std::fstream& fp, const std::string& );

bool binaryFileRead( std::fstream& fp, std::string& );

bool binaryFileWrite( std::fstream& fp, const char* str);

bool binaryFileRead( std::fstream& fp, char str[]);


// OUTPUT OF DISCRETISED CSMP VARIABLES TO FILE

/// additions by P. Lang (2012); SKM @todo explain what this is good for
namespace femDataOutputDispatch {

template<class Var>
inline void initVariable( csmp::Index, Var& ) {}

template<>
inline void initVariable( csmp::Index key, ArrayVariable& var )
{ var.Resize( key.dataDepth ); }

template<>
inline void initVariable( csmp::Index key, FlaggedArrayVariable& var )
{ var.Resize( key.dataDepth ); }

} // femDataOutputDispatch



  /// domain (Model, Region...) variables binary IO
template<class V, class D, uint32_t dim>
bool variablesOut( std::fstream& fp, const D& domain, const PropertyDatabase<dim>& pref, VARIABLE_TYPE vtype )
{
  size_t vcount( pref.VariableCount( domain.Placement(), vtype ) );
  fp.write( (char*)&vcount, sizeof( size_t ) );
  std::set<std::string> propList;
  pref.ListVariables( domain.Placement(), vtype, propList );
  for ( std::set<std::string>::const_iterator it( propList.begin() ); it != propList.end(); ++it )
  {
    V var;
    Index key( pref.StorageKey( it->c_str() ) );
    femDataOutputDispatch::initVariable( key, var );
    domain.Read( key, var );
    binaryFileWrite( fp, it->c_str() );
    if ( !var.Out( fp ) )
      return false;
  }
  return true;
}


template<class V, class D, uint32_t dim>
bool variablesIn( std::fstream& fp, D& domain, const PropertyDatabase<dim>& pref, VARIABLE_TYPE )
{
  size_t vcount( -1 );
  fp.read( reinterpret_cast<char*>(&vcount), sizeof( size_t ) );
  for ( size_t i( 0 ); i < vcount; ++i )
  {
    V var;
    char propName[NAME_STRING];
    binaryFileRead( fp, propName );
    if ( pref.IsDefined( propName ) ) {
      Index key( pref.StorageKey( propName ) );
      femDataOutputDispatch::initVariable( key, var );
      if ( !var.In( fp ) )
        return false;
      domain.Store( key, var );
    }else{
      if ( !var.In( fp ) )
        return false;
    }
  }
  return true;
}


template<class D, uint32_t dim>
bool domainVariablesOut( std::fstream& fp, const D& domain, const PropertyDatabase<dim>& pref )
{
  if ( !variablesOut<ScalarVariable>( fp, domain, pref, SCALAR ) )
    return false;
  if ( !variablesOut<VectorVariable<dim> >( fp, domain, pref, VECTOR ) )
    return false;
  if ( !variablesOut<TensorVariable<dim> >( fp, domain, pref, TENSOR ) )
    return false;
  if ( !variablesOut<ArrayVariable>( fp, domain, pref, ARRAY ) )
    return false;
  if ( !variablesOut<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY ) )
    return false;
  return true;
}


template<class D, uint32_t dim>
bool domainVariablesIn( std::fstream& fp, D& domain, const PropertyDatabase<dim>& pref )
{
  if ( !variablesIn<ScalarVariable>( fp, domain, pref, SCALAR ) )
    return false;
  if ( !variablesIn<VectorVariable<dim> >( fp, domain, pref, VECTOR ) )
    return false;
  if ( !variablesIn<TensorVariable<dim> >( fp, domain, pref, TENSOR ) )
    return false;
  if ( !variablesIn<ArrayVariable>( fp, domain, pref, ARRAY ) )
    return false;
  if ( !variablesIn<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY ) )
    return false;
  return true;
}


/**
    selective variable reader, that extracts only those variables from file whose names are contained in the target set
    @author SKM
    @date 6/9/2021
 */
template<class V, class D, uint32_t dim>
bool selectedVariablesIn( std::fstream& fp, D& domain,
                          const PropertyDatabase<dim>& pref,
                          VARIABLE_TYPE, const std::set<std::string>& selection )
{
  size_t vcount( -1 );
  fp.read( reinterpret_cast<char*>(&vcount), sizeof( size_t ) );
  for ( size_t i( 0 ); i < vcount; ++i )
    {
      V var;
      char propName[NAME_STRING];
      binaryFileRead( fp, propName );
      if ( selection.find( propName ) != selection.end() && pref.IsDefined( propName ) ) {
        Index key( pref.StorageKey( propName ) );
        femDataOutputDispatch::initVariable( key, var );
        if ( !var.In( fp ) )
          return false;
        domain.Store( key, var );
      }else{
        if ( !var.In( fp ) )
          return false;
      }
    }
  return true;
}


/**
    selective variable reader, that extracts only those variables from file whose names are contained in the target set
    @author SKM
    @date 6/9/2021
 */
template<class D, uint32_t dim>
bool selectedDomainVariablesIn( std::fstream& fp, D& domain,
                                const PropertyDatabase<dim>& pref,
                                const std::set<std::string>& selection )
{
  if ( !selectedVariablesIn<ScalarVariable>( fp, domain, pref, SCALAR, selection ) )
    return false;
  if ( !selectedVariablesIn<VectorVariable<dim> >( fp, domain, pref, VECTOR, selection ) )
    return false;
  if ( !selectedVariablesIn<TensorVariable<dim> >( fp, domain, pref, TENSOR, selection ) )
    return false;
  if ( !selectedVariablesIn<ArrayVariable>( fp, domain, pref, ARRAY, selection ) )
    return false;
  if ( !selectedVariablesIn<FlaggedArrayVariable>( fp, domain, pref, FLAGGEDARRAY, selection ) )
    return false;
  return true;
}






/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL deque to a binary file. The data segment is preceded by an
'unsigned long' number which determines the number of objects which are
written to file.

@section arguments Input Arguments

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL deque which will be written to file.

@return returns a boolean indicating whether all data have been
written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite().

@section application Application

To efficiently write deque data to a binary file.

@section messages Messages

If the file pointer is invalid, method will quit, reporting an error.
*/
template<class T>
bool binaryFileWrite( std::fstream& fp, const std::vector<T>& stl_ctner )
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(vector): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileWrite: ERROR: invalid file pointer." << std::endl;
      return false;
    }

	const size_t	bytes = sizeof(T);
	const size_t	elements(stl_ctner.size());

	// writing the size of the object
	fp.write( reinterpret_cast<const char*>(&elements), sizeof(size_t) );

	// writing all elements
  if ( elements > 0 )
    fp.write( reinterpret_cast<const char*>(&stl_ctner[0]), bytes * elements );

	return true;
}


/// deque version
template<typename T>
bool binaryFileWrite( std::fstream& fp, const std::deque<T>& stl_ctner )
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(deque): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileWrite( deque<T>& ): ERROR: invalid file pointer." << std::endl;
      return false;
    }

	const size_t bytes = sizeof(T);
	size_t       elements = stl_ctner.size();

	// writing the size of the object
	fp.write( reinterpret_cast<const char*>(&elements), sizeof(size_t) );

	// writing all elements
  if ( elements > 0 )
	  for ( typename std::deque<T>::const_iterator it = stl_ctner.begin(); it != stl_ctner.end(); it++ )
		  fp.write(reinterpret_cast<const char*>(&(*it)), bytes);

	return true;
}


/**

Reads data from a binary file into an STL deque which is erased before
adding the data, if it already contained any data. How many
records are read is specified by an 'unsigned long' number which precedes
the dataset.

@section arguments Input Arguments

A pointer to a binary file opened in read binary mode ("rb"), and a
reference to an STL deque.

@return returns the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the
file pointer is invalid or less data are read, the function returns false.

@section implementation Implementation

Uses the ANSI C function fread(). The deque which will hold the data
is erased if it is not empty, and storage is reserved for the new number
of elements which is read from file. Then the elements are read and the
deque is returned.

@section application Application

Efficiently read data into STL deques.

@section messages Messages

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this
number does not match the number of records which were actually read.
*/
template<class T>
bool binaryFileRead( std::fstream& fp, std::vector<T>& stl_ctner )
{
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileRead( vector<T> ): ERROR: invalid file pointer." << std::endl;
      return false;
    }
    
	const size_t bytes = sizeof(T);
	size_t       counter(0);
	T            val;

	// read size of the record and check it 
  const size_t elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
  if ( elements > 0 ) {
        if (!stl_ctner.empty())
          stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());
        stl_ctner.reserve(elements);
        // writing all elements
        for ( auto i = 0; i<elements; i++ ) {
          // counting the successfully read elements
          fp.read( reinterpret_cast<char*>(&val), bytes );
          counter += 1;
          stl_ctner.push_back(val);
      }
    }
	if (counter != elements) {
      std::cerr << "\nbool binaryFileRead: ERROR: incorrect number of records were read: ";
      std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counter << std::endl;
      return false;
    }
    
	return true;
}

/* ALTERNATIVE
std::vector<uint8_t> read_vector_from_disk(std::string file_path)
  {
      std::ifstream instream(file_path, std::ios::in | std::ios::binary);
      std::vector<uint8_t> data((std::istreambuf_iterator<char>(instream)), std::istreambuf_iterator<char>());
      return data;
  }
*/



/// deque version
template<typename T>
bool binaryFileRead( std::fstream& fp, std::deque<T>& stl_ctner )
{
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileRead: ERROR: invalid file pointer." << std::endl;
      return false;
    }
	const size_t bytes = sizeof(T);
  size_t counter(0);

	// read size of the record and assert this 
  size_t elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
	if ( elements > 0 )
    {
      if (!stl_ctner.empty())
        stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());
        
    	T val;
      // writing all elements
      for ( size_t i{0U}; i<elements; i++ )
        {
          // counting the successfully read elements
          fp.read( reinterpret_cast<char*>(&val), bytes );
          counter += 1;
          stl_ctner.push_back(val);
        }
    }
	if ( counter != elements )
    {
      std::cerr << "\nbool binaryFileRead: ERROR: incorrect number of records were read: ";
      std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counter << std::endl;
      return false;
    }
	return true;
}



/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL deque<vector<T> > to a binary file. The data segment is preceded
by an 'unsigned long' number which determines the number of vector objects
which are written to file.

@section arguments Input Arguments

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL vector of vectors which will be written to
file.

@return returns a boolean indicating whether all vector data have
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite().

@section application Application

To efficiently write vector data to a binary file.

@section messages Messages

If the file pointer is invalid, method will quit, reporting an error.
*/
template<class T>
bool binaryFileWrite( std::fstream& fp, const std::deque<std::vector<T> >& stl_ctner )
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(deque<vector>): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileWrite(const deque<vector<T> >&): ";
      std::cerr << "ERROR: invalid file pointer." << std::endl;
      return false;
    }

	// writing the number of vector objects
	const int64_t   elements(stl_ctner.size());
	fp.write( reinterpret_cast<const char*>(&elements), sizeof(size_t) );

	// writing all elements
  if ( elements > 0 )
    for ( typename std::deque<std::vector<T> >::const_iterator
          it = stl_ctner.begin(); it != stl_ctner.end(); it++ )
      binaryFileWrite(fp, (*it) );

	return true;
}



/**

Reads a vector of STL vector data from a binary file. The STL container
is erased before adding the data, if it already contained any data. How many
records are read is specified by an 'unsigned long' number which precedes
the dataset.

@section arguments Input Arguments

A pointer to a binary file opened in read binary mode ("rb"), and a
reference to an STL vector of vector.

@return returns the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the
file pointer is invalid or less data are read, the function returns false.

@section implementation Implementation

Uses the ANSI C function fread(). The vector<vector<T> > which will hold the
data is erased if it is not empty. Storage is reserved for the new number
of elements which is read from file. Then the elements are read and the
container is returned.

@section application Application

Efficiently read data into STL vector of vectors.

@section messages Messages

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this
number does not match the number of records which were actually read.
*/
template<class T>
bool binaryFileRead(std::fstream& fp, std::deque<std::vector<T> >& stl_ctner)
{
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileRead(deque<vector<T> >&): ";
      std::cerr << "ERROR: invalid file pointer." << std::endl;
      return false;
    }

	// 1. reading number of vector records and assert this reading
  const int64_t  elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
	size_t counter(0);
	if ( elements > 0U ) {
      if (!stl_ctner.empty())
        stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());
      //stl_ctner.reserve( elements );
      std::vector<T>  val;
      // 2. reading all the vector records
      for ( auto i = 0; i<elements; i++ ) {
        // counting the successfully read elements
        if (binaryFileRead(fp, val)) {
          counter++;
          stl_ctner.push_back(val);
        }
     }
	}
	if (counter != elements) {
      std::cerr << "\nbool binaryFileRead(deque<vector<T> >&): ";
      std::cerr << "ERROR: incorrect number of records were read: ";
      std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counter << std::endl;
      return false;
    }

	return true;
}



/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL map to a binary file. The data segment is preceded
by an 'unsigned long' number which determines the number of objects
which are written to file. Each record is preceded by the
corresponding map key.

@section arguments Input Arguments

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL map which will be written to
file.

@return a boolean indicating whether all data have
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). NOTE that the map key must
not contain any dynamically allocated data. Otherwise these data are
sliced of the binary record.

@section application Application

To efficiently write map data to a binary file.

@section messages Messages

If the file pointer is invalid, method will quit, reporting an error.
*/
template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::map<M, T>& stl_ctner)
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(map): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open())
    {
      std::cerr << "\nbool binaryFileWrite(const map<M,T,less<M> >&): ";
      std::cerr << "ERROR: invalid file pointer." << std::endl;
      return false;
    }

	const size_t bytesM = sizeof(M);
	const size_t bytesT = sizeof(T);
	size_t       elements = stl_ctner.size();
	M            key;
	T            val;

	// writing the number of vector objects
	fp.write( reinterpret_cast<const char*>(&elements), sizeof(size_t) );

	// writing all key-value pairs
  if ( elements > 0 )
    for ( typename std::map<M, T>::const_iterator it = stl_ctner.begin(); it != stl_ctner.end(); it++)
      {
        key = (*it).first;
        val = (*it).second;
        fp.write( reinterpret_cast<const char*>(&key), bytesM );
        fp.write( reinterpret_cast<const char*>(&val), bytesT );
      }
	return true;
}




/**

Reads an STL map from a binary file. The map
is erased before adding the data if it already contains data. How many
records are read is specified by an 'unsigned long' number which precedes
the dataset.

@section arguments Input Arguments

A pointer to a binary file opened in read binary mode ("rb"), and a
reference to an STL map.

@return binaryFileRead() returns the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the
file pointer is invalid or less data are read, the function returns false.

@section implementation Implementation

Uses the ANSI C function fread(). The map which will hold the
data is erased if it is not empty. Storage is reserved for the new number
of elements which is read from file. Then the elements are read and the
container is returned.

@section application Application

Efficiently read data into an STL map.

@section messages Messages

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this
number does not match the number of records which were actually read.
*/
template<class M, class T>
bool binaryFileRead(std::fstream& fp, std::map<M, T>& stl_ctner)
{
	if (!stl_ctner.empty())
		stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());

	if (!fp.is_open())
	{
		std::cerr << "\nbool binaryFileRead(map<M,T>&): ";
		std::cerr << "ERROR: invalid file pointer." << std::endl;
		return false;
	}

	const size_t bytesM = sizeof(M), bytesT = sizeof(T);
  size_t  counterM(0), counterT(0);
	M       key;
	T       val;

	// 1. read number of record in the map and assert reading
  const int64_t  elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
	if (elements > 0) {
		// 2. reading all map records
		for ( auto i = 0; i<elements; i++ )
		{
			// reading key
			if (fp.read( reinterpret_cast<char*>(&key), bytesM )) counterM++;
			// reading value
			if (fp.read( reinterpret_cast<char*>(&val), bytesT )) counterT++;
			// storing value in map after key
			stl_ctner[key] = val;
		}
	}
	if (counterM != elements || counterT != elements) {
      std::cerr << "\nbool binaryFileRead(vector<map<M,T>&): ";
      std::cerr << "ERROR: incorrect number of map records were read: ";
      std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counterM << std::endl;
      return false;
    }
	return true;
}


/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL unordered_map to a binary file. The data segment is preceded
by an 'unsigned long' number which determines the number of objects
which are written to file. Each record is preceded by the
corresponding map key.

@section arguments Input Arguments

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL map which will be written to
file.

@return a boolean indicating whether all data have
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). NOTE that the map key must
not contain any dynamically allocated data. Otherwise these data are
sliced of the binary record.

@section application Application

To efficiently write map data to a binary file.

@section messages Messages

If the file pointer is invalid, method will quit, reporting an error.
*/
template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::unordered_map<M, T>& stl_ctner )
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(unordered map): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open())
	{
		std::cerr << "\nbool binaryFileWrite(const unordered_map<M,T,less<M> >&): ";
		std::cerr << "ERROR: invalid file pointer." << std::endl;
		return false;
	}
	const size_t elements = stl_ctner.size();
	const size_t bytesM = sizeof(M);
	const size_t bytesT = sizeof(T);
	M      key;
	T      val;

	// writing the number of vector objects
	fp.write(reinterpret_cast<const char*>(&elements), sizeof(size_t));

	// writing all key-value pairs
  if ( elements > 0 )
    for ( typename std::unordered_map<M, T>::const_iterator it = stl_ctner.begin(); it != stl_ctner.end(); it++)
      {
        key = (*it).first;
        val = (*it).second;
        fp.write( reinterpret_cast<const char*>(&key), bytesM );
        fp.write( reinterpret_cast<const char*>(&val), bytesT );
      }
	return true;
}




/**

Reads an STL unordered_map from a binary file. The map
is erased before adding the data if it already contains data. How many
records are read is specified by an 'unsigned long' number which precedes
the dataset.

@section arguments Input Arguments

A pointer to a binary file opened in read binary mode ("rb"), and a
reference to an STL map.

@return binaryFileRead() returns the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the
file pointer is invalid or less data are read, the function returns false.

@section implementation Implementation

Uses the ANSI C function fread(). The map which will hold the
data is erased if it is not empty. Storage is reserved for the new number
of elements which is read from file. Then the elements are read and the
container is returned.

@section application Application

Efficiently read data into an STL map.

@section messages Messages

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this
number does not match the number of records which were actually read.
*/
template<class M, class T>
bool binaryFileRead(std::fstream& fp, std::unordered_map<M, T>& stl_ctner)
{
	if (!fp.is_open())
    {
      std::cerr << "\nbool binaryFileRead(unordered_map<M,T>&): ";
      std::cerr << "ERROR: invalid file pointer." << std::endl;
      return false;
    }
	const size_t bytesM = sizeof(M), bytesT = sizeof(T);
	size_t  counterM(0), counterT(0);
	 M      key;
	 T      val;

	// 1. read number of record in the map and assert reading
  const size_t elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
	if (elements > 0) {
      if (!stl_ctner.empty())
        stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());
      // 2. reading all map records
      for ( auto i = 0; i<elements; i++)
        {
          // reading key
          if (fp.read( reinterpret_cast<char*>(&key), bytesM)) counterM++;
          // reading value
          if (fp.read( reinterpret_cast<char*>(&val), bytesT)) counterT++;
          // storing value in map after key
          stl_ctner[key] = val;
        }
    }
	if (counterM != elements || counterT != elements) {
		std::cerr << "\nbool binaryFileRead(vector<map<M,T>&): ";
		std::cerr << "ERROR: incorrect number of map records were read: ";
		std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counterM << std::endl;
		return false;
	}
	return true;
}




/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL map of vector<T> to a binary file. The data segment is preceded
by an 'size_t' number which determines the number of vector objects
which are written to file. Each vector record is preceded by the
corresponding map key.

@section arguments Input Arguments

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL map of vectors which will be written to
file.

@return a boolean indicating whether all map data have
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). NOTE that the map key must
not contain any dynamically allocated data. Otherwise these data are
sliced of the binary record.

@section application Application

To efficiently write a map of vector data to a binary file.

@section messages Messages

If the file pointer is invalid, method will quit, reporting an error.
*/
template<class M, class T>
bool binaryFileWrite(std::fstream& fp, const std::map<M, std::vector<T> >& stl_ctner)
{
#ifdef DEBUG
  if ( stl_ctner.empty() )
    std::cerr <<"\nbinaryFileWrite(map<vector>>): WARNING: container is empty."<< std::endl;
#endif
	if (!fp.is_open()) {
		std::cerr << "\nbool binaryFileWrite(const map<M,vector<T>,less<M> >&): ";
		std::cerr << "ERROR: invalid file pointer." << std::endl;
		return false;
	}

	const size_t elements = stl_ctner.size();
	const size_t bytes = sizeof(M);
	M            val;

	// writing the number of vector objects
	fp.write(reinterpret_cast<const char*>(&elements), sizeof(size_t));

	// writing all elements
  if ( elements > 0 )
    for ( typename std::map<M, std::vector<T> >::const_iterator it = stl_ctner.begin(); it != stl_ctner.end(); it++)
      {
        val = (*it).first;
        fp.write( reinterpret_cast<const char*>(&val), bytes );
        binaryFileWrite(fp, (*it).second );
      }
	return true;
}



/**

Reads a map of STL vectors and their keys from a binary file. The map
is erased before adding the data if it already contained any data. How many
records are read is specified by an 'size_t' number which precedes
the dataset.

@section arguments Input Arguments

A pointer to a binary file opened in read binary mode ("rb"), and a
reference to an STL map of vectors.

@return the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the
file pointer is invalid or less data are read, the function returns false.

@section implementation Implementation

Uses the ANSI C function fread(). The map which will hold the
data is erased if it is not empty. Storage is reserved for the new number
of elements which is read from file. Then the elements are read and the
container is returned.

@section application Application

Efficiently read data into an STL map of vectors.

@section messages Messages

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this
number does not match the number of records which were actually read.
*/
template<class M, class T>
bool binaryFileRead(std::fstream& fp, std::map<M, std::vector<T> >& stl_ctner)
{
	if (!fp.is_open()) {
      std::cerr << "\nbool binaryFileRead(map<M,vector<T>,less<M> >&): ";
      std::cerr << "ERROR: invalid file pointer." << std::endl;
      return false;
    }
	const size_t   bytes(sizeof(M));
	size_t         counter(0);
	std::vector<T> val;
	M              key;

	// 1. read number of record in the map and assert reading
  const size_t elements = checkContainerSize( fp );

	// TODO: write the elements in one go, rather than one by one
	if (elements > 0) {
    if (!stl_ctner.empty())
      stl_ctner.erase(stl_ctner.begin(), stl_ctner.end());
		// 2. reading all map records
		for ( auto i = 0; i<elements; i++)
      {
        // reading key
        fp.read( reinterpret_cast<char*>(&key), bytes);
        // reading vector and counting the successfully read records
        if (binaryFileRead(fp, val)) {
          counter++;
          stl_ctner[key] = val;
        }
      }
	}
	if (counter != elements) {
      std::cerr << "\nbool binaryFileRead(vector<map<M,vector<T> >&): ";
      std::cerr << "ERROR: incorrect number of map records were read: ";
      std::cerr << "\nIndicated number: " << elements << ", actual number read: " << counter << std::endl;
      return false;
    }
	return true;
}

/**
@}
*/

} // csmp


#endif
