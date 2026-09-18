// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CSMP_highLevelUtilities.h"
#include "compareFloats.h"
#include "Box.h"
#include "VSet.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "ModelTime.h"
#include "TextInterface.h"
#include "Matrix.h"
#include "convertColorToPermeability.h"
#include "Triangulator.h"
#include "MeshManager.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearLineElement.h"

using namespace std;

namespace csmp {

// strings
void replaceWhiteSpaceBy( string& p, char ascii_char )
 {
    if ( !p.empty() )
      for ( uint32_t i{0U}; i<p.size(); i++ )
        if ( p[i] == ' ' || p[i] == '\t' ||
             p[i] == '\n' || p[i] == '\r' ) p[i] = ascii_char;
 }






/**
 * Returns true if file on ifstream is empty.(Aug 2014)
 * @author Julian E. Mindel
 */
bool isInputFileEmpty( ifstream& pFile )
{
   return pFile.peek() == ifstream::traits_type::eof();
}




size_t  renumberElementNodes( vector<Element<1U>*>::iterator first,
                              vector<Element<1U>*>::iterator last )
 {
    assert( first != last );

    set<uint32_t>  node_numbers;
    uint32_t       counts(0);
    
    while ( first != last ) {
         for ( auto nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              pair<set<uint32_t>::iterator,bool>
              sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes





size_t  renumberElementNodes( vector<Element<2U>*>::iterator first,
                              vector<Element<2U>*>::iterator last )
 {
    assert( first != last );

    set<size_t>  node_numbers;
    size_t       counts(0);
    
    while ( first != last ) {
         for ( auto nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              auto sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes




size_t  renumberElementNodes( vector<Element<3U>*>::iterator first,
                              vector<Element<3U>*>::iterator last )
 {
    assert( first != last );

    set<size_t>  node_numbers;
    size_t       counts(0);
    
    while ( first != last ) {
         for ( auto nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); nit++ ) {
              auto sit=node_numbers.insert(counts);
              if ( sit.second == true ) (*nit)->Idx( counts++ );
              else assert( (*nit)->Idx() == (*sit.first) );
           }
         first++;
      }
    
    return node_numbers.size();
    
 } // end renumberElementNodes





void printRangeOfVectorOfVectors( const vector<vector<double> >&  data )
 {
    assert( !data.empty() );
    double vmin(data[0][0]), 
              vmax(data[0][0]);
    
    for ( auto it=data.begin(); it!=data.end(); it++ )
      for ( auto dit=(*it).begin(); dit!=(*it).end(); dit++ ) {
           vmin = min( vmin, (*dit) );
           vmax = max( vmax, (*dit) );
        }
    
    cout <<"\nprintRangeOfVectorOfVectors: Data range: "<< vmin <<" to "<< vmax << endl;
    
 } // end printRangeOfVectorOfVectors



 
void printRangeOf( const vector<pair<double,double> >&  data )
 {
    assert( !data.empty() );
    double vmin(data[0].first), 
              vmax(data[0].second);
    
    for ( auto it=data.begin(); it!=data.end(); it++ ) {
           vmin = min( vmin, (*it).first );
           vmax = max( vmax, (*it).second );
        }
    
    cout <<"\nprintRangeOf: Data range: "<< vmin <<" to "<< vmax << endl;
    
 } // end printRangeOf
 

  

/**
   @TODO document this time-tokenisation and formatting function
*/
#if defined _MSC_VER || defined __MINGW32__
const char* strp_weekdays[] =
{ "sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday" };
const char* strp_monthnames[] =
{ "january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december" };
bool strp_atoi(const char * & s, int & result, int low, int high, int offset)
{
	bool worked = false;
	char* end(nullptr);
	unsigned long num = strtoul(s, &end, 10);
	if (num >= (unsigned long)low && num <= (unsigned long)high)
	{
		result = (int)(num + offset);
		s = end;
		worked = true;
	}
	return worked;
}
char * strptime(const char *s, const char *format, struct tm *tm)
{
	bool working = true;
	while (working && *format && *s)
	{
		switch (*format)
		{
		case '%':
		{
					++format;
					switch (*format)
					{
					case 'a':
					case 'A': // weekday name
						tm->tm_wday = -1;
						working = false;
						for (auto i = 0; i < 7; ++i)
						{
							size_t len = strlen(strp_weekdays[i]);
							if (!strnicmp(strp_weekdays[i], s, len))
							{
								tm->tm_wday = i;
								s += len;
								working = true;
								break;
							}
							else if (!strnicmp(strp_weekdays[i], s, 3))
							{
								tm->tm_wday = i;
								s += 3;
								working = true;
								break;
							}
						}
						break;
					case 'b':
					case 'B':
					case 'h': // month name
						tm->tm_mon = -1;
						working = false;
						for (auto i = 0; i < 12; ++i)
						{
							size_t len = strlen(strp_monthnames[i]);
							if (!strnicmp(strp_monthnames[i], s, len))
							{
								tm->tm_mon = i;
								s += len;
								working = true;
								break;
							}
							else if (!strnicmp(strp_monthnames[i], s, 3))
							{
								tm->tm_mon = i;
								s += 3;
								working = true;
								break;
							}
						}
						break;
					case 'd':
					case 'e': // day of month number
						working = strp_atoi(s, tm->tm_mday, 1, 31, 0);
						break;
					case 'D': // %m/%d/%y
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_mon, 1, 12, -1);
								  if (working && *s == '/')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_mday, 1, 31, 0);
									  if (working && *s == '/')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_year, 0, 99, 0);
										  if (working && tm->tm_year < 69)
											  tm->tm_year += 100;
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'H': // hour
						working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
						break;
					case 'I': // hour 12-hour clock
						working = strp_atoi(s, tm->tm_hour, 1, 12, 0);
						break;
					case 'j': // day number of year
						working = strp_atoi(s, tm->tm_yday, 1, 366, -1);
						break;
					case 'm': // month number
						working = strp_atoi(s, tm->tm_mon, 1, 12, -1);
						break;
					case 'M': // minute
						working = strp_atoi(s, tm->tm_min, 0, 59, 0);
						break;
					case 'n': // arbitrary whitespace
					case 't':
						while (isspace((int)*s))
							++s;
						break;
					case 'p': // am / pm
						if (!strnicmp(s, "am", 2))
						{ // the hour will be 1 -> 12 maps to 12 am, 1 am .. 11 am, 12 noon 12 pm .. 11 pm
							if (tm->tm_hour == 12) // 12 am == 00 hours
								tm->tm_hour = 0;
						}
						else if (!strnicmp(s, "pm", 2))
						{
							if (tm->tm_hour < 12) // 12 pm == 12 hours
								tm->tm_hour += 12; // 1 pm -> 13 hours, 11 pm -> 23 hours
						}
						else
							working = false;
						break;
					case 'r': // 12 hour clock %I:%M:%S %p
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 1, 12, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
									  if (working && *s == ':')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
										  if (working && isspace((int)*s))
										  {
											  ++s;
											  while (isspace((int)*s))
												  ++s;
											  if (!strnicmp(s, "am", 2))
											  { // the hour will be 1 -> 12 maps to 12 am, 1 am .. 11 am, 12 noon 12 pm .. 11 pm
												  if (tm->tm_hour == 12) // 12 am == 00 hours
													  tm->tm_hour = 0;
											  }
											  else if (!strnicmp(s, "pm", 2))
											  {
												  if (tm->tm_hour < 12) // 12 pm == 12 hours
													  tm->tm_hour += 12; // 1 pm -> 13 hours, 11 pm -> 23 hours
											  }
											  else
												  working = false;
										  }
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'R': // %H:%M
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'S': // seconds
						working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
						break;
					case 'T': // %H:%M:%S
					{
								  const char * s_save = s;
								  working = strp_atoi(s, tm->tm_hour, 0, 23, 0);
								  if (working && *s == ':')
								  {
									  ++s;
									  working = strp_atoi(s, tm->tm_min, 0, 59, 0);
									  if (working && *s == ':')
									  {
										  ++s;
										  working = strp_atoi(s, tm->tm_sec, 0, 60, 0);
									  }
								  }
								  if (!working)
									  s = s_save;
					}
						break;
					case 'w': // weekday number 0->6 sunday->saturday
						working = strp_atoi(s, tm->tm_wday, 0, 6, 0);
						break;
					case 'Y': // year
						working = strp_atoi(s, tm->tm_year, 1900, 65535, -1900);
						break;
					case 'y': // 2-digit year
						working = strp_atoi(s, tm->tm_year, 0, 99, 0);
						if (working && tm->tm_year < 69)
							tm->tm_year += 100;
						break;
					case '%': // escaped
						if (*s != '%')
							working = false;
						++s;
						break;
					default:
						working = false;
					}
		}
			break;
		case ' ':
		case '\t':
		case '\r':
		case '\n':
		case '\f':
		case '\v':
			// zero or more whitespaces:
			while (isspace((int)*s))
				++s;
			break;
		default:
			// match character
			if (*s != *format)
				working = false;
			else
				++s;
			break;
		}
		++format;
	}
	return (working ? (char *)s : 0);
}
#endif // defined _MSC_VER || defined __MINGW32__





/// Utility that tokenises string into substrings using the supplied delimiter(s).
vector<string> splitString( string str, char delimiter )
{
  size_t pos = 0U;
  string token, s = str;
  vector<string> items;
  while ( (pos = s.find( delimiter )) != string::npos ) {
    token = s.substr( 0, pos );
    items.push_back( token );
    s.erase( 0, pos + 1 );
  }
  if ( !s.empty() ) items.push_back( s );
  return items;
}








 /**
     Finds all possible combinations of single or multiple values in the input vector,
     after sorting it and making it unique.
     
      @param samples refers to the subset of values for which unique combinations shall be found.
      @param combinations will store the uniqe combinations that were found
      
      @attention combinations is not equal to permutations.
      
      @author SKM (modified from example on stackoverflow)
      @date 4/10/2021
      
 */
template<typename intType>
size_t createUniqueCombinations( vector<intType>& sequence, uint32_t samples,
                                 deque<vector<intType> >& combinations )
 {
    // checking the input
    if ( sequence.empty() ) return 0U;
    if ( samples > sequence.size() ) {
        cerr <<"\ncreateUniqueCombinations: can't combine more numbers than are in the input vector.\n";
        return 0U;
      }
    sort( sequence.begin(), sequence.end() );
    sequence.erase( unique( sequence.begin(), sequence.end() ), sequence.end() );
    const size_t N{sequence.size()};
    
    // generating combinations by selectively sampling sequence using 011.. pattern in bitmap
    string bitmask(samples, 1); // generating leading 1's
    bitmask.resize(N, 0);       // adding (N - samples) trailing 0's
 
    do {
        combinations.push_back( vector<intType>{} );
        combinations.back().reserve( samples );
        for ( size_t i{0U}; i < N; ++i ) { // [0..N-1] integers
             if ( bitmask[i] == 1 )
               combinations.back().push_back( sequence[i] );
          }
      }
    while ( prev_permutation( bitmask.begin(), bitmask.end() ) );
    
    return combinations.size();
    
} // end createUniqueCombinations
 
template size_t createUniqueCombinations( vector<uint32_t>&, uint32_t, deque<vector<uint32_t> >& );
template size_t createUniqueCombinations( vector<size_t>&, uint32_t, deque<vector<size_t> >& );
template size_t createUniqueCombinations( vector<int>&, uint32_t, deque<vector<int> >& );
template size_t createUniqueCombinations( vector<int64_t>& sequence, uint32_t, deque<vector<int64_t> >& combinations );



/// usage example for createUniqueCombinations()
static void test_createUniqueCombinations()
 {
   vector<int64_t>          sequence{0,123,20,43,17,5,8};
   const uint32_t           samples{2};
   deque<vector<int64_t> >  combinations;
   
   size_t n_combinations = createUniqueCombinations( sequence, samples, combinations );
   
   // printing the results
   cout <<"\nmain: combinations: "<< n_combinations <<"\n";
   for ( const auto& i : combinations ) {
        cout <<"\n\t";
        for ( const auto& j : i )
          cout <<" "<< j;
     }
   cout << endl << endl;
 }
 
 
// is the executable compiled for the Rosetta X86 replacement environment on Apple Silicon
#if defined(__APPLE__)
#include <sys/sysctl.h>

inline bool running_under_Rosetta()
{
    int translated = 0;
    size_t size = sizeof(translated);
    if (sysctlbyname("sysctl.proc_translated",
                     &translated,
                     &size,
                     nullptr,
                     0) == 0)
    {
        return translated == 1;
    }
    return false;
}
#endif
 
 
} // end namespace csmp
