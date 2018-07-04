#include "Standard_IO_Handler.h"

using namespace std;

namespace csmp {

/**Default constructor

Establishes  Standard_IO_Handler-protocol.txt as protocol file
@author S.K. Matthai
*/
Standard_IO_Handler::Standard_IO_Handler()
  : protocol_file_("Standard_IO_Handler-protocol.txt")
  {
  }
  
/** Custom Constructor 

establishes protocol filename based on argument
@param doc_file Filename for protocol file
*/
Standard_IO_Handler::Standard_IO_Handler( const char* doc_file )
  : protocol_file_(doc_file)
  {
  }

///Default destructor
Standard_IO_Handler::~Standard_IO_Handler() 
 {
 }



/** Prompts the user for yes/no input

@param question Question to be asked
@returns true if user replies 'yes'
*/
bool Standard_IO_Handler::YesNo( const char question[150] )
  {
    bool answer = false;
    char  c;
    int   i;
    
    if ( (c=static_cast<char>(cin.peek())) == '\n' ) cin.get( c );
    
    /* USER YES/NO PROMPT */
    quest : { cout<<"\n\n"<< question << " [y/n] ";
              cout.flush();
              cin.get( c );
              i = (int32)( c );
              if ( i == 121 ) answer = true;
              if ( i == 110 ) answer = false;
              if ( i != 121 && i != 110 ) 
                {
                   cin.get( c );
                   goto quest; 
                }
            } // END QUESTION

    string  communication(question);
    communication += " -> ";
    if ( c == 'y' ) communication += "true.";
    else            communication += "false.";
    input_output.push_back(communication);
    
    return ( answer );

  } // end yes_no





/** Documents outputs in logfile

@param info The text to be displayed if file was opened successfully

*/
void Standard_IO_Handler::Protocol( const char* info )
 {
    ofstream  ofs;
    string    s(protocol_file_);

    ofs.open ( s.c_str(), ios::out|ios::trunc );
    if ( !ofs.is_open() )
      {
         cout <<"\nStandard_IO_Handler::Protocol: ";
         cout <<"Cannot open documentation file !"<< endl;
         cout <<"\nPlease enter new filename: ";
         cin >> s;
       }
    ofs.open ( s.c_str(), ios::out|ios::trunc );
    if ( !ofs.is_open() ) 
      throw invalid_argument( (s="Standard_IO_Handler::Protocol: cannot open protocol file.") );
    else
    ofs << info << endl; 

    ofs.close();
   
 } // end protocoll





/** Documents outputs in logfile

@param info The text to be displayed if file was opened successfully
@param parameter Number to be displayed after successfully opening file

*/
void Standard_IO_Handler::Protocol( const char* info, double64 parameter )
  {
    ofstream  ofs;
    string    s(protocol_file_);

    ofs.open ( s.c_str(), ios::out|ios::trunc );
    if ( !ofs.is_open() )
      {
         cout <<"\nStandard_IO_Handler::Protocol: ";
         cout <<"Cannot open documentation file !"<< endl;
         cout <<"\nPlease enter new filename: ";
         cin >> s;
       }
    ofs.open ( s.c_str(), ios::out|ios::trunc );
    if ( !ofs.is_open() ) 
      throw invalid_argument( (s="Standard_IO_Handler::Protocol: cannot open protocol file.") );
    else
    ofs << info <<"\t\t" << parameter << endl; 

    ofs.close();
   } // end protocoll






/** Records user choice to the protocol cache

@param question The question that is to be prompted
@return The user input
*/

double64 Standard_IO_Handler::RecordChoice( const char* question )
 {
    cout <<"\n"<< question <<" ";
    double64 result;
    cin >> result;

    string  communication(question);
    char    num[50];
    sprintf( num, "%lf", result );
    communication += " -> ";
    communication += num;
    input_output.push_back(communication);
    
    return result;
 }                               





 /** Records user choice to the protocol cache

@param question The question that is to be prompted
@return The user input
*/

long Standard_IO_Handler::RecordIntChoice( const char* question )
 {
    cout <<"\n"<< question <<" ";
    long   result;
    cin >> result;

    string  communication(question);
    char    num[50];
    sprintf( num, "%ld", result );
    communication += " -> ";
    communication += num;
    input_output.push_back(communication);
    
    return result;
 }                         




/** Records logical user choice to the protocol cache

@param question The question that is to be prompted
@return The boolean user input, true if 'yes'

*/
bool Standard_IO_Handler::RecordLogicalChoice( const char* question )
 {
    bool yes_or_no = YesNo(question);

    string  communication(question);
    communication += " -> ";
    if ( yes_or_no ) communication += "true.";
    else             communication += "false.";
    input_output.push_back(communication);
    
    return yes_or_no;
 }                               





/** Records literal user choice to the protocol cache

@param question The question that is to be prompted
@return The user's response

*/
string Standard_IO_Handler::RecordLiteralChoice( const char* question )
 {
    cout <<"\n"<< question <<" ";
    string result;
    cin >> result;

    string  communication(question);
    communication += " answer -> ";
    communication += result;
    input_output.push_back(communication);
    
    return result;
 }                               


/** Records argument to the protocol cache

@param info The information that is to be protocolled
*/

void  Standard_IO_Handler::RecordInformation( const string& info )
 {
    input_output.push_back(info);
 }





/** CSMP output feature

Out() returns all the internal data that Index needs to access 
variables. These include the type and placement of a variable in the 
finite element mesh and the associated data storage index.

*/
void Standard_IO_Handler::Out() const
 {
    ofstream  ofs( protocol_file_.c_str() );
    ofs <<"\nStandard_IO_Handler::Out: File '"<< protocol_file_ <<"'";
    
    ofs.setf(ios::scientific);
    long prec = ofs.precision(15U);
   
    for ( list<string>::const_iterator 
          it=input_output.begin(); it!=input_output.end(); it++ )
      ofs <<"\n"<< (*it);
    ofs << endl;

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
 }


} // end namespace csmp




