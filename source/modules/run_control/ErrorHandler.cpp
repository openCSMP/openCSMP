#include "ErrorHandler.h"
#include "license.h"

using namespace std;

namespace csmp {

/**

Default constructor initializes message counters to zero, sets the maximum
error limit to a gratuitous 100 such that you can stumble quite far into a 
first draft of a simulation, and sets the message reporting mode to verbose.
 
@section application Application

The default constructor is the only constructor for the error handler and 
it is called before the simulation enters main(). 
*/
ErrorHandler::ErrorHandler()
   : verbose(VERBOSE), debug(false), errors(0), warnings(0), total(0), error_limit(10000)
 {
    if ( !(aus.duration==0) and aus.secs_remaining <= 0 ) {
         cerr <<"\n\nErrorHandler: CSMP++ "<< parse(aus.duration);
         cerr <<" license has expired since "<< -aus.secs_remaining <<" secs.";
         cerr <<" Contact: 'stephan.matthai@unileoben.ac.at' for an extension.\n";
         terminate();
      }
 }
 


/**

The destructor of the ErrorHandler has a special role in that it dumps
all the recorded messages in sequential order into the text file
'ErrorHandler.log' for examination of the conditions which caused a 
fatal error or which occurred until the normal termination of the run. 

@section application Application

Destructors are called automatically when a C++ program exits normally.
If the Terminate() method ends the run it prompts the generation of
the same logfile. 

 */
ErrorHandler::~ErrorHandler()
 {
    ofstream ofs;
    
    // before exiting errors are logged to file
    // ----------------------------------------
    // create output file (overwrite earlier file)
    // log messages
    if ( !error_sequence.empty() )
      {
         ofs.open( "Runtime.log", ios::out|ios::trunc );
         if ( !ofs )
           {
              cout << "ErrorHandler (Destructor): Cannot open error documentation file !\n" << endl;
              return;
           }
         ofs << "ErrorHandler (Destructor): logfile: Runtime messages in chronological sequence " << endl;
         ofs <<"(if this starts with an info or a warning, there were no errors):"<< endl;
         for ( map<string,string>::const_iterator
               err_it=error_sequence.begin(); err_it!=error_sequence.end(); err_it++ )
           ofs << (*err_it).first <<": "<< (*err_it).second << endl;
          
         if ( !merrors.empty() ) {  
             ofs << "\nErrorHandler (Destructor): logfile: CSMP ERROR messages in sequence of occurrence:" << endl;
             multimap<string,string>::const_iterator  wit;
             for ( wit=merrors.begin(); wit!=merrors.end(); wit++ )
               ofs << (*wit).first <<": "<< (*wit).second << endl;
           }
         if ( !mwarnings.empty() ) {  
             ofs << "\nErrorHandler (Destructor): logfile: CSMP WARNING messages in sequence of occurrence:" << endl;
             multimap<string,string>::const_iterator  wit;
             for ( wit=mwarnings.begin(); wit!=mwarnings.end(); wit++ )
               ofs << (*wit).first <<": "<< (*wit).second << endl;
           }
         if ( !minfos.empty() ) {  
             ofs << "\nErrorHandler (Destructor): logfile: CSMP INFO messages in sequence of occurrence:" << endl;
             multimap<string,string>::const_iterator  wit;
             for ( wit=minfos.begin(); wit!=minfos.end(); wit++ )
               ofs << (*wit).first <<": "<< (*wit).second << endl;
           }
           
         ofs.close();
      }
      
 } // end






/**

Method to suppress runtime error messages. The latter can still be
output selectively by listing the message maps using the Print...()
methods.  */
void ErrorHandler::Verbose( bool verb )
{
    verbose = ( verb ? VERBOSE : SILENT );
}
void ErrorHandler::Verbose( size_t verb )
{
    verbose = verb;
}
size_t ErrorHandler::Verbose( ) const
{
    return verbose;
}

void ErrorHandler::Debug( bool debug_flag )
{
    debug = debug_flag;
}
bool ErrorHandler::Debug( ) const
{
    return debug;
}

bool ErrorHandler::operator>( size_t level )  const
{
    return ( verbose > level );
}
bool ErrorHandler::operator>=( size_t level ) const
{
    return ( verbose >= level );
}
bool ErrorHandler::operator==( size_t level ) const
{
    return ( verbose == level );
}
bool ErrorHandler::operator!=( size_t level ) const
{
    return ( verbose != level );
}

 
/**

Various versions of the notice() method are provided such that CSP objects
can report their messages to the ErrorHandler. These are distinguished
by the number and type of arguments. 

@section arguments Input Arguments

The first argument is an enumeration named csmp_error with four optional
values INFO, WARNING, ERROR, and FATAL_ERROR. The last value will prompt
the error handler to terminate the CSP run. An INFO just refers to 
information that may be helpful for the user. In the case of a WARNING
something looks as if it may have gone wrong but it may have no 
consequences for the run or it may have been done on purpose. An ERROR
has occurred if something went wrong such that it corrupts part of the 
computation. Finally, in the case of a FATAL_ERROR, the error is so
severe that there is no point to continue the simulation.  

The second argument must give the object class and the name of the method
which was being executed when the message was send. This argument helps to
track down where a problem did arise. 

The third and other arguments are used for information strings and/or to
describe the condition under which a problem did arise. Variable values 
can be added here through using the number-to-string operators of the 
string class before the string is supplied as a method argument. 

@section return Return Arguments

The input arguments are used by notice() to compose the message that is
either sent to the 'cout', or 'cout' streams or output in a message box. 

@section implementation Implementation

notice() counts the messages of the different types and stores them in the
respective maps for later output. 

@section application Application

Especially, to direct messages to a target output device and to create a
protocol of a simulation run. 
*/
void ErrorHandler::notice( CSMP_MESSAGE err_type, 
                           const string& source, const string& msg )
 {
    string message(source);
   
   // 1. logging the message
    switch( err_type )
      {
         case INFO:        minfos.insert( make_pair(source,msg) );      
                           message += " Info: ";
           break;
         case WARNING:     mwarnings.insert( make_pair(source,msg) );
                           message += " Warning: ";
                           warnings++;
           break;
         case ERROR:       merrors.insert( make_pair(source,msg) );
                           message += " Error: ";
                           errors++;
           break;
         case EXCEPTION:   merrors.insert( make_pair(source,msg) );
                           message += " Exception: ";
                           errors++;
                           WriteErrorsToFile();
                           throw Exception( err_type, source.c_str(), msg );
           break;
         case FATAL_ERROR: merrors.insert( make_pair(source,msg) );
                           message += " Fatal Error: ";
                           cout <<"Fatal Error: "<< source <<" "<< msg << endl;
                           cout.flush();
                           errors++;
                           throw Exception( err_type, source.c_str(), msg );
           break;
         default: cerr <<"\nErrorHandler::notice: error of unknown type was detected."<< endl;
                  terminate();
      }
   message += msg;
   total++;
   
   // 2. recording the messages in chronological order 
   error_sequence[ timer_.ascii() ] = message;
   
   // reporting the error
   if ( verbose ) {
        if ( err_type == INFO ) cout <<"\n"<< message << endl;
        else cout <<"\n"<< message << endl; 
        cout.flush();
#ifndef NDEBUG 
        if ( err_type < WARNING ) {
             cout <<"\nHit return to continue."<< endl;
             getchar();
             getchar();
          }
#endif
     }
   // terminating the run if too many errors occured
   if ( errors > error_limit ) {
       WriteErrorsToFile();
       throw Exception( err_type, source.c_str(), msg );
   }

 } // end notice



void ErrorHandler::notice( CSMP_MESSAGE err_type, const string& source,
                           const string& message1, const string& message2 )
 {
     string  msg(message1);
     msg +=", ";
     msg += message2;
     notice( err_type, source, msg );
 } // end

void ErrorHandler::notice( CSMP_MESSAGE err_type, const char*  source, const string& message )
 {
     notice( err_type, string(source), message );
 }

void ErrorHandler::notice( CSMP_MESSAGE err_type, const char*  source, const string& message1, const string& message2 )
 {
     notice( err_type, string(source), message1, message2 );
 }


void ErrorHandler::notice( CSMP_MESSAGE err_type, const char* source, const char* message )
 {
     notice( err_type, string(source), string(message) );
     
 } // end



void ErrorHandler::notice( CSMP_MESSAGE err_type, const char* source, 
                           const char* message1, const char* message2 )
 {
     string  msg(message1);
     msg     +=", ";
     if ( message2 != NULL ) msg += message2;
     notice( err_type, string(source), msg );

 } // end



/**

To print messages of a specific type, such as INFO, as a list ordered by
the object and method names from which the messages originated. This is 
useful, for instance, to test whether a new object performed in the 
expected way during a complex simulation. 

Apart from PrintInfos() there are methods which print warnings 
(PrintWarnings()), and for errors. 
*/
void ErrorHandler::PrintInfos() const
 {
    cout <<"\nErrorHandler::PrintInfos:"<< endl;
    Print( minfos );
 }
 
 
void ErrorHandler::PrintWarnings() const
 {
    cout <<"\nErrorHandler::PrintWarnings:"<< endl;
    Print( mwarnings );
 }
 
 
void ErrorHandler::PrintErrors() const
 {
    cout <<"\nErrorHandler::PrintErrors:"<< endl;
    Print( merrors );
 }



/**

The standard CSmP interface used to list the state of the object.
Here, all recorded messages are printed on 'stdout' in chronological
sequence. 
*/
void ErrorHandler::Out() const
 {
    cerr <<"\nErrorHandler::Out:"<< endl;
    cerr <<"errors: "<< errors <<", warnings: "<< warnings << endl;
     
    for ( map<string,string>::const_iterator
          it=error_sequence.begin(); it!=error_sequence.end(); it++ )
      cerr << (*it).first << endl << (*it).second << endl << endl;
 }


void ErrorHandler::Print( const multimap<string,string>& m ) const
 {
     for ( multimap<string,string>::const_iterator
           it=m.begin(); it!=m.end(); it++ )
       cerr << (*it).first << endl << (*it).second << endl << endl;
 }





/**

OS-specific program termination which prompted by an error specification
as FATAL_ERROR. Before program terminates, this method writes all the 
messages which occurred until then to the text file 'ErrorHandler.log'.
This is similar to calling the destructor of the error handler.  

@section messages Messages

Terminate() will report if it cannot open the output file and prompt the
user for a new file name such that the event sequence until the error
occurred can still be recorded. 
*/
void ErrorHandler::WriteErrorsToFile( const char* err_file )
 {
    // before exiting errors are logged to file
    // ----------------------------------------
    // create output file (overwrite earlier file)
    // log messages
    map<string,string>::const_iterator  err_it;
    ofstream   ofs;
    char       ofile[50];
    
    if ( !error_sequence.empty() )
      {
         ofs.open ( err_file, ios::out|ios::trunc );
         if ( !ofs )
           {
              cerr << "ErrorHandler::WriteErrorsToFile: Cannot open error documentation file !" << endl;
              cerr << "\nPlease enter a new valid output file name: ";
              cin >> ofile;
              ofs.open ( ofile, ios::out|ios::trunc );
              return;
           }
         ofs <<"ErrorHandler::Terminate: Logfile created before abnormal program termination."<< endl;
         ofs <<"CSMP messages in chronological sequence:" << endl;
         for ( err_it=error_sequence.begin(); err_it!=error_sequence.end(); err_it++ )
           ofs << (*err_it).first <<": "<< (*err_it).second << endl;
         ofs.close();
      }
    
 } // end WriteErrorsToFile


} // end namespace csmp












 
