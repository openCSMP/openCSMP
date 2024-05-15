#include "BE_Time.h"

using namespace std;

// exercise: tokenizing timestring
void BE_Time::Time( long& yr, long& hr, long& mi, long& sec )
 {
    // get the current time
    updateAscii();
 
    // get tools to break timestring up
    char *token;
    const char* const delims =" :";

    // Tokenizing the time-string: Fri Apr 23 09:19:43 1999
    // ------------------------------------------------------
    // day
    token = strtok( Ascii, delims );
    // month
    token = strtok(NULL, delims);
    // date
    token = strtok(NULL, delims);
    // hour
    token = strtok(NULL, delims);
    hr    = atol(token);
    // minute
    token = strtok(NULL, delims);
    mi    = atol(token);
    // second
    token = strtok(NULL, delims);
    sec   = atol(token);
    // get the year 
    token = strtok(NULL, delims );
    yr    = atol(token);
     
    // get again the current time string
    updateAscii();

 } // end time


// exercise: tokenizing timestring
// tested: O. K. SKM
void BE_Time::Time( long& year, long& month, long& day )
 {
    // get the current time
    updateAscii();
 
    // get tools to break timestring up
    char *token;
    const char* const delims =" :";

    // Tokenizing the time-string: Fri Apr 23 09:19:43 1999
    // ------------------------------------------------------
    // day literal
    token = strtok( Ascii, delims );
    // month literal
    token = strtok(NULL, delims);
    if      ( strcmp(token,"Jan") == 0 ) month = 1;
    else if ( strcmp(token,"Feb") == 0 ) month = 2;
    else if ( strcmp(token,"Mar") == 0 ) month = 3;
    else if ( strcmp(token,"Apr") == 0 ) month = 4;
    else if ( strcmp(token,"May") == 0 ) month = 5;
    else if ( strcmp(token,"Jun") == 0 ) month = 6;
    else if ( strcmp(token,"Jul") == 0 ) month = 7;
    else if ( strcmp(token,"Aug") == 0 ) month = 8;
    else if ( strcmp(token,"Sep") == 0 ) month = 9;
    else if ( strcmp(token,"Oct") == 0 ) month = 10;
    else if ( strcmp(token,"Nov") == 0 ) month = 11;
    else if ( strcmp(token,"Dec") == 0 ) month = 12;
    else
    cout <<"\nBE_Time::Time: Could not parse month: "<< token << endl;
    // day number
    token = strtok(NULL, delims);
    day   = atoi(token);
    // hour
    token = strtok(NULL, delims);
    // minute
    token = strtok(NULL, delims);
    // second
    token = strtok(NULL, delims);
    // year 
    token = strtok(NULL, delims );
    year  = atol(token);
    // get again the current time string
    updateAscii();

 } // end time
