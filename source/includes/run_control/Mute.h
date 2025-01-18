#ifndef MUTE_H
#define MUTE_H

#include <iostream>
#include <fstream>

namespace csmp {

  /** @brief Allows redirection of the output streams cout and cerr from the screen to null 
      (output goes nowhere) and file (ascii file)
      
      Class is used to suppress output on systems where the user has no access to the console
      and/or where sceen output would dramatically slow the computations.

      Only a single stream can be handled by an instance, so return values will reflect that.
      
      Application example:
      @code
      Mute mute;
      cout << "\nMute mode...\n";
      mute.Null(cout); // or i.e. mute.File( cout, "Output.txt");

      // ... here goes the main code with supressed output

      mute.Reinstate(cout);
      @endcode

      a note 'just in case': if you want to supress samg output too, call following for SAMG_Settings instance (here 'settings')
      @code
      settings.Set_idmp(-1);
      settings.Set_iout1(-1);
      settings.Set_iout2(-1);
      @endcode

  */
  class Mute {
    public:
      Mute();
      ~Mute();

      /// inactivates output to given stream
      bool Null( std::ostream& );
      
      /// redirects output to text file
      bool File( std::ostream&, const char* text_file );
      
      /// brings back original settings
      bool Reinstate( std::ostream& );
      
      /// releases redirection to file
      bool Redirected() const;
      
    private:
      std::streambuf*       cachebf_;
      std::ofstream*        fileRedirect_;
      std::ofstream* const  nullOut_;
 };


} // csmp


#endif