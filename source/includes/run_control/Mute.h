#ifndef MUTE_H
#define MUTE_H

#include <iostream>
#include <fstream>

namespace csmp{


  /// Allows redirection of ostream (i.e. cout, cout...) from screen to null (output goes nowhere) and file (ascii file)
  class Mute
  {
    std::streambuf*       cachebf_;
    std::ofstream*        fileRedirect_;
    std::ofstream* const  nullOut_;

  public:
    Mute();
    ~Mute();

    bool Null( std::ostream& );
    bool File( std::ostream&, const char* );
    bool Reinstate( std::ostream& );
    bool Redirected() const;
  };


  Mute::Mute()
    : cachebf_(nullptr), fileRedirect_(nullptr), nullOut_( new std::ofstream("/dev/null") )
  {

  }


  /// Dtor does not reinstate ostream
  Mute::~Mute()
  {

  }


  /// Redirects to null stream, returns false if stream already redirected, true if redirection successful
  bool Mute::Null( std::ostream& os )
  {
    if(cachebf_)
      return false;

    cachebf_ = os.rdbuf();
    os.rdbuf( nullOut_->rdbuf() );
    return true;
  }


  /// Redirects to file, returns false if stream already redirected, true if redirection successful
  bool Mute::File( std::ostream& os, const char* fileName )
  {
    if( cachebf_ || fileRedirect_ )
      return false;

    fileRedirect_ = new ofstream;
    fileRedirect_->open(fileName);
    if( !fileRedirect_->is_open() )
      return false;

    cachebf_ = os.rdbuf();
    os.rdbuf( fileRedirect_->rdbuf() );
    return true;
  }


  /// Reinstates screen output, returns false if no redirection current, true if reinstation successful
  bool Mute::Reinstate( std::ostream& os )
  {
    if(!cachebf_)
      return false;

    if(fileRedirect_)
      fileRedirect_->close();

    os.rdbuf(cachebf_);
    cachebf_ = nullptr;
    fileRedirect_ = nullptr;
    return true;
  }


  bool Mute::Redirected() const
  {
    return static_cast<bool>(cachebf_);
  }


  /** 
  @class Mute
  @brief Allows redirection of ostream (i.e. cout, cout...) from screen to null (output goes nowhere) and file (ascii file)
  @author P Lang

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


} // csmp


#endif