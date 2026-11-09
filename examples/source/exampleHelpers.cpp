//
//  exampleHelpers.h
//  openCSMP_github
//
//  Created by Stephan Matthai on 28/10/2026.
//  Copyright © 2026 Stephan Matthai. All rights reserved.
//

#include "exampleHelpers.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

using namespace std;

namespace csmp {

void divider( char c ) {
    cout << string( LINE_WIDTH, c ) << endl;
}

void blank() {
    cout << "\n";
}

void heading( uint32_t number, uint32_t total, const string& title ) {
    blank();
    divider();
    cout << "  Step " << number << " of " << total
         << ":  " << title << "\n";
    divider();
    blank();
    cout.flush();
}

void subheading( const string& title ) {
    blank();
    divider('-');
    cout << "  " << title << "\n";
    divider('-');
    blank();
    cout.flush();
}

void pause( const string& prompt ) {
    blank();
    cout << prompt << endl;
    cin.ignore( numeric_limits<streamsize>::max(), '\n' );
}

void code( const string& snippet ) {
    blank();
    divider('.');
    istringstream ss( snippet );
    string line;
    while ( getline( ss, line ) )
        cout << "    " << line << "\n";
    divider('.');
    blank();
    cout.flush();
}

void output( const string& text ) {
    cout << "  >>  " << text << "\n";
    cout.flush();
}

void explain( const string& text ) {
    istringstream ss( text );
    string line;
    while ( getline( ss, line ) )
        cout << "  " << line << "\n";
    cout.flush();
}

void keypoint( const string& text ) {
    blank();
    divider('*');
    cout << "  KEY POINT:\n";
    istringstream ss( text );
    string line;
    while ( getline( ss, line ) )
        cout << "  " << line << "\n";
    divider('*');
    cout.flush();
}

string dbl_str( double v, int prec ) {
    ostringstream oss;
    oss << fixed << setprecision(prec) << v;
    return oss.str();
}

/// Format a double in scientific notation with controlled precision
std::string sci_str( double v, int prec )
{
    ostringstream oss;
    oss << scientific << setprecision(prec) << v;
    return oss.str();
}


} // anonymous namespace
