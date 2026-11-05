//
//  exampleHelpers.h
//  openCSMP_github
//
//  Created by Stephan Matthai on 28/10/2026.
//  Copyright © 2026 Stephan Matthai. All rights reserved.
//

#include <iostream>
#include <string>


namespace csmp {

constexpr int LINE_WIDTH{ 70 };

void divider( char c = '=' );

void blank();

void heading( uint32_t number, uint32_t total, const std::string& title );

void subheading( const std::string& title );

void pause( const std::string& prompt = "  [ Press ENTER to continue ]" );

void code( const std::string& snippet );

void output( const std::string& text );

void explain( const std::string& text );

void keypoint( const std::string& text );

std::string dbl_str( double v, int prec = 4 );

std::string sci_str( double v, int prec = 3 );

} // anonymous namespace
