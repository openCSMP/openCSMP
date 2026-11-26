// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only
//
//  Created by Stephan Matthai on 28/10/2026.
//

#include <iostream>
#include <string>
#include <cstdint>


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

} // csmp
