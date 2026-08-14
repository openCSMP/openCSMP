/// This is a demonstration of how to use NewtonRaphson.h<br />
/// to compile:  g++ -I. -lm NewtonRaphson.cpp<br />
/// description: this test the algorithm in 3 situations:<br />
///    1) calculating df at bouth points, to find the average df, for each iteration, with no boundaries <br />
///    2) calculating d2f to get second order approximation at each iteration, with no boundaries<br />
///    3) calculating d2f to get second order approximation at each iteration, and a boundary solution<br />

/***************************************************************************
 *   Copyright (C) 2009 by Clark Sims                                      *
 *   http://AcumenSoftwareInc.com/WhoWeAre/Clark_Sims.html                 *
 *   ClarkSims@AcumenSoftwareInc.com                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "NewtonRaphson.h"
#include <cmath>

/// \brief This functor is provided as a demonstration of how to create a functor for the template class NetwonRaphson0. It encapsulates the function: exp(x) - offset.  The function variable is x. Offset is a hidden variable. This class is defined in the demonstration file NewtonRaphson.cpp.

/*! NewtonRaphson

<p>
This library contains a template function NewtonRaphsonSolve0 which is an implentation of the Newton Raphson algorithm. The template uses two type, functor, and real. "real" can be any implentation of one dimensional floating point arithmetic. "functor" is a class, wich must have a member function pointer, "_f", which will be the functions which is solveved. There are several other arguments, which can be used to tweek the algorithm. These arguments are described in the documentation for the constructor of the template class NewtonRaphsonSolve0::NewtonRaphsonSolve0. The actual iteration is performed by the member function NewtonRaphsonSolve0::do_iteration.
</p>

<p>
One can view or download the original source code for NewtonRaphson.h <a href="NewtonRaphson.h"> here </a>. This header file contains the code which defines the template class NewtonRaphsonSolve0. <br/>
One can view or download the original source code for NewtonRaphson.cpp <a href="NewtonRaphson.cpp"> here </a>. This source file contains three examples which call call objects of class NewtonRaphsonSolve0.
</p>

<p>
This computer code is being released under the GNU general public license:
</p>

<pre>
Copyright (C) 2009 by Clark Sims                                      
http://AcumenSoftwareInc.com/WhoWeAre/Clark_Sims.html                 
ClarkSims@AcumenSoftwareInc.com                                       
                                                                          
This program is free software; you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by  
the Free Software Foundation; either version 2 of the License, or     
(at your option) any later version.                                   
                                                                          
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU General Public License for more details.                          
                                                                          
You should have received a copy of the GNU General Public License     
along with this program; if not, write to the                         
Free Software Foundation, Inc.,                                       
59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             
</pre>

*/
