// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "AP_algebraUtilities.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {


/** Area center of mass of a planar, closed, convex quadrilateral defined by vecPoints.

@section arguments Input Arguments 

The four vertices defining the quadrilateral. 

@param vecCenter The center of area.

@section implementation Implementation

Algorithm: 
-suppose the quadrilateral has vertices: p0,p1,p2,p3 
-suppose center of mass (p0,p1,p2)=c1; center of mass (p0,p2,p3)=c2 
-suppose center of mass (p3,p0,p1)=c3; center of mass (p1,p3,p2)=c4 
The area centroid of mass will be the intersection between c1c2 and c3c4; 

@section application Application

The main idea is to calculate the center of the area of the given quadrilateral.
*/
template<uint32_t dim>
void areaCenterOfMass(const vector<Point<dim> >& vecPoints, Point<dim>& vecCenter)
{
	//if vector<double> of points is empty, return an empty centroid
	assert(vecPoints.size() == 4);
	
	Point<dim> pt1, pt2, pt3, pt4;
	vertexCenterOfMass3Vertices(vecPoints[0],vecPoints[1],vecPoints[3], pt1);
	vertexCenterOfMass3Vertices(vecPoints[1],vecPoints[2],vecPoints[3], pt2);
	vertexCenterOfMass3Vertices(vecPoints[0],vecPoints[1],vecPoints[2], pt3);
	vertexCenterOfMass3Vertices(vecPoints[0],vecPoints[2],vecPoints[3], pt4);
	
	//pt1 and pt2 delimit line 1
	//pt3 and pt4 delimit line 2
	intersection(pt1, pt2, pt3, pt4, vecCenter);
}

template void areaCenterOfMass( const vector<Point<3U> >&, Point<3U>& );







/** Calculates the intersection of two lines in 3D.

@section arguments Input Arguments 

Four delimiting point of the intersecting lines. Pt1 and pt2 delimit the first 
line while pt3 and pt4 delimit the second line.

@return true if the lines do intersect, false if they are parallel or equal.
If intersecting, it returns the point where they intersect.

@section implementation Implementation

Based on Goldman 1990 - Mathworld Line-Line Intersection 
Finds the inetrsection between the two lines, defined each by
the given pair of 3d coordinates.

@section application Application

Determine the point where two lines intersect.
*/
template<uint32_t dim>
bool intersection(const Point<dim>& c1, const Point<dim>& c2, const Point<dim>& c3, const Point<dim>& c4, Point<dim>& vecIntersection)
{
	//assume that c1,c2,c3, and c4 are coplanar.
	 
	//we define a = c2-c1
	//          b = c4-c3
	//          c = c3-c1
  Point<dim> a(c2-c1), 
                b(c4-c3), 
                c(c3-c1),
                axb(crossProduct(a,b)),
                cxb(crossProduct(c,b));
  
	//the intersection point p is: p= p1 + a*( (cxb).(axb)/ |axb|^2 )
	
	vecIntersection = c1 + (a* (dotProduct(cxb,axb)/axb.SquaredLength()) );
	
	return true;
}




/// dimensionless cross product
vector<double> crossProduct( const vector<double>& vector1, const vector<double>& vector2 )
{
	const size_t iSize(vector1.size());

	assert(iSize == vector2.size());
	assert(iSize == 2 || iSize == 3);
	
	vector<double> vecReturn;
	
	if(iSize == 2)
	{
		vecReturn.reserve(1U);
		vecReturn.push_back(vector1[0]*vector2[1]-vector1[1]*vector2[0]);
	}
	else if(iSize == 3)
	{
		vecReturn.reserve(3U);
		vecReturn.push_back(vector1[1]*vector2[2]-vector1[2]*vector2[1]);
		vecReturn.push_back(vector1[2]*vector2[0]-vector1[0]*vector2[2]);
		vecReturn.push_back(vector1[0]*vector2[1]-vector1[1]*vector2[0]);
	}
	
	return vecReturn;
}





/// 3x3 cross product
vector<double> crossProduct3by3(const double& f1_0, const double& f1_1, const double& f1_2, const double& f2_0, const double& f2_1, const double& f2_2)
{
	vector<double> vecReturn(3U);
		
	vecReturn[0] = f1_1*f2_2-f1_2*f2_1;
	vecReturn[1] = f1_2*f2_0-f1_0*f2_2;
	vecReturn[2] = f1_0*f2_1-f1_1*f2_0;
	
	return vecReturn;
}


bool normalOfPolygon(const vector< Point<3U> >& vecPolygon, const Point<3U>& vecNormalAt, Point<3U>& vecNormal)
{
  //assumes that the polygon is planar and convex
	
	//if empty, return false
	if(vecPolygon.empty())
		return false;
	
	//get the dimension of the points, and the number of points	
	const size_t iNrPts(vecPolygon.size());
	
	if ( iNrPts == 1 )
	{
		vecNormal[0] = 1.;
	}			
	else if ( iNrPts == 2 )
	{
		//if the polygon has 2 points only -we call it line
		//we calculate the orthogonal vector<double> to the line and to the 0,0,1 vector<double> 
		vecNormal = crossProduct(vecPolygon[1]-vecPolygon[0],Point<3U>(0,0,1));
		vecNormal.NormalizeLengthTo(1.);
	}
	else if ( iNrPts >= 3 )
	{
		localSurfaceNormal(vecPolygon, vecNormalAt, 2, vecNormal);
	}
	
	else
		return false;
		
	return true;	
}

/**

normal approximation (non-planar surface)
using the "local surface normal" per approximation
*/

template<uint32_t dim>
void localSurfaceNormal(const vector< Point<dim> >& vecPoints, const Point<dim>& vecNormalAt, const size_t& iLevelOfRefinement, Point<dim>& vecNormal)
{
	// vecPolygon is of the form -> p1, p2, p3, p4, ... 
	assert(!vecPoints.empty());
		
	//vecNormalAt is where the normal is to be calculated
	
	//refine
	vector< Point<dim> > vecRefPoints; // := p1, p1/2, p2, p3/2, p3, ...
	
	if (iLevelOfRefinement <= 1.)
	  {
	    vecRefPoints = vecPoints;
	  }
	else
	  {
	   Point<dim> vecTemp;
	   
	    //refine!
	   //for vecPoints on the surface, generate pi/2, such that for each pi-pi+1 pairs of vertices pi/2 = //(pi + pi+1)/2
	   const size_t iSize(vecPoints.size());
	   for(size_t iPt = 0U; iPt < iSize; iPt++)
	   {
		   for(auto i{0U}; i < iLevelOfRefinement; i++)	
		   {
			  const size_t iOffset((iPt+1)%iSize);
			  
			  //transform(vecPoints[iOffset].begin(), vecPoints[iOffset].end(), vecPoints[iPt].begin(), vecTemp.begin(), minus());
			  vecTemp = vecPoints[iOffset]-vecPoints[iPt];
			  //transform(vecTemp.begin(), vecTemp.end(), vecTemp.begin(), bind2nd(multiplies(),static_cast<double>(i)/static_cast<double>(iLevelOfRefinement)));
			  vecTemp *= static_cast<double>(i)/static_cast<double>(iLevelOfRefinement);
			  //transform(vecTemp.begin(), vecTemp.end(), vecPoints[iPt].begin(), vecTemp.begin(), plus());
			  vecTemp += vecPoints[iPt];
			  
			  vecRefPoints.push_back(vecTemp);
		   }
	   } 
	  }
		
	Point<dim> vecTemp1, vecTemp2;
		
	const size_t iSizeRef(vecRefPoints.size());
	for(size_t iPt = 0U; iPt < iSizeRef; iPt++)
   	{
   		vecTemp1 = vecRefPoints[iPt]-vecNormalAt;
   		vecTemp2 = vecRefPoints[(iPt+1)%iSizeRef]-vecNormalAt;  
   		//transform(vecRefPoints[iPt].begin(), vecRefPoints[iPt].end(), vecNormalAt.begin(), vecTemp1.begin(), minus());
   		//transform(vecRefPoints[(iPt+1)%iSizeRef].begin(), vecRefPoints[(iPt+1)%iSizeRef].end(), vecNormalAt.begin(), vecTemp2.begin(), minus());
   		vecNormal += crossProduct(vecTemp1, vecTemp2);
   	}	   	
	
	
	vecNormal.NormalizeLengthTo(1.);
}

/// Reference: Ronald Goldman, "Area of Planar Polygons and Volume of Polyhedra" in Graphics Gems II (1994)

template<uint32_t dim>
bool areaOfPolygon(const vector< Point<dim> >& vecPolygon, const size_t& iNrOfFacetPoints, double& fArea)
{
	//if the polygon has only one point, we define the area as 1
	if ( iNrOfFacetPoints < 2 ) {
		fArea = 1.;
		return true;
	  }
	
    if(vecPolygon.empty())
    	return false;
    
    //if the polygon is a line, the area is actually a distance
	if ( iNrOfFacetPoints == 2 ) {
	    fArea = vecPolygon[1].DistanceTo(vecPolygon[0]);
	    return true;
	  }
	
	//else we define area as:
    //2 A(P) = abs(N . (sum_{i=0}^{n-1} (v_i x v_{i+1})))
	
	//get the normal, or approximation of the normal
	Point<dim> vecNormal, vecCentroid;
	vertexCenterOfMass(vecPolygon, vecCentroid);
	
	if(!normalOfPolygon(vecPolygon, vecCentroid, vecNormal))
	{
		//clear contents of area variable
		fArea = 0.;
		return false; //no normal, no area	
  }
    
	//vecSum will hold the sum of the cross products
	
	//calculate summ of all cross products
	const size_t iSize(vecPolygon.size());
	Point<dim> vecSum;
	for(size_t iPt = 0U; iPt < iSize; iPt++)
   	vecSum += crossProduct(vecPolygon[iPt],vecPolygon[(iPt+1)%iSize]);
    
	//calculate area
	fArea = 0.5 * fabs( dotProduct(vecNormal, vecSum));
	
	//successful computation
	return true;	
}


template bool areaOfPolygon(const vector< Point<3U> >&, const size_t&, double& );



/**

Computes the distance between two points (in any dimension).
The dimension of the points is given by the size of the input vectors.
*/
double distanceBetweenPoints(const vector<double>& v1, const vector<double> & v2)
{
	double fDist(0.);
	
	//just calculate the distance
	assert(v1.size() == v2.size());
	
	vector<double>::const_iterator vIter2(v2.begin());
	const vector<double>::const_iterator vIterEnd(v1.end());
	for(vector<double>::const_iterator vIter1 = v1.begin(); vIter1 != vIterEnd; vIter1++, vIter2++)
	{
		fDist += square( *vIter1 - *vIter2 );
	}
			
	fDist = sqrt(fDist);
    
    return fDist;
}

} // csmp

