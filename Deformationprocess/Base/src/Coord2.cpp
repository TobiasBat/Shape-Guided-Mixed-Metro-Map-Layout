//******************************************************************************
// Coord2.cc
//	: program file for 2D coordinatse
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Sun Sep 16 15:02:45 2012
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cctype>
#include <cmath>
#include <algorithm>
using namespace std;

#include "Common.h"
#include "Coord2.h"


//------------------------------------------------------------------------------
//	Macro Definitions
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//	Protected Functions
//------------------------------------------------------------------------------

//
//  Coord2::_init --	$B$9$Y$F$NMWAG$r%<%m$K$9$k(B
//
//  $B0z?t(B
//	$B$J$7(B
//
//  $BJV$jCM(B
//	$B$J$7(B
//
void Coord2::_init( void )
{
    _element[ 0 ] = _element[ 1 ] = 0.0;
}

//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//	Constuructors
//------------------------------------------------------------------------------

//
//  Coord2::Coord2 --	$B%3%s%9%H%i%/%?(B($B%G%U%)%k%H(B)
//
//  $B0z?t(B
//	$B$J$7(B
//
//  $BJV$jCM(B
//	$B$J$7(B
//
Coord2::Coord2()
{
    _init();
}


//
//  Coord2::Coord2 --	$B:BI8$rF~NO(B
//
//  $B0z?t(B
//	x, y :	2$B<!85:BI8(B
//
//  $BJV$jCM(B
//	$B$J$7(B
//
Coord2::Coord2( const double x, const double y )
{
    _element[ 0 ]	= x;
    _element[ 1 ]	= y;
}

//
//  Coord2::Coord2 --	$B%3%T!<!&%3%s%9%H%i%/%?(B
//
//  $B0z?t(B
//	v	: 2$B<!85%Y%/%H%k!&%*%V%8%'%/%H(B
//
//  $BJV$jCM(B
//	$B$J$7(B
//
Coord2::Coord2( const Coord2 & v )
{
    _element[ 0 ]	= v._element[ 0 ];
    _element[ 1 ]	= v._element[ 1 ];
 }


//------------------------------------------------------------------------------
//	Assignment opereators
//------------------------------------------------------------------------------

//
//  Coord2::operator = --	assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Coord2 & Coord2::operator = ( const Coord2 & v )
{
    if ( this != &v ) {
	_element[ 0 ]	= v._element[ 0 ];
	_element[ 1 ]	= v._element[ 1 ];
    } 
    return *this;
}


//
//  Coord2::operator += --	addition + assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Coord2 & Coord2::operator += ( const Coord2 & v )
{
    _element[ 0 ]	+= v._element[ 0 ];
    _element[ 1 ]	+= v._element[ 1 ];
    return *this;
}


//
//  Coord2::operator -= --	subtraction + assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Coord2 & Coord2::operator -= ( const Coord2 & v )
{
    _element[ 0 ]	-= v._element[ 0 ];
    _element[ 1 ]	-= v._element[ 1 ];
    return *this;
}


//
//  Coord2::operator -= --	scalar product + assignment
//
//  Inputs
//	d	: scalar value
//
//  Outputs
//	reference to this object
//
Coord2 & Coord2::operator *= ( const double d )
{
    _element[ 0 ]	*= d;
    _element[ 1 ]	*= d;
    return *this;
}


//
//  Coord2::operator /= --	scalar division$B>&(B + assignment
//
//  Inputs
//	d	: scalar value
//
//  Outputs
//	reference to this object
//
Coord2 & Coord2::operator /= ( const double d )
{
    double d_inv = 1./d;
    _element[ 0 ]		*= d_inv;
    _element[ 1 ]		*= d_inv;
    return *this;
}


//
//  Coord2::operator [] --	reference to an element
//
//  Inputs
//	i	: index of the coordinate
//
//  Outputs
//	the corresponding coordinate
//
const double & Coord2::operator [] ( int i ) const
{
#ifdef VEC2_INDEX_CHECK
    const char theName[] = "Coord2::operator [] : ";
    if ( ( i < 0 ) || ( i > 1 ) ) {
	cerr << theName << " index = " << i << endl;
	assert( ( 0 <= i ) && ( i <= 1 ) );
    }
#endif	// VEC2_INDEX_CHECK
    return _element[ i ];
}


//
//  Coord2::operator [] --	reference to an element
//
//  Inputs
//	i	: index of the coordinate
//
//  Outputs
//	the corresponding coordinate
//
double & Coord2::operator [] ( int i )
{
#ifdef VEC2_INDEX_CHECK
    const char theName[] = "Coord2::operator [] : ";
    if ( ( i < 0 ) || ( i > 1 ) ) {
	cerr << theName << " index = " << i << endl;
	assert( ( 0 <= i ) && ( i <= 1 ) );
    }
#endif	// VEC2_INDEX_CHECK
    return _element[ i ];
}


//
//  Coord2::set --	set all the coordinates
//
//  Inputs
//	x, y	: x and y coordinates
//
//  Returns
//	none
//
void Coord2::set( const double x, const double y )
{
    _element[ 0 ]	= x;
    _element[ 1 ]	= y;
}


//------------------------------------------------------------------------------
//	$B%Y%/%H%kFCM-$N4X?t(B
//	Special functions
//------------------------------------------------------------------------------

//
//  Coord2::norm --	return the norm of the vector
//
//  Inputs
//	none
//
//  Outputs
//	norm of this vector
//
double Coord2::norm( void ) const
{
    return sqrt( squaredNorm() );
}


//
//  Coord2::squaredNorm --	return the squared norm of the vector
//
//  Inputs
//	none
//
//  Outputs
//	squared norm of this vector
//
double Coord2::squaredNorm( void ) const
{
    return ( _element[ 0 ]*_element[ 0 ] + _element[ 1 ]*_element[ 1 ] );
}


//
//  Coord2::normalize --	normalize the vector
//
//  Inputs
//	none
//
//  Returns
//	reference to this object (after normalization)
//
Coord2 & Coord2::normalize( void )
    // division-by-zero should be cared by the caller side
{
    double l = norm();
    *this /= l;
    return *this;
}


//
//  Coord2::unit --	compute the unit vector
//
//  Inputs
//	none
//
//  Outputs
//	return the unit vector
//
Coord2 Coord2::unit( void ) const
    // division-by-zero should be cared by the caller side
{
    Coord2 ret( *this );
    ret.normalize();
    return ret;
}




//
//  Coord2::manhattan --	return the Manhattan norm of the vector
//
//  Inputs
//	none
//
//  Outputs
//	Manhattan norm
//
double Coord2::manhattan( void ) const
{
    return ( fabs( _element[ 0 ] ) + fabs( _element[ 1 ] ) );
}

//------------------------------------------------------------------------------
//	Intersection check
//------------------------------------------------------------------------------

//
//  Coord2::cross --	cross product
//
//  Inputs
//	a, b	: two input vectors
//
//  Outputs
//	cross product
//
double crossProd( const Coord2 & a, const Coord2 & b )
{
    return ( a.x() * b.y() - b.x() * a.y() );
}


//
//  Coord2::doubleArea --	return the twice of the triangle area
//
//  Inputs
//	a, b, c :	coordinate of the three corner points
//
//  Outputs
//	twice of the area
//
double doubleArea( const Coord2 & a, const Coord2 & b, const Coord2 & c )
{
    return crossProd( ( a - b ), ( a - c ) );
}


//
//  isSeparate --	return true if the bounding boxes spanned by "a-b" and
//			"c-d" has no overlap.
//
//  Inputs
//	a, b, c :	coordinate of the three points
//
//  Outputs
//	boolean value
//
bool isSeparate( const Coord2 & a, const Coord2 & b, const Coord2 & c, const Coord2 & d )
{
    double xminAB = min( a.x(), b.x() );
    double xmaxAB = max( a.x(), b.x() );
    double yminAB = min( a.y(), b.y() );
    double ymaxAB = max( a.y(), b.y() );
    double xminCD = min( c.x(), d.x() );
    double xmaxCD = max( c.x(), d.x() );
    double yminCD = min( c.y(), d.y() );
    double ymaxCD = max( c.y(), d.y() );

    if ( xmaxCD < xminAB ) return true;
    if ( xmaxAB < xminCD ) return true;
    if ( ymaxCD < yminAB ) return true;
    if ( ymaxAB < yminCD ) return true;

    return false;
}


bool isCollinear ( const Coord2 & a, const Coord2 & b, const Coord2 & c )
{
    return ( fabs( doubleArea( a, b, c ) ) < 1.0e-8 );
}

bool isLeft( const Coord2 & a, const Coord2 & b, const Coord2 & c )
{
    return ( doubleArea( a, b, c ) > 1.0e-8 );
}

bool isCCW( const Coord2 & a, const Coord2 & b, const Coord2 & c )
{
    return isLeft( a, b, c );
}

bool isLeftOn( const Coord2 & a, const Coord2 & b, const Coord2 & c )
{
    return ( doubleArea( a, b, c ) > -USER_EPS );
}


// 
//
//  Coord2::isIntersected -	Find intersection between the segments ab and cd. 
//				If they are intersected, p will be the
//				interecting point and s and t are internal
//				ratios for the two segments, respectively. 
//
//  Inputs
//	none
//
//  Outputs
//	boolean value acoording to the intersection between two line segments
//
bool isIntersected( const Coord2 & a, const Coord2 & b, const Coord2 & c, const Coord2 & d ) 
{
    if ( isSeparate( a, b, c, d ) ) return false;

    if ( isCollinear( a, b, c ) ||
	 isCollinear( a, b, d ) ||
	 isCollinear( c, d, a ) ||
	 isCollinear( c, d, a ) ) return false;

    return ( ( doubleArea( a, b, c ) * doubleArea( a, b, d ) < 0.0 ) &&
	     ( doubleArea( c, d, a ) * doubleArea( c, d, b ) < 0.0 ) ); 
}


// 
//
//  Coord2::isIntersected -	Find intersection between the segments ab and cd. 
//				If they are intersected, p will be the
//				interecting point and s and t are internal
//				ratios for the two segments, respectively. 
//
//  Inputs
//	intersection	: intersection point
//
//  Outputs
//	boolean value acoording to the intersection between two line segments
//
bool isIntersected( const Coord2 & a, const Coord2 & b, const Coord2 & c, const Coord2 & d,
		    Coord2 & intersection ) 
{
#ifdef DEBUG
    cerr << " a = " << a.x << " , " << a.y << endl;
    cerr << " b = " << b.x << " , " << b.y << endl;
    cerr << " c = " << c.x << " , " << c.y << endl;
    cerr << " d = " << d.x << " , " << d.y << endl;
#endif	// DEBUG

    double denominator =
        ( double )a.x() * ( double )( d.y() - c.y() ) +
        ( double )b.x() * ( double )( c.y() - d.y() ) +
        ( double )d.x() * ( double )( b.y() - a.y() ) +
        ( double )c.x() * ( double )( a.y() - b.y() );

    // If denominator vanishes, the two segments are parallel to each other
    // In that case, return false
    if ( fabs( denominator ) < USER_EPS ) return false;

    double s =  (
		 ( double )a.x() * ( double )( d.y() - c.y() ) +
		 ( double )c.x() * ( double )( a.y() - d.y() ) +
		 ( double )d.x() * ( double )( c.y() - a.y() )
		 ) / denominator;
    double t = -( 
		 ( double )a.x() * ( double )( c.y() - b.y() ) +
		 ( double )b.x() * ( double )( a.y() - c.y() ) +
		 ( double )c.x() * ( double )( b.y() - a.y() )
		  ) / denominator;

#ifdef DEBUG
    cerr << " deno = " << denominator << " s = " << s << " t = " << t << endl;
#endif  // DEBUG

    // If the intersecting point exists within the segments, return true
    if ( ( -USER_EPS <= s ) && ( s <= 1.0 + USER_EPS ) &&
         ( -USER_EPS <= t ) && ( t <= 1.0 + USER_EPS ) ) {
#ifdef NONEED
        s = MIN2( 1.0, MAX2( 0.0, s ) );
        t = MIN2( 1.0, MAX2( 0.0, t ) );
#endif	// NONEED
#ifdef DEBUG
	cerr << " True" << endl;
#endif	// DEBUG
	intersection.setX( ( double )a.x() + s * ( double )( b.x() - a.x() ) );
	intersection.setY( ( double )a.y() + s * ( double )( b.y() - a.y() ) );
        return true;
    }
    else {
#ifdef DEBUG
	cerr << " False" << endl;
#endif	// DEBUG
        return false;
    }
}


// 
//
//  Coord2::isIntersected -	Find intersection between the segments ab and cd. 
//				If they are intersected, p will be the
//				interecting point and s and t are internal
//				ratios for the two segments, respectively. 
//
//  Inputs
//	none
//
//  Outputs
//	boolean value acoording to the intersection between two line segments
//
bool doConflict( const Coord2 & a, const Coord2 & b, const Coord2 & c, const Coord2 & d ) 
{
    const double gap = 0.01;
    Coord2 cc =   (1.0+gap) * c - gap * d;
    Coord2 dd =  -gap * c + (1.0+gap) * d;

    if ( isSeparate( a, b, cc, dd ) ) return false;

    if ( isCollinear( a, b, cc ) ||
	 isCollinear( a, b, dd ) ||
	 isCollinear( cc, dd, a ) ||
	 isCollinear( cc, dd, a ) ) return false;

    return ( ( doubleArea( a, b, cc ) * doubleArea( a, b, dd ) < 0.0 ) &&
	     ( doubleArea( cc, dd, a ) * doubleArea( cc, dd, b ) < 0.0 ) ); 
}


//------------------------------------------------------------------------------
//	Friend functions
//------------------------------------------------------------------------------

//
//  operator - --	sign change
//
//  Inputs
//	a	: 2D coordinates
//
//  Outputs
//	2D coordinates in the opposite direction
//
Coord2 operator - ( const Coord2 & a )
{
    return Coord2( -a._element[ 0 ], -a._element[ 1 ] );
}


//
//  operator + --	addition
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	addition of the two 2D coordinates
//
Coord2 operator + ( const Coord2 & a, const Coord2 & b )
{
    return Coord2( a._element[ 0 ] + b._element[ 0 ],
		   a._element[ 1 ] + b._element[ 1 ] );
}


//
//  operator - --	difference
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	difference of the two 2D coordinates
//
Coord2 operator - ( const Coord2 & a, const Coord2 & b )
{
    return Coord2( a._element[ 0 ] - b._element[ 0 ],
		   a._element[ 1 ] - b._element[ 1 ] );
}


//
//  operator * --	scalar product
//
//  Inputs
//	d	: scalar value
//	a	: 2D coordinates
//
//  Outputs
//	scalar product 
//
Coord2 operator * ( const double d, const Coord2 & a )
{
    return Coord2( d * a._element[ 0 ], d * a._element[ 1 ] );
}


//
//  operator * --	inner product
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	inner product
//
double operator * ( const Coord2 & a, const Coord2 & b )
{
    return ( a._element[ 0 ] * b._element[ 0 ] +
	     a._element[ 1 ] * b._element[ 1 ] );
}


//
//  operator / --	scalar division
//
//  Inputs
//	a	: 2D coordinates
//	d	: scalar value
//
//  Outputs
//	scalar division
//
Coord2 operator / ( const Coord2 & a, const double d )
{
    double d_inv = 1./d;
    return Coord2( a._element[ 0 ] * d_inv, a._element[ 1 ] * d_inv );
}


//
//  operator == --	equivalence
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	boolean value
//
int operator == ( const Coord2 & a, const Coord2 & b )
{
    return ( ( a._element[ 0 ] == b._element[ 0 ] ) &&
	     ( a._element[ 1 ] == b._element[ 1 ] ) );
}


//
//  operator < --	comparison (less than)
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	boolean value
//
int operator < ( const Coord2 & a, const Coord2 & b )
{
    if ( a._element[ 0 ] < b._element[ 0 ] ) return true;
    else if ( a._element[ 0 ] > b._element[ 0 ] ) return false;
    else {
	if ( a._element[ 1 ] < b._element[ 1 ] ) return true;
	else if ( a._element[ 1 ] > b._element[ 1 ] ) return false;
	else return false;
    }
}


//
//  operator > --	comparison (more than)
//
//  Inputs
//	a, b	: 2D coordinates
//
//  Outputs
//	boolean value
//
int operator > ( const Coord2 & a, const Coord2 & b )
{
    if ( a._element[ 0 ] > b._element[ 0 ] ) return true;
    else if ( a._element[ 0 ] < b._element[ 0 ] ) return false;
    else {
	if ( a._element[ 1 ] > b._element[ 1 ] ) return true;
	else if ( a._element[ 1 ] < b._element[ 1 ] ) return false;
	else return false;
    }
}


//------------------------------------------------------------------------------
//	I/O functions
//------------------------------------------------------------------------------

//
//  operator << --	output
//
//  Inputs
//	s	: reference to output stream
//	v	: 2D coordinates
//
//  Outputs
//	reference to output stream
//
ostream & operator << ( ostream & stream, const Coord2 & obj )
{
    int i;		// loop counter
    // set the output formatting
    //stream << setiosflags( ios::showpoint );
    //stream << setprecision( 8 );
    //int width = 16;
    // print out the elements
    for ( i = 0; i < 2; i++ ) {
	//stream << setw( width ) << obj._element[ i ];
	stream << setw( 4 ) << obj._element[ i ];
	if ( i != 1 ) stream << "\t";
    }
    stream << endl;

    return stream;
}


//
//  operator >> --	input
//
//  Inputs
//	s	: reference to input stream
//	v	: 2D coordinates
//
//  Outputs
//	reference to input stream
//
istream & operator >> ( istream & stream, Coord2 & obj )
{
    int i;		// loop counter
    // reading the elements
    for ( i = 0; i < 2; i++ )
	stream >> obj._element[ i ];
    return stream;
}




