//******************************************************************************
// Grid2.cc
//	: program file for 2D grid coordinatse
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Mon Mar 14 02:16:23 2011
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cctype>
#include <cmath>
using namespace std;

#include "Grid2.h"


//------------------------------------------------------------------------------
//	Macro Definitions
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//	Protected Functions
//------------------------------------------------------------------------------

//
//  Grid2::_init --	すべての要素をゼロにする
//
//  引数
//	なし
//
//  返り値
//	なし
//
void Grid2::_init( void )
{
    _element[ 0 ] = _element[ 1 ] = 0;
}

//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//	Constuructors
//------------------------------------------------------------------------------

//
//  Grid2::Grid2 --	コンストラクタ(デフォルト)
//
//  引数
//	なし
//
//  返り値
//	なし
//
Grid2::Grid2()
{
    _init();
}


//
//  Grid2::Grid2 --	座標を入力
//
//  引数
//	x, y :	2次元座標
//
//  返り値
//	なし
//
Grid2::Grid2( const int x, const int y )
{
    _element[ 0 ]	= x;
    _element[ 1 ]	= y;
}

//
//  Grid2::Grid2 --	コピー・コンストラクタ
//
//  引数
//	v	: 2次元ベクトル・オブジェクト
//
//  返り値
//	なし
//
Grid2::Grid2( const Grid2 & v )
{
    _element[ 0 ]	= v._element[ 0 ];
    _element[ 1 ]	= v._element[ 1 ];
 }


//------------------------------------------------------------------------------
//	Assignment opereators
//------------------------------------------------------------------------------

//
//  Grid2::operator = --	assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Grid2 & Grid2::operator = ( const Grid2 & v )
{
    if ( this != &v ) {
	_element[ 0 ]	= v._element[ 0 ];
	_element[ 1 ]	= v._element[ 1 ];
    }
    return *this;
}


//
//  Grid2::operator += --	addition + assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Grid2 & Grid2::operator += ( const Grid2 & v )
{
    _element[ 0 ]	+= v._element[ 0 ];
    _element[ 1 ]	+= v._element[ 1 ];
    return *this;
}


//
//  Grid2::operator -= --	subtraction + assignment
//
//  Inputs
//	v	: 2D coordinates
//
//  Outputs
//	reference to this object
//
Grid2 & Grid2::operator -= ( const Grid2 & v )
{
    _element[ 0 ]	-= v._element[ 0 ];
    _element[ 1 ]	-= v._element[ 1 ];
    return *this;
}


//
//  Grid2::operator -= --	scalar product + assignment
//
//  Inputs
//	d	: scalar value
//
//  Outputs
//	reference to this object
//
Grid2 & Grid2::operator *= ( const int d )
{
    _element[ 0 ]	*= d;
    _element[ 1 ]	*= d;
    return *this;
}


//
//  Grid2::operator [] --	reference to an element
//
//  Inputs
//	i	: index of the coordinate
//
//  Outputs
//	the corresponding coordinate
//
const int & Grid2::operator [] ( int i ) const
{
#ifdef GRID2_INDEX_CHECK
    const char theName[] = "Grid2::operator [] : ";
    if ( ( i < 0 ) || ( i > 1 ) ) {
	cerr << theName << " index = " << i << endl;
	assert( ( 0 <= i ) && ( i <= 1 ) );
    }
#endif	// GRID2_INDEX_CHECK
    return _element[ i ];
}


//
//  Grid2::operator [] --	reference to an element
//
//  Inputs
//	i	: index of the coordinate
//
//  Outputs
//	the corresponding coordinate
//
int & Grid2::operator [] ( int i )
{
#ifdef GRID2_INDEX_CHECK
    const char theName[] = "Grid2::operator [] : ";
    if ( ( i < 0 ) || ( i > 1 ) ) {
	cerr << theName << " index = " << i << endl;
	assert( ( 0 <= i ) && ( i <= 1 ) );
    }
#endif	// GRID2_INDEX_CHECK
    return _element[ i ];
}


//
//  Grid2::set --	set all the coordinates
//
//  Inputs
//	x, y	: x and y coordinates
//
//  Returns
//	none
//
void Grid2::set( const int x, const int y )
{
    _element[ 0 ]	= x;
    _element[ 1 ]	= y;
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
Grid2 operator - ( const Grid2 & a )
{
    return Grid2( -a._element[ 0 ], -a._element[ 1 ] );
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
Grid2 operator + ( const Grid2 & a, const Grid2 & b )
{
    return Grid2( a._element[ 0 ] + b._element[ 0 ],
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
Grid2 operator - ( const Grid2 & a, const Grid2 & b )
{
    return Grid2( a._element[ 0 ] - b._element[ 0 ],
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
Grid2 operator * ( const int d, const Grid2 & a )
{
    return Grid2( d * a._element[ 0 ], d * a._element[ 1 ] );
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
int operator * ( const Grid2 & a, const Grid2 & b )
{
    return ( a._element[ 0 ] * b._element[ 0 ] +
	     a._element[ 1 ] * b._element[ 1 ] );
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
int operator == ( const Grid2 & a, const Grid2 & b )
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
int operator < ( const Grid2 & a, const Grid2 & b )
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
int operator > ( const Grid2 & a, const Grid2 & b )
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
ostream & operator << ( ostream & stream, const Grid2 & obj )
{
    // set the output formatting
    int width = 16;
    // print out the elements
    for ( int i = 0; i < 2; i++ ) {
	stream << setw( width ) << obj._element[ i ];
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
istream & operator >> ( istream & stream, Grid2 & obj )
{
    // reading the elements
    for ( int i = 0; i < 2; i++ )
	stream >> obj._element[ i ];
    return stream;
}




