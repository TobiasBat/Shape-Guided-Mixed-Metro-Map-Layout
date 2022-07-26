//******************************************************************************
// Grid2.h
//	: header file for 2D grid coordinaes
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Mon Mar 14 20:13:35 2011
//
//******************************************************************************

#ifndef	_Grid2_H
#define _Grid2_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
using namespace std;



//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Defining Classes
//------------------------------------------------------------------------------

class Grid2 {

  protected:

    int			_element[ 2 ];	// x, y coordinates

    virtual void	_init( void );	// initialize all coordinates to zero

  public:

//------------------------------------------------------------------------------
//	Constuructors
//------------------------------------------------------------------------------
    Grid2();				// constructor (default)
    Grid2( const int x, const int y );
					// 2D coordinates as input
    Grid2( const Grid2 & v );		// copy constructor
    virtual ~Grid2() {}		        // destructor

//------------------------------------------------------------------------------
//	Assignment opereators
//------------------------------------------------------------------------------
    Grid2 &		operator = ( const Grid2 & v );
				// assignment
    Grid2 &		operator += ( const Grid2 & v );
				// addition + assignment
    Grid2 &		operator -= ( const Grid2& v );
				// subtraction + assignment
    Grid2 &		operator *= ( const int d );
				// scalar product + assignment
//------------------------------------------------------------------------------
//	Reference to elements
//------------------------------------------------------------------------------
    void		init( void )		{ _init(); }
    void		zero( void )		{ _init(); }
				// initialze all the coordinates to zero
    const int &		operator [] ( int i ) const;
    int &		operator [] ( int i );
				// reference to a specific coordinate
    const int *		element( void ) const	{ return _element; }
				// pointer to an array of coordinates
    int &		p( void )		{ return _element[ 0 ]; }
    int &		q( void )		{ return _element[ 1 ]; }
    const int &		p( void ) const		{ return _element[ 0 ]; }
    const int &		q( void ) const		{ return _element[ 1 ]; }
    const int &		getP( void ) const	{ return _element[ 0 ]; }
    const int &		getQ( void ) const	{ return _element[ 1 ]; }
				// reference to a specific coordinate
    void		set( const int p, const int q );
    void		setP( const int p )	{ _element[ 0 ] = p; }
    void		setQ( const int q )	{ _element[ 1 ] = q; }
				// set the coordinate(s)

//------------------------------------------------------------------------------
//	Friend functions
//------------------------------------------------------------------------------
    friend Grid2	operator - ( const Grid2 & v );
				// sign change
    friend Grid2	operator + ( const Grid2 & a, const Grid2 & b );
				// addition
    friend Grid2	operator - ( const Grid2 & a, const Grid2 & b );
				// subtraction
    friend Grid2	operator * ( const int d, const Grid2 & a );
				// scalar product
    friend int		operator * ( const Grid2 & a, const Grid2 & b );
				// inner product

    friend int		operator == ( const Grid2 & a, const Grid2 & b );
				// equivalence
    friend int		operator != ( const Grid2 & a, const Grid2 & b ) {
	return ( ! ( a == b ) );
    }				// inequivalence
    friend int		operator < ( const Grid2 & a, const Grid2 & b );
				// comparison (less than)
    friend int		operator > ( const Grid2 & a, const Grid2 & b );
				// comparison (more than)
    friend int		operator <= ( const Grid2 & a, const Grid2 & b ) {
	return ( ( a == b ) || ( a < b ) );
    }				// comparison (equal to or less than)
    friend int		operator >= ( const Grid2 & a, const Grid2 & b ) {
	return ( ( a == b ) || ( a > b ) );
    }				// comparison (equal to or more than)


//------------------------------------------------------------------------------
//	I/O functions
//------------------------------------------------------------------------------
    friend ostream &	operator << ( ostream & s, const Grid2 & v );
				// 出力
    friend istream &	operator >> ( istream & s, Grid2 & v );
				// 入力
    virtual const char * className( void ) const { return "Grid2"; }
				// クラス名

};


#endif // _Grid2_H
