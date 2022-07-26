//******************************************************************************
// BaseGraphProperty.h
//	: header file for base graph property
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Dec 27 23:16:12 2018
//
//******************************************************************************

#ifndef	_Graph_BaseGraphProperty_H
#define _Graph_BaseGraphProperty_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>

using namespace std;

#include "Coord2.h"
#include "Common.h"


//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace Graph {

    //------------------------------------------------------------------------------
    //	Defining Classes
    //------------------------------------------------------------------------------
    class BaseGraphProperty {

    private:

    protected:

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void		    _init( double &width, double &height );

    public:

        Coord2 *                    centerPtr;
        const double *              widthPtr;
        const double *              heightPtr;

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        BaseGraphProperty( void );
        // copy constructor
        BaseGraphProperty( const BaseGraphProperty & c ) {}
        // destructor
        virtual ~BaseGraphProperty( void ) {}

        //------------------------------------------------------------------------------
        //	Assignment operators
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Reference to elements
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------

        void        init( double &__width, double &__height ) {
            _init( __width, __height );
        }

        //------------------------------------------------------------------------------
        //	Friend functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const BaseGraphProperty & v );
        // input
        friend istream &	operator >> ( istream & s, BaseGraphProperty & v );
        // class name
        virtual const char * className( void ) const { return "BaseGraphProperty"; }

    };

} // namespace Graph

#endif // _Graph_BaseGraphProperty_H
