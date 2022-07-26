//******************************************************************************
// BaseEdgeProperty.h
//	: header file for base edge property
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Dec 27 23:16:12 2018
//
//******************************************************************************

#ifndef	_Graph_BaseEdgeProperty_H
#define _Graph_BaseEdgeProperty_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <vector>

using namespace std;

#include "Coord2.h"

//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace Graph {

    //------------------------------------------------------------------------------
    //	Defining Classes
    //------------------------------------------------------------------------------
    class BaseEdgeProperty {

    private:

    protected:

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void		    _init( void );

    public:

        unsigned int                id;

        double                      angle;
        double                      weight;
        bool                        visit;
        int                         visitedTimes;

        bool                        isFore;
        bool                        isBack;
        bool                        isShow;
        bool                        isCurved;
        bool                        isMetro;
        bool                        matchPath;
        bool                        isSmoothPath;
        bool                        originalMetroEdge = true;
        bool                        isVisible = true;

        double                      target;
        bool                        targetSet = false;
        double                      geoAngle;               // geographical angle
        double                      smoAngle;               // smooth angle
        double                      curAngle;               // current angle
	    vector< unsigned int >      lineID;                 // line id
        vector< string >            nameStationsOnEdge;
	    int                         stationsOnEdge = 0;

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        BaseEdgeProperty( void );
        // copy constructor
        BaseEdgeProperty( const BaseEdgeProperty & e ) {
            id	    = e.id;
            weight	= e.weight;
        }
        // destructor
        virtual ~BaseEdgeProperty( void ) {}

        //------------------------------------------------------------------------------
        //	Assignment operators
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Reference to elements
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------

        void        init( void )		      { _init(); }

        //------------------------------------------------------------------------------
        //	Friend functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const BaseEdgeProperty & v );
        // input
        friend istream &	operator >> ( istream & s, BaseEdgeProperty & v );
        // class name
        virtual const char * className( void ) const { return "BaseEdgeProperty"; }

    };

} // namespace Graph

#endif // _Graph_BaseEdgeProperty_H
