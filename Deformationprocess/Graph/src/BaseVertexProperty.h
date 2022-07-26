//******************************************************************************
// BaseVertexProperty.h
//	: header file for base vertex property
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Aug 13 23:16:12 2020
//
//******************************************************************************

#ifndef	_Graph_BaseVertexProperty_H
#define _Graph_BaseVertexProperty_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <vector>

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
    class BaseVertexProperty {

    private:

    protected:

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void		    _init( void );

    public:

        unsigned int                id;
        string *                    namePtr;
        double *                    namePixelWidthPtr;      // pixel width of the name
        double *                    namePixelHeightPtr;     // pixel height of the name
        double                      weight;                 // vertex weight

        Coord2 *                    geoPtr;                 // center coordinates
	    Coord2 *                    smoothPtr;              // center coordinates
	    Coord2 *                    octilinearPtr;          // center coordinates
	    Coord2 *                    coordPtr;               // center coordinates
        double *                    widthPtr;               // vertex width
        double *                    heightPtr;              // vertex height
        double *                    areaPtr;                // vertex area
       


	
	    vector< unsigned int >      lineID;                 // line id
        int                         color;                  // color type
        bool                        flag;                   // flag
        bool                        metroShape;             // true if vertex part of metro Shape 
        bool                        validPos; 
        bool                        inflectionPoint;
        bool                        intersection = false;
        int                         intersectionsCount = 0;
        double                      smoothDistance = 0;
        double                      initDistance = 0;
        Coord2                      closestPointOnEdge;
        bool                        autoPath = false;       // the automatically found path in geo map
        bool                        smoothPath = false;     // path after the smooth stage
        bool                        isStation = true;
        bool                        collision = false;
        Coord2                      closOnMetro = Coord2(0 , 0);
        bool                        interstate;
        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        BaseVertexProperty( void );
        // copy constructor
        BaseVertexProperty( const BaseVertexProperty & c ) {}
        // destructor
        virtual ~BaseVertexProperty( void ) {}

        //------------------------------------------------------------------------------
        //	Assignment operators
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Reference to elements
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------

        void        init( void )		{ _init(); }

        //------------------------------------------------------------------------------
        //	Friend functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const BaseVertexProperty & v );
        // input
        friend istream &	operator >> ( istream & s, BaseVertexProperty & v );
        // class name
        virtual const char * className( void ) const { return "BaseVertexProperty"; }

    };

} // namespace Graph

#endif // _Graph_BaseVertexProperty_H
