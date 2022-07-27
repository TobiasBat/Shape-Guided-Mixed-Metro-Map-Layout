/*
*	Smooth.h
*   Based on FocusContext/Smooth.h
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   Calculating Smooth integrating guide Layout 
*/

#ifndef _Smooth_H        // begining of header file
#define _Smooth_H        // notifying that this file is included


//----------------------------------------------------------------------
//  Including header files
//----------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

#include "Metro.h"
#include "Guide.h"

//------------------------------------------------------------------------------
//	Defining data types
//------------------------------------------------------------------------------
#define SMOOTH_CONFLICT
#define SMOOTH_BOUNDARY

#define ALONGEDGE
#define SHORTESTANGLE
#define PARALLEL
#define NORMALPAR
#define OVERLAPPING
#define WEIGHT_SHAPEAPPROXIMATION_BYDEGREE

//----------------------------------------------------------------------
//	Defining macros
//----------------------------------------------------------------------

class Smooth
{
private:

    Metro           * _metro;
    Guide           * _guide; 
    Eigen::VectorXd _var;           // x  
    Eigen::VectorXd _output;        // b  output = (var * weights)
    Eigen::MatrixXd _coef;          // A
    double          _half_width;    // window_width
    double          _half_height;   // window_height

    unsigned int    _nVars;
    unsigned int    _nConstrs;

    double          _w_labelangle;
    double          _d_Alpha;
    double          _targetLength;
    bool            _pathMode;


    Coord2 closestPointOnEdge( Coord2 source, Coord2 target, Coord2 vertex );
    int calculateIntersections( Coord2 vertexCord, Coord2 closestPointOnGuide );
    double magnitude( Coord2 v );

    void _calculateAlongEdgeParameter();
    void _handleHigherDegreeVertex();
    void _removeSingleNonItersectingVertex();
    void _asignSmoothEdge();

protected:

    void            _setVars        ( unsigned int & nRows );
    void            _setConstraints ( unsigned int & nRows );
    void            _initVars       ( void );
    void            _initCoefs      ( void );
    void            _initOutputs    ( void );
    void            _updateCoefs    ( void );
    void            _updateOutputs  ( void );
    virtual void    _init           ( Metro * __metro, Guide * __guide, double __width, double __height );


    
public:
    double          _w_alongEdge = 4.;                          // w_c
    double          _w_angle = 2.0;                             // w_a
    double          _w_position = sqrt(0.025);                  // w_p
    double          _w_contextlength = 1.0; // sqrt(3.0);       // w_l
    double          _w_boundary = sqrt(4.0);
    double          _w_crossing = sqrt(10.0);
    double          _gamma = 100.0;                             // legacy Not Used Anymore
    double          _constSmooth = 0.08;                        // legacy Not Used Anymore
    double          _w_overlapping = 0.;                        // legacy Not Used Anymore
    double          _w_parallel = .0;                           // legacy Not Used Anymore
    bool            _exclusivlyPathMode = false;
    
    void            _resetSmoothPtr( void );
    
    Smooth();                     // default constructor
    Smooth( const Smooth & obj ); // Copy constructor
    virtual ~Smooth();            // Destructor

//------------------------------------------------------------------------------
//  Reference to members
//------------------------------------------------------------------------------
    
//------------------------------------------------------------------------------
//  Specific functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//      Initialization functions
//------------------------------------------------------------------------------
    void prepare( Metro * __metro, Guide * __guide ,double __width, double __height ) {
        _init(__metro, __guide, __width, __height ); 
    }

//------------------------------------------------------------------------------
//  File I/O
//------------------------------------------------------------------------------
    void prepare( void );
    void clear( void );
    void retrieve( void );
    double LeastSquare( unsigned int iter );
    double ConjugateGradient( unsigned int iter );

//------------------------------------------------------------------------------
//      I/O
//------------------------------------------------------------------------------
    friend ostream & operator << ( ostream & stream, const Smooth & obj );
                                // Output
    friend istream & operator >> ( istream & stream, Smooth & obj );
                                // Input

    virtual const char * className( void ) const { return "Smooth"; }
                                // Class name
};

#endif // _Smooth_H

// end of header file
// Do not add any stuff under this line.

