/*
*	Mixedlayout.h
*   Based on FocusContext/Octilinear.h
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   To calculate the Mixed layout
*/

#ifndef _Octilinear_H        // begining of header file
#define _Octilinear_H        // notifying that this file is included

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
#define OCTILINEAR_CONFLICT
#define OCTILINEAR_BOUNDARY
#define GUID_POSITION
#define DEGREE_TWO_HEURISTIC
// #define OVERLAPPIN_CONST

//----------------------------------------------------------------------
//	Defining macros
//----------------------------------------------------------------------

class Mixedlayout
{
private:

    Metro               * _metro;
    Guide               * _guide; 
    Eigen::VectorXd     _var;           // x
    Eigen::VectorXd     _output;        // b
    Eigen::MatrixXd     _coef;          // A
    double              _half_width;    // window_width
    double              _half_height;   // window_height    

    unsigned int        _nVars;
    unsigned int        _nConstrs;
    
    
    double              _d_Alpha;
    double              _d_Beta;
    vector< double >    _theta;         // closest octilinear theta

    // Sector Assignment
    double getOffsetToSector(EdgeDescriptor edge, int sector);
    bool continioueAssigning(vector<vector<EdgeDescriptor>> assignments);
   
protected:

    void                _setVars        ( unsigned int & nRows );
    void                _setConstraints ( unsigned int & nRows );
    void                _initCoefs      ( void );
    void                _initVars       ( void );
    void                _initOutputs    ( void );
    void                _updateCoefs    ( void );
    void                _updateOutputs  ( void );
    virtual void        _init           ( Metro * __metro, Guide * __guid,double __width, double __height );
    void                _setTargetAngle( void );
    void                _setTargetAngle2( void );
    void                _updateEdgeCurAngle( void );
    double              magnitude( Coord2 v );
    void                _checkSameAnge();
    vector<pair<EdgeDescriptor, double>> _solveTargetAnglesConflict(vector<EdgeDescriptor> edges, vector<double> curAngles);
    void                removeStation(VertexDescriptor vertex);
    bool containsASetEdge(vector<EdgeDescriptor> edges );



public:
    double              _w_mixed            = 2.0; // was 5     // w_o
    double              _w_position         = .1;   // w_p
    double              _w_position_path    = 10.0; // was 20 //w_c
    double              _w_boundary         = sqrt( 20.0 );
    double              _w_crossing         = 14.0; //sqrt( 20.0 );
    double              _w_overlap          = 4.0;  // legacy not Used anymore
    double              _gama               = 4.0;  // legacy not Used anymore
    double              _constMixed         = 2.0;                   // C_m // legacy not used anymore
    double              _minAngle           = M_PI_4 *0.5;            // legacy not used anymore



    Mixedlayout();                     // default constructor
    Mixedlayout( const Mixedlayout & obj ); // Copy constructor
    virtual ~Mixedlayout();            // Destructor
//------------------------------------------------------------------------------
//      Initialization functions
//------------------------------------------------------------------------------
    void prepare( Metro * __metro, Guide * __guid,double __half_width, double __half_height ) {
        _init( __metro, __guid ,__half_width, __half_height );
    }

//------------------------------------------------------------------------------
//  File I/O
//------------------------------------------------------------------------------
    void prepare( void );
    void clear( void );
    void retrieve( void );
    double LeastSquare( unsigned int iter );
    double ConjugateGradient( unsigned int iter );
    void postProcessing( void );
    void complexPostProcessing( void );
    void preparePostProcessing();

//------------------------------------------------------------------------------
//      I/O
//------------------------------------------------------------------------------
    friend ostream & operator << ( ostream & stream, const Mixedlayout & obj );
                                // Output
    friend istream & operator >> ( istream & stream, Mixedlayout & obj );
                                // Input

    virtual const char * className( void ) const { return "Mixedlayout"; }
                                // Class name
};

#endif // _Octilinear_H

// end of header file
// Do not add any stuff under this line.

