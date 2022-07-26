/*
*	Stationmatching.h
*	@author Tobias Batik
*	@version 1.0  14/03/2021
*
*   Creating matching, and transforming guide
*   based on matching
*/

//----------------------------------------------------------------------
//  Including header files
//----------------------------------------------------------------------
#ifndef _Stationmatching_H        // begining of header file
#define _Stationmatching_H        // notifying that this file is included

#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector> 

using namespace std;

#include "Metro.h"
#include "Guide.h"

class Stationmatching {
    private: 
        Metro*                        _metro; 
        Guide*                         _guide; 
        double                        _half_width; 
        double                        _half_height; 
        double                        distTol = 20; 

        Eigen::VectorXd                 _var;       // x 
        Eigen::VectorXd                 _output;    // b  output = (var * weights)
        Eigen::MatrixXd                 _coef;      // A
        unsigned int                    _nVars; 
        unsigned int                    _nVertices; 
        unsigned int                    _nConstrs; 
        unsigned int                    _nDifConstrs;
        double                          _d_Alpha; 
        double                          _d_Beta; 

        vector<VertexDescriptor>      _guideDisc; 
        vector<VertexDescriptor>      _stationDisc; 
        vector<Coord2>                _guideCoord; 
        vector<Coord2>                _stationCoord; 
        
        double                          _wTransX; 
        double                          _wTransY; 
        double                          _wRotate;
        double                          _wScale; 

        void _initVars( void ); 
        void _initCoefs( void ); 
        void _initOutputs( void ); 
        void _updateOutputs( void ); 
        
        double magnitude(Coord2 v); 
        Coord2 getAsPolar(Coord2 center, Coord2 point); 
        Coord2 polarAsCart(Coord2 center, Coord2 point);
        Coord2 getAllCenters( vector<Coord2> points); 
        Coord2 getCenterOfGuideNodes( void ); 

    protected: 
        
        
    public: 
        void init( void ); 
        void prepare( Metro * __metro, Guide * __guide,  double __half_width, double __half_height ); 
        void retrieve( void ); 

        int getNumberMatchings( void ) { return _stationDisc.size(); }
        bool addNewNode(Coord2 cord); 
        bool removeNode(Coord2 cord);

        vector<Coord2> getStationCoord() { return _stationCoord; }
        vector<Coord2> getGuideCoord() { return _guideCoord; }
        vector<VertexDescriptor> getGuides() { return _guideDisc; }
        vector<VertexDescriptor> getStations() { return _stationDisc; }
        void distortMap( void ); 

        void conjugateGradient( void );

        Stationmatching();
        virtual ~Stationmatching();            // Destructor
};

#endif