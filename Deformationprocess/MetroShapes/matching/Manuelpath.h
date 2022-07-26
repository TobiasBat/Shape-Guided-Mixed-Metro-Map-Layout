//
// Created by tobias batik on 26.09.21.
//

#ifndef ZUKAI_MANUELPATH_H
#define ZUKAI_MANUELPATH_H

#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Metro.h"
#include "Guide.h"

class Manuelpath {

private:
    Metro                       * _metro;
    Guide                       * _guide;
    vector<VertexDescriptor>      _path;
    vector<int>                   _pathId;
    double                        _maxDistanceVertex;
    bool                          _deformUnuniformly = true;

    Coord2                      getCenterPath();

public:
    void init( Metro * __metro, Guide * __guide );
    bool addVertex( int id );
    bool addVertex( Coord2 coord );
    bool addVertices( vector<int> ids );
    void alingGuide();
    void scaleGuide();
    void scaleMetroNonUniformly();

};


#endif //ZUKAI_MANUELPATH_H
