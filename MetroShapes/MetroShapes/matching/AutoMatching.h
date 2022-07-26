//
// Created by tobias batik on 28.04.21.
//

#ifndef ZUKAI_AUTOMATCHING_H
#define ZUKAI_AUTOMATCHING_H

#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>

using namespace std;

#include "Metro.h"
#include "Guide.h"

#define FREDIST
// #define PARTMATCH
// #define TESTSIMILARITYFORALLPATHS
#define NORMALIZEPARTMATCHING
#define LESSSECTIONS
// #define DEBUG

#define REFPOINT


class AutoMatching {
private:
    Metro*                        _metro;
    Guide*                         _guide;

    unsigned int                    _nVertices;
    double                          MAX_COST = 10e6;
    double                          _tolerance;
    double                          _costMetro = 1.0;
    double                          _costAdd = 3.0;
    double                          _weightColor = .3;
    double                          _tolerancePartInprovment = .5;

    bool                            _deformUnuniformly = false;

    vector<VertexDescriptor> _foundMetroPath;
    vector<VertexDescriptor> _guideVds;

    // Direction Based Frechet distance
    vector<vector<double>> M;
    vector<Coord2> metroPathCd;
    vector<VertexDescriptor> _metroPathVD;
    vector<Coord2> guideCD;
    vector<double> metroVertPar;
    vector<double> guideVertPar;

    pair<VertexDescriptor, double> getShortestDistanceInQuee(set<VertexDescriptor> vertexQuee, map<VertexDescriptor, double> distances );
    map<VertexDescriptor, EdgeDescriptor> getAllNeighborsInQuee(VertexDescriptor selectedVertex, set<VertexDescriptor> vertexQuee);
    map<VertexDescriptor, EdgeDescriptor> getAllNeighboarsNotRemoved(VertexDescriptor selectedVertex, set<VertexDescriptor> vertexQuee ,set<VertexDescriptor> notSetVertex);
    map<VertexDescriptor, EdgeDescriptor> getAllNeighboarsNotRemovedOrStart(VertexDescriptor selectedVertex, set<VertexDescriptor> notSetVertex, VertexDescriptor vi);
    vector<VertexDescriptor> getPathVertex(VertexDescriptor endVert,
                                           map<VertexDescriptor, VertexDescriptor> previousSelected);
    vector<double> convertToIPR(vector<VertexDescriptor> vds, UndirectedGraph graph);
    vector<double> convertToIPR(vector<VertexDescriptor> vds, UndirectedGraph graph, bool final);
    vector<Coord2> convertVertDisToCoord(vector<VertexDescriptor> vds, UndirectedGraph graph);
    pair<double, double> pathMatchingIPRScore(vector<VertexDescriptor> mapPathVD, vector<VertexDescriptor> shapeVD, double tolerance);

    double computeDFD(vector<VertexDescriptor> mapPathVD, vector<VertexDescriptor> guideVD);
    pair<double, int> computePartDFD(vector<VertexDescriptor> mapPathVD, vector<VertexDescriptor> guideVD, int minMatch);
    // double computeDFDIterative(vector<Coord2> metroPathCd, vector<double> metroPathPar, vector<Coord2> guideCd, vector<double> guidePar);
    double computeSegmentDFD(int i, int j); // c(i,j)
    double costSegmentDFD(int i, int j);
    int pathNumberOfColors(vector<VertexDescriptor> path);
    Coord2 getCenterOfPath(vector< VertexDescriptor> path, UndirectedGraph g);


    bool pathSharePrefix(vector<double> mapPathIPR, vector<double> shapeIPR, double tolerance);
    int longestPrefix(vector<double> mapPathIPR, vector<double> shapeIPR, double tolerance);

    // Boundingbox
    Coord2 getMinValues(vector<VertexDescriptor> vds, UndirectedGraph graph);
    Coord2 getMaxValues(vector<VertexDescriptor> vds, UndirectedGraph graph);


protected:

public:
    void init( Metro * __metro, Guide * __guide, double __half_width, double __half_height);
    void findPath();
    void translateGuide();
    pair<vector<VertexDescriptor>, double> spaceDijgstar(VertexDescriptor vi);
    double getCostMetro() { return _costMetro; }
    double getCostAdd() { return _costAdd; }
    double getTolerancePartInprovment() { return _tolerancePartInprovment; }
    double getWeightColor() { return _weightColor; }

    void setCostMetro(double c) { _costMetro = c; }
    void setCostAdd( double c) { _costAdd = c; }
    void setWeightColor( double w ) { _weightColor = w; }
    void setTolerancePartInprovment( double t ) { _tolerancePartInprovment = t; }
    void pathFromMetroLine( int l);



    AutoMatching(){};
    virtual ~AutoMatching() {};

};


#endif //ZUKAI_AUTOMATCHING_H
