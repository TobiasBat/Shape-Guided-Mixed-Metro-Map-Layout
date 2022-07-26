//==============================================================================
// Metro.cc
//  : program file for the metro network
//
//------------------------------------------------------------------------------
//
//              Date: Mon Dec 10 04:28:26 2012
//
//==============================================================================

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cmath>

using namespace std;

#include "Metro.h"

//------------------------------------------------------------------------------
//	Protected functions
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------

//
//  Metro::Metro -- default constructor
//
//  Inputs
//  none
//
//  Outputs
//  none
//
Metro::Metro( void )
{
    clearGraph( graph );

    _shortestPathM.clear();
    _line.clear();
    _lineSta.clear();
    // _lineColor.clear();
    // _lineName.clear();

    _nLines = 0;
    _nStations  = 0;
    _nEdges = 0;
    _meanVSize = 0.0;

    _removedVertices.clear();
    _removedEdges.clear();
    _removedWeights.clear();

    _VEconflict.clear();
    _ratioR.clear();
    _ratioGeo.clear();
}

//
//  Metro::Metro -- copy constructor
//
//  Inputs
//  obj : object of this class
//
//  Outputs
//  none
//
Metro::Metro( const Metro & obj )
{
    graph = obj.g();

    _shortestPathM = obj.spM();
    _line = obj.line();
    _lineSta = obj.lineSta();

    for ( unsigned int k = 0; k < _nLines; ++k ) {
        for ( unsigned i = 0; i < 3; ++i )
            _lineColor[ k ][ i ] = obj._lineColor[ k ][ i ];
        strcpy( _lineName[ k ], obj._lineName[ k ] );
    }

    _nLines     = obj._nLines;
    _nStations  = obj._nStations;
    _nEdges     = obj._nEdges;
    _meanVSize = obj._meanVSize;

    _removedVertices = obj.removedVertices();
    _removedEdges   = obj.removedEdges();
    _removedWeights = obj.removedWeights();

    _VEconflict = obj.VEconflict();
    _ratioR     = obj.ratioR();
    _ratioGeo  = obj.ratioGeo();
}


//------------------------------------------------------------------------------
//	Destructor
//------------------------------------------------------------------------------

//
//  Metro::~Metro --    destructor
//
//  Inputs
//  none
//
//  Outputs
//  none
//
Metro::~Metro( void )
{
    clearGraph( graph );

    _shortestPathM.clear();
    _line.clear();
    _lineSta.clear();

    _nLines = 0;
    _nStations  = 0;
    _nEdges  = 0;
    _meanVSize = 0.0;

    _removedVertices.clear();
    _removedEdges.clear();
    _removedWeights.clear();

    _VEconflict.clear();
    _ratioR.clear();
    _ratioGeo.clear();
}

//
//  Metro::clear --    clear the current metro information
//
//  Inputs
//  none
//
//  Outputs
//  none
//
void Metro::clear( void )
{
    clearGraph( graph );

    _shortestPathM.clear();
    _line.clear();
    _lineSta.clear();

    _nLines = 0;
    _nStations  = 0;
    _nEdges  = 0;
    _meanVSize = 0.0;

    _removedVertices.clear();
    _removedEdges.clear();
    _removedWeights.clear();

    _VEconflict.clear();
    _ratioR.clear();
    _ratioGeo.clear();
}

void Metro::removeUnecessaryTempVertex() {
    double sector[ 9 ] = { -M_PI, -3.0*M_PI/4.0, -M_PI/2.0, -M_PI/4.0, 0.0,
                           M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI };
    vector<VertexDescriptor> stationsToRemove;
    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
        DegreeSizeType degree = out_degree( vd, graph );
        if (!graph[vd].isStation and degree == 2) {
            OutEdgeIterator e, e_end;
            float sumOffset = 0;
            vector<VertexDescriptor> newEdgeVertices;

            for ( tie( e, e_end ) = out_edges( vd, graph ); e != e_end; ++e ) {
                EdgeDescriptor ed = *e;
                double angle = graph[ed].curAngle;
                float minOffset = 10e6;
                for (auto s : sector) {
                    if (fabs(angle - s) < minOffset)
                        minOffset = fabs(angle - s);
                }
                sumOffset += minOffset;

                if (graph[source(ed, graph)].id != graph[vd].id) {
                    newEdgeVertices.push_back(source(ed, graph));
                } else if (graph[target(ed, graph)].id != graph[vd].id) {
                    newEdgeVertices.push_back(target(ed,graph));
                }
            }

            // calc angle when vertex not here
            double minWithoutOffset = 10e6;
            if (newEdgeVertices.size() == 2) {
                Coord2 vjPos = *graph[newEdgeVertices[0]].coordPtr;
                Coord2 vkPos = *graph[newEdgeVertices[1]].coordPtr;
                Coord2 delta = Coord2( vjPos.x() - vkPos.x(),
                                       vjPos.y() - vjPos.y());

                double angle = atan2(delta.y(), delta.x());
                for (auto s : sector) {
                    if (fabs(angle - s) < minWithoutOffset) {
                        minWithoutOffset = fabs(angle - s);
                    }
                }


                // if better remove vertex
                if (minWithoutOffset < sumOffset * 0.5) {
                    stationsToRemove.push_back(vd);
                }else {
                    cout << "Will not be removed" << endl;
                }

            }




        }

    }

    for (auto vd : stationsToRemove) {
        removeStation(vd);
    }
    reorderID();
    checkVEConflicts();
}

int Metro::nonPlanarIntersections() {
    int result = 0;
    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        if (intersectingEdge(edge).size() > 0)
            result++;
    }

    return result;
}

Coord2 Metro::intersectionPoint(Coord2 vi, Coord2 vj, Coord2 ui, Coord2 uj)
{


    return Coord2(0, 0);
}

vector<pair<EdgeDescriptor, Coord2> > Metro::intersectingEdge( EdgeDescriptor inputEdge ) {

    vector<pair<EdgeDescriptor, Coord2> > intersectingEdges;

    VertexDescriptor v1D = source( inputEdge, graph );
    VertexDescriptor v2D = target( inputEdge, graph );
    Coord2 v1 = *graph[ v1D ].coordPtr;
    Coord2 v2 = *graph[ v2D ].coordPtr;

    double x1 = v1.x();
    double y1 = v1.y();
    double x2 = v2.x();
    double y2 = v2.y();

    double A1 = y2 - y1;
    double B1 = x1 - x2;
    double C1 = A1 * x1 +B1 * y1;

    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        Coord2 edgeSource = *graph[source(edge, graph)].coordPtr;
        Coord2 edgeTarget = *graph[target(edge, graph)].coordPtr;

        if (edgeSource.x() != v1.x() && edgeTarget.x() != v1.x()
        // && graph[inputEdge].id < graph[edge].id
        ) { // in result every intersection only once
            double x3 = edgeSource.x();
            double y3 = edgeSource.y();
            double x4 = edgeTarget.x();
            double y4 = edgeTarget.y();

            double A2 = y4 - y3;
            double B2 = x3 - x4;
            double C2 = A2 * x3 + B2 * y3;

            double x,y;
            double det = A1 * B2 - A2 * B1;
            if ( fabs(det) < 0.0001 ) { //Lines are parallel
                x = y = 0;
            } else {
                x = (B2 * C1 - B1 * C2) / det;
                y = (A1 * C2 - A2 * C1) / det;
                Coord2 interseionPoint = Coord2(x,y);
                double distSource = (interseionPoint - edgeSource).norm();
                double distTarget = (interseionPoint - edgeTarget).norm();
                double minDist = min(distSource, distTarget);
                double tol = 0.00001;
                if (    (min(x1,x2) - tol <= x && x <= max(x1, x2) + tol ) &&
                (min(x3,x4) - tol <= x && x <= max(x3, x4) + tol ) &&
                (min(y1,y2) - tol <= y && y <= max(y1, y2) + tol ) &&
                (min(y3,y4) - tol <= y && y <= max(y3, y4) + tol) &&
                minDist > 0.001
                ) {
                    pair<EdgeDescriptor, Coord2> intersectionData = pair<EdgeDescriptor, Coord2>(edge, Coord2(x,y));
                    if (intersectingEdges.size() < 1)
                        intersectingEdges.push_back(intersectionData);
                    // graph[inputEdge].isMetro = false;
                    // graph[edge].isMetro = false;
                }
            }
        }
    }

    return intersectingEdges;
}

void Metro::simplifyMixedLayout()
{

    bool removedStation = true;
    int i = 0;
    while (removedStation) {
        removedStation = false;
        vector<VertexDescriptor> removeStations;
        vector<EdgeDescriptor> edgesWillRemoved;
        BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
            DegreeSizeType degree = out_degree( vd, graph );
            if (degree == 2) {
                OutEdgeIterator e, e_end;
                double target = M_PI * 3.0;
                bool remove = true;
                for ( tie( e, e_end ) = out_edges( vd, graph ); e != e_end; ++e ) {
                    EdgeDescriptor ed = *e;
                    if (target == M_PI * 3.0) {
                        target = graph[ed].target;
                    } else if (target != graph[ed].target and target != graph[ed].target + M_PI and target != graph[ed].target - M_PI) {
                        remove = false;
                    } else if (count(edgesWillRemoved.begin(), edgesWillRemoved.end(), ed)) { // edge will be removed not add vertex to remove list
                        remove = false;
                    } else if (graph[ed].isSmoothPath)
                        remove = false;
                }

                if (remove) {
                    removeStations.push_back(vd);
                    removedStation = true;
                }

            }


        }

        for (auto vd : removeStations)
            removeStation(vd);
        i += 1;
    }

    reorderID();
    checkVEConflicts();
}

void Metro::removeStation(VertexDescriptor vd)
{
    OutEdgeIterator e, e_end;
    vector<EdgeDescriptor> edgesRemove;
    vector<VertexDescriptor> newEdgesVertex;
    int stationsOnNewEdge = 0;
    for ( tie( e, e_end ) = out_edges( vd, graph ); e != e_end; ++e ) {
        EdgeDescriptor ed = *e;
        VertexDescriptor vS = source( ed, graph );
        VertexDescriptor vT = target( ed, graph );
        if (vS != vd)
            newEdgesVertex.push_back(vS);
        else if (vT !=  vd)
            newEdgesVertex.push_back(vT);
        edgesRemove.push_back(ed);
    }
    vector< unsigned int > lineID;
    double target, geoAngle, smoAngle, curAngle;
    string stationName = *graph[vd].namePtr;

    vector<string> prevNames;

    for (auto ed : edgesRemove) {
        lineID = graph[ed].lineID;
        target = graph[ed].target;
        geoAngle = graph[ed].geoAngle;
        smoAngle = graph[ed].smoAngle;
        curAngle = graph[ed].curAngle;
        stationsOnNewEdge += graph[ed].stationsOnEdge;
        for (string n : graph[ed].nameStationsOnEdge) {
            prevNames.push_back(n);
        }
        remove_edge(ed, graph);
    }
    if (graph[vd].isStation)
        stationsOnNewEdge += 1;
    remove_vertex(vd, graph);

    
    if (newEdgesVertex.size() == 2) {

        add_edge(newEdgesVertex[0], newEdgesVertex[1], graph);
        if ( boost::edge(newEdgesVertex[0], newEdgesVertex[1], graph).second) {
            auto edge = boost::edge(newEdgesVertex[0], newEdgesVertex[1], graph).first;
            _nEdges++;
            graph[ edge ].id          = _nEdges;
            graph[ edge ].weight      = 1.0 + stationsOnNewEdge;
            graph[ edge ].geoAngle    = geoAngle;
            graph[ edge ].smoAngle    = smoAngle;
            graph[ edge ].curAngle    = curAngle;
            graph[ edge ].lineID      = lineID;
            graph[ edge ].isCurved    = false;
            graph[ edge ].isMetro     = true;
            graph[ edge ].matchPath   = false;
            graph[ edge ].isSmoothPath = false;
            graph[ edge ].target      = target;
            graph[ edge ].stationsOnEdge = stationsOnNewEdge;
            for (auto n : prevNames) {
                graph[ edge ].nameStationsOnEdge.push_back(n);
            }
            graph[ edge ].nameStationsOnEdge.push_back(stationName);
        }
    }
}

void Metro::subdivideComplexEdges() {
    // Subdivide edges that probably can not rotatet to octolinear later
    vector<EdgeDescriptor> subdivideEdges;
    BGL_FORALL_EDGES( edge, graph, UndirectedGraph )
    {
        VertexDescriptor sourceVd = source(edge, graph);
        VertexDescriptor targetVd = target(edge, graph);
        DegreeSizeType degreeSource = out_degree( sourceVd, graph );
        DegreeSizeType degreeTarget = out_degree( targetVd, graph );
        if (degreeSource > 3 and degreeTarget > 3 and
        (!graph[sourceVd].autoPath or !graph[targetVd].autoPath) and
        (graph[sourceVd].isStation and graph[targetVd].isStation)
        // and (*graph[sourceVd].coordPtr - *graph[targetVd].coordPtr).norm() > 7.0
        )// _distBeta * .5)
            subdivideEdges.push_back(edge);
    }

    for (auto edge : subdivideEdges) {
        subdivideEdge(edge);
    }
    for (auto edge : subdivideEdges) {
        remove_edge(edge, graph);
    }
    reorderID();
    cout << "Subdiving Complex Edges done " << endl;


    // Label instates
    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ) {
        DegreeSizeType degree = out_degree( vd, graph );
        if (degree == 0)
            cerr << "degree 0 exists" << *graph[vd].namePtr << endl;
    }
}

void Metro::subdivideEdge(EdgeDescriptor edge) {
    cout << "   " << "subdividing an Edge" << endl;
    VertexDescriptor sVd = source(edge, graph);
    VertexDescriptor tVd = target(edge, graph);
    Coord2 sourcePos = *graph[sVd].coordPtr;
    Coord2 targetPos = *graph[tVd].coordPtr;

    Coord2 pos =  Coord2(0.5 * sourcePos.x(), 0.5 * sourcePos.y());

    pos += Coord2(0.5 * targetPos.x(), 0.5 * targetPos.y());


    // create new station
    _nStations++;
    VertexDescriptor vDQ = add_vertex( graph );
    graph[ vDQ ].id               = _nStations;
    graph[ vDQ ].coordPtr         = new Coord2( pos.x(), pos.y() );
    graph[ vDQ ].geoPtr           = new Coord2( pos.x(), pos.y() );
    graph[ vDQ ].smoothPtr        = new Coord2( pos.x(), pos.y() );
    graph[ vDQ ].octilinearPtr    = new Coord2( pos.x(), pos.y() );
    graph[ vDQ ].namePtr          = new string( "complexity d Sta");
    graph[ vDQ ].metroShape       = true;
    graph[ vDQ ].inflectionPoint  = false;
    if (graph[sVd].autoPath and graph[tVd].autoPath)
        graph[vDQ].autoPath = true;
    else
        graph[vDQ].autoPath = false;
    if (graph[sVd].smoothPath and graph[tVd].smoothPath)
        graph[vDQ].smoothPath = true;
    else
        graph[ vDQ ].smoothPath       = false;
    graph[ vDQ ].interstate       = false;
    graph[ vDQ ].isStation        = false;
    graph[ vDQ ].weight           = 1.0;
    graph[ vDQ ].lineID = graph[ sVd ].lineID;

    double angle; // = graph[edge].geoAngle;
    double diffX = sourcePos.x() - pos.x();
    double diffY = sourcePos.y() - pos.y();
    angle = atan2(diffY, diffX);

    add_edge( sVd, vDQ, graph);
    auto eI = boost::edge( sVd, vDQ, graph).first;
    _nEdges++;
    graph[ eI ].id          = _nEdges;
    graph[ eI ].weight      = 1.0;
    graph[ eI ].geoAngle    = angle;
    graph[ eI ].smoAngle    = angle;
    graph[ eI ].curAngle    = angle;
    graph[ eI ].lineID      = graph[edge].lineID;
    graph[ eI ].isCurved    = false;
    graph[ eI ].isMetro     = true;
    graph[ eI ].matchPath   = graph[edge].matchPath;
    graph[ eI ].isSmoothPath = graph[edge].isSmoothPath;
    graph[ eI ].target      = graph[edge].target;

    if (graph[edge].stationsOnEdge > 0) {
        graph[ eI ].stationsOnEdge = ((int) graph[edge].stationsOnEdge / 2) + 1;
    }
    else {
        graph[ eI ].stationsOnEdge = 0;
    }

    add_edge(tVd, vDQ, graph);
    auto eJ = boost::edge(tVd, vDQ, graph).first;

    diffX = targetPos.x() - pos.x();
    diffY = targetPos.y() - pos.y();
    angle = atan2(diffY, diffX);

    // cout << "     " << boost::edge(vDQ, tVd, graph).second << endl;
    _nEdges++;
    graph[ eJ ].id          = _nEdges;
    graph[ eJ ].weight      = .5;
    // double geoAngle = graph[edge].geoAngle;
    angle += 0.000001;
    graph[ eJ ].geoAngle    = angle;
    graph[ eJ ].smoAngle    = angle;
    graph[ eJ ].curAngle    = angle;
    graph[ eJ ].lineID      = graph[edge].lineID;
    graph[ eJ ].isCurved    = false;
    graph[ eJ ].isMetro     = true;
    graph[ eJ ].matchPath   = graph[edge].matchPath;
    graph[ eI ].isSmoothPath = graph[edge].isSmoothPath;
    graph[ eJ ].target      = graph[edge].target;

    if (graph[edge].stationsOnEdge > 0) {
        graph[ eJ ].stationsOnEdge = (int) graph[edge].stationsOnEdge / 2;
    }
    else {
        graph[ eJ ].stationsOnEdge = 0;
    }

    // remove_edge(edge, graph);
}

void Metro::subdiveToOctilinear() {
    cout << "Subdividing until octolinear" << endl;
    vector<EdgeDescriptor> edgesRemove;
    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        if (!graph[edge].isSmoothPath) {
            auto targetAngle = graph[edge].target;
            auto curAngle = graph[edge].curAngle;
            if (fabs(curAngle - targetAngle) > 0.04) { // thresshold what people can see?
                cout << "line is not octolinear" << endl;
                edgesRemove.push_back(edge);
                auto sourceDes = source(edge, graph);
                auto targetDes = target(edge, graph);
                if (graph[sourceDes].id > graph[targetDes].id ) {
                    // sourceDes = targetDes;
                    // targetDes = source(edge, graph);
                }
                Coord2 vI = *graph[sourceDes].coordPtr;
                Coord2 vJ = *graph[targetDes].coordPtr;
                // TODO make dependent on number of vertices
                Coord2 center = Coord2(vI.x() * 0.5 + vJ.x() * 0.5,
                                       vI.y() * 0.5 + vJ.y() * 0.5);

                Coord2 vICenter = (center - vI);
                double alpha = targetAngle - curAngle;
                Coord2 vICenterRot = Coord2(vICenter.x() * cos(alpha) - vICenter.y() * sin(alpha),
                                            vICenter.x() * sin(alpha) + vICenter.y() * cos(alpha));
                Coord2 pos = vICenterRot + vI;

                double xI = pos.x();
                double yI = pos.y();

                Coord2 vIdelta = vICenterRot - vICenter;
                Coord2 posJ = vICenter - vIdelta + vI;
                double xJ = posJ.x();
                double yJ = posJ.y();

                // create new station
                _nStations++;
                VertexDescriptor curVD = add_vertex( graph );
                graph[ curVD ].id               = _nStations;
                graph[ curVD ].coordPtr         = new Coord2( xI, yI );
                graph[ curVD ].geoPtr           = new Coord2( xI, yI );
                graph[ curVD ].smoothPtr        = new Coord2( xI, yI );
                graph[ curVD ].octilinearPtr    = new Coord2( xI, yI );
                graph[ curVD ].namePtr          = new string( "re-inserted station");
                graph[ curVD ].metroShape       = true;
                graph[ curVD ].inflectionPoint  = false;
                graph[ curVD ].smoothPath       = false;
                graph[ curVD ].interstate       = false;
                graph[ curVD ].isStation        = false;
                graph[ curVD ].weight           = 1.0;
                graph[ curVD ].lineID = graph[ sourceDes ].lineID;

                add_edge(sourceDes, curVD, graph);
                auto eI = boost::edge(sourceDes, curVD, graph).first;
                _nEdges++;
                graph[ eI ].id          = _nEdges;
                graph[ eI ].weight      = 1.0;
                graph[ eI ].geoAngle    = targetAngle;
                graph[ eI ].smoAngle    = targetAngle;
                graph[ eI ].curAngle    = targetAngle;
                graph[ eI ].lineID      = graph[edge].lineID;
                graph[ eI ].isCurved    = false;
                graph[ eI ].isMetro     = true;
                graph[ eI ].matchPath   = false;
                graph[ eI ].isSmoothPath = false;
                graph[ eI ].target      = targetAngle;
                if (graph[edge].stationsOnEdge > 0) {
                    graph[ eI ].stationsOnEdge = ((int) graph[edge].stationsOnEdge / 2) + 1;
                }
                else {
                    graph[ eI ].stationsOnEdge = 0;
                }

                // create new station
                _nStations++;
                VertexDescriptor vDQ = add_vertex( graph );
                graph[ vDQ ].id               = _nStations;
                graph[ vDQ ].coordPtr         = new Coord2( xJ, yJ );
                graph[ vDQ ].geoPtr           = new Coord2( xJ, yJ );
                graph[ vDQ ].smoothPtr        = new Coord2( xJ, yJ );
                graph[ vDQ ].octilinearPtr    = new Coord2( xJ, yJ );
                graph[ vDQ ].namePtr          = new string( "re-inserted station");
                graph[ vDQ ].metroShape       = true;
                graph[ vDQ ].inflectionPoint  = false;
                graph[ vDQ ].smoothPath       = false;
                graph[ vDQ ].interstate       = false;
                graph[ vDQ ].isStation        = false;
                graph[ vDQ ].weight           = 1.0;
                graph[ vDQ ].lineID = graph[ sourceDes ].lineID;


                add_edge(vDQ, targetDes, graph);
                auto eJ = boost::edge(vDQ, targetDes, graph).first;
                _nEdges++;
                graph[ eJ ].id          = _nEdges;
                graph[ eJ ].weight      = 1.;
                graph[ eJ ].geoAngle    = targetAngle;
                graph[ eJ ].smoAngle    = targetAngle;
                graph[ eJ ].curAngle    = targetAngle;
                graph[ eJ ].lineID      = graph[edge].lineID;
                graph[ eJ ].isCurved    = false;
                graph[ eJ ].isMetro     = true;
                graph[ eJ ].matchPath   = false;
                graph[ eJ ].isSmoothPath = false;
                graph[ eJ ].target      = targetAngle;
                if (graph[edge].stationsOnEdge > 0) {
                    graph[ eJ ].stationsOnEdge = (int) graph[edge].stationsOnEdge / 2;
                }
                else {
                    graph[ eJ ].stationsOnEdge = 0;
                }


                add_edge(vDQ, curVD, graph);
                eJ = boost::edge(vDQ, curVD, graph).first;
                _nEdges++;
                graph[ eJ ].id          = _nEdges;
                graph[ eJ ].weight      = 1.0;
                graph[ eJ ].geoAngle    = targetAngle;
                graph[ eJ ].smoAngle    = targetAngle;
                graph[ eJ ].curAngle    = targetAngle;
                graph[ eJ ].lineID      = graph[edge].lineID;
                graph[ eJ ].isCurved    = false;
                graph[ eJ ].isMetro     = true;
                graph[ eJ ].matchPath   = false;
                graph[ eJ ].isSmoothPath = false;
                graph[ eJ ].target      = targetAngle;
                graph[ eJ ].stationsOnEdge = 0;
            }
        }
    }
    for (auto edge : edgesRemove)
        remove_edge(edge, graph);

    reorderID();

}

void Metro::reInsertRemovedStations() {
    cout << "reInsert Removed Stations " << endl;
    vector<EdgeDescriptor> edgesToRemove;
    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        if (graph[edge].stationsOnEdge > 0) {
            edgesToRemove.push_back(edge);
        }
    }

    for (EdgeDescriptor edge : edgesToRemove) {
        VertexDescriptor svd = source(edge, graph);
        VertexDescriptor tvd = target(edge, graph);
        if (graph[svd].id > graph[tvd].id ) {
            svd = tvd;
            tvd = source(edge, graph);
        }
        Coord2 source = *graph[ svd ].coordPtr;
        Coord2 target = *graph[ tvd ].coordPtr;
        VertexDescriptor prevVD = svd;
        for (int i = 1; i <= graph[edge].stationsOnEdge; i++) {
            double k = double(i) / double(graph[edge].stationsOnEdge + 1);
            double x = (1.0 - k) * source.x() + k * target.x();
            double y = (1.0 - k) * source.y() + k * target.y();

            // create new station
            _nStations++;
            VertexDescriptor curVD = add_vertex( graph );
            graph[ curVD ].id               = _nStations;
            graph[ curVD ].coordPtr         = new Coord2( x, y );
            graph[ curVD ].geoPtr           = new Coord2( x, y );
            graph[ curVD ].smoothPtr        = new Coord2( x, y );
            graph[ curVD ].octilinearPtr    = new Coord2( x, y );
            graph[ curVD ].namePtr          = new string( "re-inserted station");
            graph[ curVD ].namePtr          = new string(graph[edge].nameStationsOnEdge[0]);
            
            graph[edge].nameStationsOnEdge.erase(graph[edge].nameStationsOnEdge.begin());
            graph[ curVD ].metroShape       = true;
            graph[ curVD ].inflectionPoint  = false;
            graph[ curVD ].smoothPath       = false;
            graph[ curVD ].interstate       = false;
            graph[ curVD ].weight           = 1.0;
            graph[ curVD ].lineID = graph[ svd ].lineID;

            // connect with prev
            add_edge(curVD, prevVD, graph);
            auto e = boost::edge(curVD, prevVD, graph).first;
            _nEdges++;
            graph[ e ].id          = _nEdges;
            graph[ e ].weight      = 1.0;
            graph[ e ].geoAngle    = graph[edge].geoAngle;
            graph[ e ].smoAngle    = graph[edge].smoAngle;
            graph[ e ].curAngle    = graph[edge].curAngle;
            graph[ e ].lineID      = graph[edge].lineID;
            graph[ e ].isCurved    = false;
            graph[ e ].isMetro     = true;
            graph[ e ].matchPath   = false;
            graph[ e ].isSmoothPath = false;
            graph[ e ].target      = graph[edge].target;
            graph[ e ].stationsOnEdge = 0;

            prevVD = curVD;

            // check if connect with last
            if ( i == graph[edge].stationsOnEdge) {
                add_edge(tvd, curVD, graph);
                e = boost::edge(tvd, curVD, graph).first;
                _nEdges++;
                graph[ e ].id          = _nEdges;
                graph[ e ].weight      = 1.0;
                graph[ e ].geoAngle    = graph[edge].geoAngle;
                graph[ e ].smoAngle    = graph[edge].smoAngle;
                graph[ e ].curAngle    = graph[edge].curAngle;
                graph[ e ].lineID      = graph[edge].lineID;
                graph[ e ].isCurved    = false;
                graph[ e ].isMetro     = true;
                graph[ e ].matchPath   = false;
                graph[ e ].isSmoothPath = false;
                graph[ e ].target      = graph[edge].target;
                graph[ e ].stationsOnEdge = 0;
            }
        }
    }

    for (auto ed : edgesToRemove) {
        graph[ed].isVisible = false;
        remove_edge(ed, graph);
        // _nEdges--;
    }
    reorderID();
}

void Metro::enhanceMetro(double distTol) {
    cout << "Enhance Metro by tolerance " << distTol << endl;
    int addedEdges = 0;
    _nLines++;

    _lineColor[ _nLines ][ 0 ] = 170.0/255.0;
    _lineColor[ _nLines ][ 1 ] = 170.0/255.0;
    _lineColor[ _nLines ][ 2 ] = 170.0/255.0;

    BGL_FORALL_VERTICES(vertI, graph, UndirectedGraph) {
        Coord2 locI = *graph[vertI].geoPtr;
        BGL_FORALL_VERTICES(vertJ, graph, UndirectedGraph) {
            if (vertJ != vertI && !boost::edge(vertI, vertJ, graph).second) {
                Coord2 locJ = *graph[vertJ].geoPtr;
                Coord2 delta = locJ - locI;
                if (delta.norm() < distTol) {
                    add_edge(vertI, vertJ, graph);
                    auto edge = boost::edge(vertI, vertJ, graph).first;
                    addedEdges++;
                    _nEdges++;
                    double angle = atan2( delta.x(), delta.y() );

                    graph[ edge ].id          = _nEdges;
                    graph[ edge ].weight      = 3.0;
                    graph[ edge ].geoAngle    = angle;
                    graph[ edge ].smoAngle    = angle;
                    graph[ edge ].curAngle    = angle;
                    graph[ edge ].lineID.push_back( _nLines );
                    graph[ edge ].isCurved    = false;
                    graph[ edge ].isMetro     = false;
                    graph[ edge ].matchPath   = false;
                    graph[ edge ].isSmoothPath = false;
                }
            }
        }
    }
    cout << "Enhanced Metro by " << addedEdges << " edges" << endl;
}

void Metro::removeAdditionalEdges() {
    cout << "Removing Additional Edges" << endl;
    int removedEdges = 0;
    vector<EdgeDescriptor> edgesToRemove;
    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        if (!graph[edge].isMetro) {
            // remove_edge(edge, graph);
            edgesToRemove.push_back(edge);
            removedEdges++;
            // _nEdges--;
        }
    }
    for (auto edge : edgesToRemove) {
        remove_edge(edge, graph);
        _nEdges--;
    }

    reorderID();
    cout << "Removed Edges " << removedEdges << endl;
}

Coord2 Metro::closestPointOnEdge(EdgeDescriptor edge, Coord2 vertex ) {
    VertexDescriptor svd = source(edge, graph);
    VertexDescriptor tvd = target(edge, graph);
    Coord2 s = *graph[ svd ].coordPtr;
    Coord2 t = *graph[ tvd ].coordPtr;

    Coord2 e = t - s;
    e = e / e.norm();
    Coord2 sv = vertex - s;
    Coord2 p = (sv.x() * e.x() + sv.y() * e.y()) * e;
    p = s + p;

    // check x
    if (p.x() < min(s.x(), t.x())) {
        if (s.x() < t.x()) p = s;
        else p = t;
    } else if (p.x() > max(s.x(), t.x())) {
        if (s.x() > t.x()) p = s;
        else p = t;
    }
    // check y
    else if ( p.y() < min(s.y(), t.y()) ) {
        if ( s.y() < t.y() ) p = s;
        else p = t;
    } else if ( p.y() > max(s.y(), t.y()) ) {
        if ( s.y() > t.y() ) p = s;
        else p = t;
    }

    return p;
}

Coord2 Metro::closestPointOnGraph(VertexDescriptor vertexDis) {
    double minDist = 10e+10;
    Coord2 closestPoint = Coord2(0,0);
    Coord2 vertex = *graph[ vertexDis ].coordPtr;
    unsigned int id = graph[ vertexDis ].id;

    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        VertexDescriptor svd = source(edge, graph);
        VertexDescriptor tvd = target(edge, graph);

        // if ( graph[ svd ].id != id && graph[ tvd ].id != id ) {
        if ( svd != vertexDis && tvd != vertexDis ) {
            Coord2 p = closestPointOnEdge(edge, vertex);
            double d = (vertex - p).norm();
            if (d < minDist) {
                minDist = d;
                closestPoint = p;
            }
        } else {
           //  cout << "same id" << endl;
        }
    }
    // graph[ vertexDis ].closOnMetro = closestPoint;
    graph[ vertexDis ].closOnMetro = closestPoint;
    return closestPoint;
}

Coord2 Metro::closestPointOnGraph(VertexDescriptor vertexDis, Coord2 newCoord) {
    double minDist = 10e+10;
    Coord2 closestPoint = Coord2(0,0);
    Coord2 vertex = newCoord; //*graph[ vertexDis ].coordPtr;
    unsigned int id = graph[ vertexDis ].id;

    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        VertexDescriptor svd = source(edge, graph);
        VertexDescriptor tvd = target(edge, graph);

        if ( svd != vertexDis && tvd != vertexDis ) {
            Coord2 p = closestPointOnEdge(edge, vertex);
            // maybe remove
            if ( ( vertex - *graph[svd].coordPtr ).norm() < ( vertex - *graph[tvd].coordPtr ).norm() )
                p = *graph[svd].coordPtr;
            else
                p = *graph[tvd].coordPtr;
            double d = (vertex - p).norm();
            if (d < minDist) {
                minDist = d;
                closestPoint = p;
            }
        } else {
            //  cout << "same id" << endl;
        }
    }
    // graph[ vertexDis ].closOnMetro = closestPoint;
    graph[ vertexDis ].closOnMetro = closestPoint;
    return closestPoint;
}

void Metro::updateEffectedMetroLines( void )
{
    BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph ) {
        OutEdgeIterator e, e_end;
        graph[vertex].metroShape = false;
        for ( tie( e, e_end ) = out_edges( vertex, graph ); e != e_end; ++e ) {
            EdgeDescriptor ed = *e;
            vector< unsigned int > lineIDs = graph[ ed ].lineID;
            for (unsigned int i = 0; i < lineIDs.size(); i++) {
                unsigned int id = lineIDs[i];
                if (_shapeLines[id]) {
                    graph[vertex].metroShape = true;
                    break;
                }
            }
        }
    }
}

void Metro::scale(Coord2 s) {
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 v = *graph[vertex].coordPtr;
        v.x() *= s.x();
        v.y() *= s.y();
        // graph[vertex].coordPtr = new Coord2(v);
        // graph[vertex].geoPtr = new Coord2(v);

        graph[vertex].coordPtr->setX(v.x());
        graph[vertex].coordPtr->setY(v.y());
        graph[vertex].geoPtr->setX(v.x());
        graph[vertex].geoPtr->setY(v.y());
        graph[vertex].smoothPtr->setX(v.x());
        graph[vertex].smoothPtr->setY(v.y());
    }
    calculateDistBeta();

}

void Metro::translate(Coord2 trans) {
    cerr << "TRANSLATE: " << trans << endl;
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 v = *graph[vertex].coordPtr;
        graph[vertex].coordPtr->setX(v.x() + trans.x());
        graph[vertex].coordPtr->setY(v.y() + trans.y());
        graph[vertex].geoPtr->setX(v.x() + trans.x());
        graph[vertex].geoPtr->setY(v.y() + trans.y());
        graph[vertex].smoothPtr->setX(v.x() + trans.x());
        graph[vertex].smoothPtr->setY(v.y() + trans.y());
        // graph[vertex].geoPtr = new Coord2(0,0);
    }

}

double Metro::magnitude( Coord2 v )
{
    return sqrt(v.x() * v.x() + v.y() * v.y());
}

//
//  Metro::~checkVEConflicrs --    detect vertex-edge pair that is close to each other
//
//  Inputs
//  none
//
//  Outputs
//  none
//
bool Metro::checkVEConflicts( void )
{
    bool                doConflict      = true;
    _VEconflict.clear();
    // find the vertex-edge pair that has smaller distance than the threshold
    BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
    {
        double              distThreshold   = MIN_DISTANCE;

        Coord2 coord = *graph[vertex].coordPtr;//*graph[vertex].smoothPtr;
        unsigned idV = graph[vertex].id;

        double minDist = 1000000;
        BGL_FORALL_EDGES( edge, graph, UndirectedGraph )
        {
            unsigned idE = graph[edge].id;
            VertexDescriptor vdS = source( edge, graph );
            VertexDescriptor vdT = target( edge, graph );

            if ( ( vertex != vdS ) && ( vertex != vdT ) ) {
                Coord2 coordS = *graph[vdS].coordPtr;
                Coord2 coordT = *graph[vdT].coordPtr;

                double m = ( coordS.y() - coordT.y() ) / ( coordS.x() - coordT.x() );
                double k = coordS.y() - m * coordS.x();
                double dist = fabs( m * coord.x() - coord.y() + k ) / sqrt( SQUARE( m ) + 1.0 );
                Coord2 vec, p;
                vec.x() = ( coordS - coordT ).y();
                vec.y() = -( coordS - coordT ).x();
                vec = vec / vec.norm();
                p = coord + dist * vec;
                if( fabs( p.y() - m * p.x() - k ) > 1.0e-8 ){
                    p = coord - dist * vec;
                }
                double r = ( p - coordT ).x() / ( coordS - coordT ).x() ;
                if (!graph[vdT].isStation or !graph[vdS].isStation)
                    distThreshold *= 0.5;
                 if( ( dist < distThreshold ) && ( 0.0 <= r ) && ( r <= 1.0 ) ) {
                    _VEconflict.insert( pair< Grid2, VEPair >( Grid2( idV, idE ), VEPair( vertex, edge ) ) );
                 }
            }
        }
    }
//
//    // update ratioR
    _ratioR.clear();
    _ratioR.resize( _VEconflict.size() );
    _ratioGeo.clear();
    _ratioGeo.resize(_VEconflict.size() );
    for (int i = 0; i < _ratioR.size(); i++ ) {
        _ratioR[i] = 0.0;
        _ratioGeo[i] = 0.0;
    }
    unsigned int countVE = 0;
    // for ( VEMap::iterator it = _VEconflict.begin();
    //        it != _VEconflict.end(); ++it ) {
    for (auto it : _VEconflict) {
        VertexDescriptor vdV = it.second.first;
        EdgeDescriptor ed = it.second.second;
        VertexDescriptor vdS = source( ed, graph );
        VertexDescriptor vdT = target( ed, graph );
        Coord2 pointV = *graph[vdV].coordPtr;
        Coord2 pointS = *graph[vdS].coordPtr;
        Coord2 pointT = *graph[vdT].coordPtr;

        double pointSPointT = pointS.x() - pointT.x();
        if (fabs(pointSPointT) < 0.0000001){
            cerr << "Devide by zero" << endl;
            if (pointSPointT < 0.0) {
                pointSPointT -= 0.001;
            } else {
                pointSPointT += 0.001;
            }
        }


        double m = ( pointS.y() - pointT.y() ) /  pointSPointT;
        double k = pointS.y() - m * pointS.x();
        double dist = fabs( m * pointV.x() - pointV.y() + k ) / sqrt( SQUARE( m ) + 1.0 );
        Coord2 vec, p;
        vec.x() = ( pointS - pointT ).y();
        vec.y() = -( pointS - pointT ).x();
        vec = vec / vec.norm();
        p = pointV + dist * vec;
        if( fabs( p.y() - m * p.x() - k ) > 1.0e-8 ){
            p = pointV - dist * vec;
        }
        double r;
        if (( pointS - pointT ).x() > 0.00001) {
            r = ( p - pointT ).x() / ( pointS - pointT ).x() ;
        } else {
            r = 1.;
        }

        _ratioR[ countVE ] = r;

        Coord2 pointV_ = *graph[vdV].geoPtr;
        Coord2 pointS_ = *graph[vdS].geoPtr;
        Coord2 pointT_ = *graph[vdT].geoPtr;
        double poitS_PointT = ( pointS_.x() - pointT_.x() );
        if (fabs(poitS_PointT) < 0.0000001) {
            cerr << "Devide by zero" << endl;
            if (poitS_PointT < 0) {
                poitS_PointT -= 0.0001;
            } else {
                poitS_PointT += 0.0001;
            }
        }

        double m_ = ( pointS_.y() - pointT_.y() ) / ( pointS_.x() - pointT_.x() );

        double k_ = pointS_.y() - m_ * pointS_.x();
        double dist_ = fabs( m_ * pointV_.x() - pointV_.y() + k_ ) / sqrt( SQUARE( m_ ) + 1.0 );
        Coord2 vec_, p_;
        vec_.x() = ( pointS_ - pointT_ ).y();
        vec_.y() = -( pointS_ - pointT_ ).x();
        vec_ = vec_ / vec_.norm();
        p_ = pointV_ + dist_ * vec_;
        if( fabs( p_.y() - m_ * p_.x() - k_ ) > 1.0e-8 ){
            p_ = pointV_ - dist_ * vec_;
        }
        double r_;
        if (( pointS_ - pointT_ ).x() > 0.00001 ) {
            r_ = ( p_ - pointT_ ).x() / ( pointS_ - pointT_ ).x() ;
            r_ = ( p_ - pointT_ ).x() / ( pointS_ - pointT_ ).x() ;
        } else {
            r_ = 1.;
        }
        _ratioGeo[ countVE ] = r_;
        countVE++;
    }

    return doConflict;
}


//
//  Metro::adjustsize --  adjust size of the layout of the metro network
//
//  Inputs
//      width: window width
//      height: window height
//
//  Outputs
//      none
//
void Metro::adjustsize( const int & width, const int & height )
{
//    cerr << "width = " << width << " height = " << height << endl;
    double xMin =  INFINITY;
    double xMax = -INFINITY;
    double yMin =  INFINITY;
    double yMax = -INFINITY;
    double aspect = ( double )width/( double )height;

    // Scan all the vertex coordinates first
    BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
    {
        Coord2 & coord = *graph[ vertex ].coordPtr;
        if ( coord.x() < xMin ) xMin = coord.x();
        if ( coord.x() > xMax ) xMax = coord.x();
        if ( coord.y() < yMin ) yMin = coord.y();
        if ( coord.y() > yMax ) yMax = coord.y();
    }

    // double range = 0.5 * MAX2( xMax - xMin, yMax - yMin );
    double xRange;
    double yRange;
    double xMid;
    double yMid;
    if( ( xMax - xMin ) / width > ( yMax - yMin ) / height ) {
        xRange  = 1.0 * ( xMax - xMin );
        yRange  = 1.0 * ( xMax - xMin ) * ( 1.0/ aspect );
    }
    else {
        xRange  = 1.0 * ( yMax - yMin ) * aspect;
        yRange  = 1.0 * ( yMax - yMin );
    }

    xRange *= 1.05;
    yRange *= 1.05;
    xMid    = 0.5 * ( xMin + xMax );
    yMid    = 0.5 * ( yMin + yMax );

    // Normalize the coordinates
    BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
    {
        Coord2 geo          = * graph[ vertex ].geoPtr;
        Coord2 smooth       = * graph[ vertex ].smoothPtr;
        Coord2 octilinear   = * graph[ vertex ].octilinearPtr;
	    Coord2 coord        = * graph[ vertex ].coordPtr;
//        Coord2 label = vertexExternal[ vertex ].curSite();

        geo.setX( width  * ( geo.x() - xMid ) / xRange );
        geo.setY( height * ( geo.y() - yMid ) / yRange );
        smooth.setX( width  * ( smooth.x() - xMid ) / xRange );
        smooth.setY( height * ( smooth.y() - yMid ) / yRange );
        coord.setX( width  * ( coord.x() - xMid ) / xRange );
        coord.setY( height * ( coord.y() - yMid ) / yRange );
//        label.setX( width  * ( label.x() - xMid ) / xRange );
//        label.setY( height * ( label.y() - yMid ) / yRange );

        graph[ vertex ].geoPtr->x() = geo.x();
        graph[ vertex ].geoPtr->y() = geo.y();
        graph[ vertex ].smoothPtr->x() = smooth.x();
        graph[ vertex ].smoothPtr->y() = smooth.y();
	    graph[ vertex ].octilinearPtr->x() = octilinear.x();
	    graph[ vertex ].octilinearPtr->y() = octilinear.y();
        graph[ vertex ].coordPtr->x() = coord.x();
        graph[ vertex ].coordPtr->y() = coord.y();
//        vertexExternal[ vertex ].geoSite().x() = label.x();
//        vertexExternal[ vertex ].curSite().x() = label.x();
//        vertexExternal[ vertex ].geoSite().y() = label.y();
//        vertexExternal[ vertex ].curSite().y() = label.y();

    }
    calculateDistBeta();

}


void Metro::calculateDistBeta() {
    // compute the unit length of an edge (ratio)
    int nAlpha = 0;
    int nBeta = 0;
    double totallength = 0.0;
    BGL_FORALL_EDGES( edge, graph, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, graph );
        VertexDescriptor vdT = target( edge, graph );
        unsigned int idS = graph[ vdS ].id;
        unsigned int idT = graph[ vdT ].id;

        Coord2 coord = *graph[ vdT ].geoPtr - *graph[ vdS ].geoPtr;
        totallength += coord.norm();
        double w = graph[ edge ].weight;
        if( w == 1.0 ) nBeta++;
        else nAlpha++;
    }
    double magLength = 2.0; // for Tokyo 1.0
    cerr << "calculating _distBeta " << endl;
    _distBeta = totallength / ( magLength * nAlpha + nBeta );
    _distBeta *= 1.2;
    _distAlpha = magLength * _distBeta;

    // cout << "dist beta old: " << _distBeta << endl;
    // cout << "dist alpha old: " << _distAlpha << endl;

    // Take the mean instead
    // vector<double> distances;
    // BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
    //         VertexDescriptor vdS = source( edge, graph );
    //         VertexDescriptor vdT = target( edge, graph );
    //         unsigned int idS = graph[ vdS ].id;
    //         unsigned int idT = graph[ vdT ].id;
//
    //         Coord2 coord = *graph[ vdT ].geoPtr - *graph[ vdS ].geoPtr;
    //         distances.push_back(coord.norm());
    // }
    // // distances.sort();
    // sort(distances.begin(), distances.end());
    // // _distBeta = _distAlpha = distances[(int) (distances.size() * .5) ];
    // cout << "dist beta new: " << _distBeta << endl;
}

//
//  Metro::simplifyLayout --    delete nodes with small distance to the edges bounded by their adjacent nodes
//
//  Inputs
//  none
//
//  Outputs
//  none
//
void Metro::simplifyLayout( void )
{
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    EdgeIndexMap        edgeIndex       = get( edge_index, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeGeoAngleMap     edgeGeoAngle    = get( edge_mygeoangle, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeLineIDMap       edgeLineID      = get( edge_mylineid, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
//
//    //double threshold = 2.5;
//    double threshold = 10.0;
//    int count = 0;
//    _removedVertices.clear();
//    _removedEdges.clear();
//    _removedWeights.clear();
//    //printGraph( graph );
//
//    while( true ){
//        // find minimal distance of middle node to the line composed by its neighbors
//        VertexDescriptor mvd, mvdP, mvdN;
//        bool noDegreeTwo = true;
//        double minDist = INFINITY;
//        BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
//        {
//            DegreeSizeType degrees = out_degree( vertex, graph );
//
//            if( degrees == 2 && ( vertexIsStation[ vertex ] == true ) ){
//
//                noDegreeTwo = false;
//                OutEdgeIterator e, e_end;
//                tie( e, e_end ) = out_edges( vertex, graph );
//                EdgeDescriptor edP = *e;
//                VertexDescriptor vTP = target( edP, graph );
//                e++;
//                EdgeDescriptor edN = *e;
//                VertexDescriptor vTN = target( edN, graph );
//
//                Coord2 p = vertexCoord[ vTP ];
//                Coord2 c = vertexCoord[ vertex ];
//                Coord2 n = vertexCoord[ vTN ];
//                double m = ( p.y() - n.y() ) / ( p.x() - n.x() );
//                double k = p.y() - m * p.x();
//                double dist = fabs( m * c.x() - c.y() + k ) / sqrt( SQUARE( m ) + 1.0 );
//
//#ifdef  DEBUG
//                cerr << "V(" << vertexIndex[ vertex ] << "), d = " << degrees
//                     << ", dist =  " << dist << endl;
//                cerr << "p = " << p;
//                cerr << "c = " << c;
//                cerr << "n = " << n;
//#endif  // DEBUG
//                if( dist < minDist ) {
//                    minDist = dist;
//                    mvd = vertex;
//                    mvdP = vTP;
//                    mvdN = vTN;
//                }
//            }
//        }
//
//        // check condition
//        // cerr << "count = " << count << ", minID = " << vertexIndex[ mvd ] << ", minDist = " << minDist << endl;
//        // if( minDist > threshold || noDegreeTwo == true || vertexIsStation[ mvd ] == false ) break;
//        if( minDist > threshold || noDegreeTwo == true ) break;
//
//        // delete vertex and corresponding neighbor edges
//        EdgeDescriptor      edP, edN;
//        bool                found = false;
//        tie( edP, found ) = edge( mvd, mvdP, graph );
//        assert( found == true );
//        tie( edN, found ) = edge( mvd, mvdN, graph );
//        assert( found == true );
//        double weight = edgeWeight[ edP ] + edgeWeight[ edN ];
//
//        // remove stations and edges
//        remove_edge( edP, graph );
//        _removedEdges.push_back( VVIDPair( vertexIndex[ mvd ], vertexIndex[ mvdP ] ) );
//        _removedWeights.push_back( edgeWeight[ edP ] );
//        remove_edge( edN, graph );
//        _removedEdges.push_back( VVIDPair( vertexIndex[ mvd ], vertexIndex[ mvdN ] ) );
//        _removedWeights.push_back( edgeWeight[ edN ] );
//        clear_vertex( mvd, graph );
//        remove_vertex( mvd, graph );
//        _removedVertices.push_back( vertexIndex[ mvd ] );
//
//        // cerr << "mv = " << vertexID[ mvd ] << " mvdP = " << vertexID[ mvdP ] << " mvdN = " << vertexID[ mvdN ] << endl;
//        // cerr << "wS = " << edgeWeight[ edP ] << " wT = " << edgeWeight[ edN ] << endl;
//
//        // add new edges
//        pair< EdgeDescriptor, unsigned int > addE = add_edge( mvdP, mvdN, 1, graph );
//        EdgeDescriptor addED = addE.first;
//
//        // calculate geographical angle
//        Coord2 coordO;
//        Coord2 coordD;
//        if( vertexIndex[ mvdP ] < vertexIndex[ mvdN ] ){
//            coordO = vertexCoord[ mvdP ];
//            coordD = vertexCoord[ mvdN ];
//        }
//        else{
//            coordO = vertexCoord[ mvdN ];
//            coordD = vertexCoord[ mvdP ];
//        }
//        double diffX = coordD.x() - coordO.x();
//        double diffY = coordD.y() - coordO.y();
//        double angle = atan2( diffY, diffX );
//
//        edgeWeight[ addED ] = weight;
//        edgeGeoAngle[ addED ] = angle;
//        edgeSmoAngle[ addED ] = angle;
//        edgeCurAngle[ addED ] = angle;
//        edgeIsLine[ addED ] = true;
//#ifdef  DEBUG
//        cerr << "addWeight = " << weight <<endl;
//        cerr << "count = " << count << ", edgeWeight[ addED ] = " << edgeWeight[ addED ] << endl;
//#endif  // DEBUG
//
//        reorderID();
//        count++;
//    }
//
//    assert( ( 2*_removedVertices.size() ) == _removedEdges.size() );
#ifdef  DEBUG
    for( unsigned int i = 0; i < _removedVertices.size(); i++ ){
        cerr << "removedV(" << _removedVertices[ i ] << "), ";
    }
    cerr << endl;
    for( unsigned int i = 0; i < _removedEdges.size(); i++ ){
        cerr << "removedE(" << _removedEdges[ i ].first << ", "
             << _removedEdges[ i ].second << "), ";
    }
    cerr << endl;
    for( unsigned int i = 0; i < _removedWeights.size(); i++ ){
        cerr << "removedW(" << _removedWeights[ i ] << "), ";
    }
    cerr << endl;
#endif  // DEBUG
}

//
//  Metro::reorderID --    reorder vertex and edge ID
//
//  Inputs
//  none
//
//  Outputs
//  none
//
void Metro::reorderID( void )
{
    _nStations = 0;
   BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
   {
        graph[vertex].id = _nStations;
       _nStations++;
   }
    _nEdges = 0;
    BGL_FORALL_EDGES( edge, graph, UndirectedGraph )
    {
        graph[edge].id = _nEdges;
        _nEdges++;
    }
}

//
//  Metro::movebackSmooth --    move the deleted nodes back to the layout
//
//  Inputs
//  obj: the original layout
//
//  Outputs
//  none
//
bool Metro::movebackSmooth( const Metro & obj )
{
//    if( ( _removedVertices.size() == 0 ) && ( _removedEdges.size() == 0 ) )
//        return true;
//
//    UndirectedGraph oGraph = obj.g();
//
//    // simplified graph
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexHomeMap       vertexHome      = get( vertex_myhome, graph );
//    VertexGeoMap        vertexGeo       = get( vertex_mygeo, graph );
//    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeGeoAngleMap     edgeGeoAngle    = get( edge_mygeoangle, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
//
//    // original graph
//    VertexHomeMap       oVertexHome     = get( vertex_myhome, oGraph );
//    VertexGeoMap        oVertexGeo      = get( vertex_mygeo, oGraph );
//
//    //cerr << "_removedVerteces = " << _removedVertices.size() << endl;
//    //cerr << "_removedEdges = " << _removedEdges.size() << endl;
//
//    unsigned int idV = _removedVertices.back();
//    _removedVertices.pop_back();
//    VVIDPair pair = _removedEdges.back();
//    _removedEdges.pop_back();
//    unsigned int idT = pair.second;
//    double wT = _removedWeights.back();
//    _removedWeights.pop_back();
//    pair = _removedEdges.back();
//    _removedEdges.pop_back();
//    unsigned int idS = pair.second;
//    double wS = _removedWeights.back();
//    _removedWeights.pop_back();
//
//
//    // find VertexDescriptor in simplified graph
//    VertexDescriptor vdS, vdT;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//        // cerr << "vertexIndex[ vd ] = " << vertexIndex[ vd ] << endl;
//        if ( vertexIndex[ vd ] == idS ) vdS = vd;
//        if ( vertexIndex[ vd ] == idT ) vdT = vd;
//    }
//
//    bool found = false;
//    EdgeDescriptor ed;
//    tie( ed, found ) = edge( vdS, vdT, graph );
//    // remove edge
//    remove_edge( ed, graph );
//    // add vertex
//    VertexDescriptor vdNew = add_vertex( graph );
//    vertexIndex[ vdNew ] = idV;
//    VertexDescriptor vdO = vertex( idV, oGraph );
//    Coord2 coord = ( wT * vertexSmooth[ vdS ] + wS * vertexSmooth[ vdT ] ) / ( wS + wT );
//    vertexHome[ vdNew ] = oVertexHome[ vdO ];
//    vertexGeo[ vdNew ] = oVertexGeo[ vdO ];
//    vertexSmooth[ vdNew ] = coord;
//    vertexCoord[ vdNew ] = coord;
//    vertexIsStation[ vdNew ] = true;
//    // add edge
//    ed = add_edge( vdNew, vdS, 1, graph ).first;
//    edgeWeight[ ed ] = wS;
//    edgeIsLine[ ed ] = true;
//    // add  edge
//    ed = add_edge( vdNew, vdT, 1, graph ).first;
//    edgeWeight[ ed ] = wT;
//    edgeIsLine[ ed ] = true;
//
//    // reorder vertexID
//    reorderID();
//
//    // update edge angle of simplified metro layout
//    BGL_FORALL_EDGES( edge, graph, UndirectedGraph ){
//
//        VertexDescriptor vdS = source( edge, graph );
//        VertexDescriptor vdT = target( edge, graph );
//
//        // calculate geographical angle
//        Coord2 coordO;
//        Coord2 coordD;
//        if( vertexID[ vdS ] < vertexID[ vdT ] ){
//            coordO = vertexCoord[ vdS ];
//            coordD = vertexGeo[ vdT ];
//        }
//        else{
//            coordO = vertexGeo[ vdT ];
//            coordD = vertexGeo[ vdS ];
//        }
//        double diffX = coordD.x() - coordO.x();
//        double diffY = coordD.y() - coordO.y();
//        double angle = atan2( diffY, diffX );
//
//        edgeGeoAngle[ edge ] = angle;
//
//        if( vertexID[ vdS ] < vertexID[ vdT ] ){
//            coordO = vertexCoord[ vdS ];
//            coordD = vertexCoord[ vdT ];
//        }
//        else{
//            coordO = vertexCoord[ vdT ];
//            coordD = vertexCoord[ vdS ];
//        }
//        diffX = coordD.x() - coordO.x();
//        diffY = coordD.y() - coordO.y();
//        angle = atan2( diffY, diffX );
//
//        edgeSmoAngle[ edge ] = angle;
//        edgeCurAngle[ edge ] = angle;
//
//    }

    return false;
}

//
//  Metro::movebackOctilinear --    move the deleted nodes back to the layout
//
//  Inputs
//  obj: the original layout
//
//  Outputs
//  none
//
bool Metro::movebackOctilinear( const Metro & obj )
{
//    if( ( _removedVertices.size() == 0 ) && ( _removedEdges.size() == 0 ) )
//        return true;
//
//    UndirectedGraph oGraph = obj.g();
//
//    // simplified graph
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexHomeMap       vertexHome      = get( vertex_myhome, graph );
//    VertexGeoMap        vertexGeo       = get( vertex_mygeo, graph );
//    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
//
//    // original graph
//    VertexHomeMap       oVertexHome     = get( vertex_myhome, oGraph );
//    VertexGeoMap        oVertexGeo      = get( vertex_mygeo, oGraph );
//    VertexSmoothMap     oVertexSmooth   = get( vertex_mysmooth, oGraph );
//
//    //cerr << "_removedVerteces = " << _removedVertices.size() << endl;
//    //cerr << "_removedEdges = " << _removedEdges.size() << endl;
//
//    unsigned int idV = _removedVertices.back();
//    _removedVertices.pop_back();
//    VVIDPair pair = _removedEdges.back();
//    _removedEdges.pop_back();
//    unsigned int idT = pair.second;
//    double wT = _removedWeights.back();
//    _removedWeights.pop_back();
//    pair = _removedEdges.back();
//    _removedEdges.pop_back();
//    unsigned int idS = pair.second;
//    double wS = _removedWeights.back();
//    _removedWeights.pop_back();
//
//
//    // find VertexDescriptor in simplified graph
//    VertexDescriptor vdS, vdT;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//        // cerr << "vertexIndex[ vd ] = " << vertexIndex[ vd ] << endl;
//        if ( vertexIndex[ vd ] == idS ) vdS = vd;
//        if ( vertexIndex[ vd ] == idT ) vdT = vd;
//    }
//
//    bool found = false;
//    EdgeDescriptor ed;
//    tie( ed, found ) = edge( vdS, vdT, graph );
//    // remove edge
//    remove_edge( ed, graph );
//    // add vertex
//    VertexDescriptor vdNew = add_vertex( graph );
//    vertexIndex[ vdNew ] = idV;
//    VertexDescriptor vdO = vertex( idV, oGraph );
//    Coord2 coord = ( wT * vertexCoord[ vdS ] + wS * vertexCoord[ vdT ] ) / ( wS + wT );
//    vertexHome[ vdNew ] = oVertexHome[ vdO ];
//    vertexGeo[ vdNew ] = oVertexGeo[ vdO ];
//    vertexSmooth[ vdNew ] = oVertexSmooth[ vdO ];
//    vertexCoord[ vdNew ] = coord;
//    vertexIsStation[ vdNew ] = true;
//    // add edge
//    ed = add_edge( vdNew, vdS, 1, graph ).first;
//    edgeWeight[ ed ] = wS;
//    edgeIsLine[ ed ] = true;
//    // add  edge
//    ed = add_edge( vdNew, vdT, 1, graph ).first;
//    edgeWeight[ ed ] = wT;
//    edgeIsLine[ ed ] = true;
//
//    // reorder vertexID
//    reorderID();
//
//    // update edge angle of simplified metro layout
//    BGL_FORALL_EDGES( edge, graph, UndirectedGraph ){
//
//        VertexDescriptor vdS = source( edge, graph );
//        VertexDescriptor vdT = target( edge, graph );
//
//        // calculate geographical angle
//        Coord2 coordO;
//        Coord2 coordD;
//        if( vertexID[ vdS ] < vertexID[ vdT ] ){
//            coordO = vertexCoord[ vdS ];
//            coordD = vertexCoord[ vdT ];
//        }
//        else{
//            coordO = vertexCoord[ vdT ];
//            coordD = vertexCoord[ vdS ];
//        }
//        double diffX = coordD.x() - coordO.x();
//        double diffY = coordD.y() - coordO.y();
//        double angle = atan2( diffY, diffX );
//
//        //edgeSmoAngle[ edge ] = angle;
//        edgeCurAngle[ edge ] = angle;
//
//    }

    return false;
}

//
//  Metro::addAdditionalEdges --    add edges if the label node is too close to a station
//
//  Inputs
//  none
//
//  Outputs
//  none
//
void Metro::addAdditionalEdges( void )
{
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexNameMap       vertexName      = get( vertex_myname, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    VertexExternalMap   vertexExternal  = get( vertex_myexternal, graph );
//    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, graph );
//    EdgeIndexMap        edgeIndex       = get( edge_index, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeGeoAngleMap     edgeGeoAngle    = get( edge_mygeoangle, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
/*
    // if the stations and labels are too close to each other, then add an edge
    BGL_FORALL_VERTICES( vdO, graph, UndirectedGraph ){

        if( vertexIsStation[ vdO ] == false ){

            Coord2 coordO = vertexCoord[ vdO ];
            BGL_FORALL_VERTICES( vdI, graph, UndirectedGraph ){

                Coord2 coordI = vertexCoord[ vdI ];
                Coord2 diff = coordO - coordI;
                double leaderW = 0.0;

                OutEdgeIterator e, e_end;
                for ( tie( e, e_end ) = out_edges( vdO, graph ); e != e_end; ++e ) {
                    EdgeDescriptor ed = *e;
                    VertexDescriptor vS = source( ed, graph );
                    VertexDescriptor vT = target( ed, graph );
                    if( vertexName[ vS ] == vertexName[ vT ] ) {
                        //leaderW = edgeWeight[ ed ];
                        leaderW = ( vertexCoord[ vS ] - vertexCoord[ vT ] ).norm()/_distBeta;
#ifdef  DEBUG
                        cerr << "id = " << vertexID[ vS ] << ", " << vertexID[ vT ]
                             << ", w = " << leaderW << endl;
#endif  // DEBUG
                    }
                }

                if( ( diff.norm() < _distBeta * leaderW ) && vertexID[ vdO ] != vertexID[ vdI ] ) {
                    vertexSelectMag[ vdI ] = true;
                }
            }
        }
    }
*/
}

//
//  Metro::cloneLayout --    copy all information of the input metro
//
//  Inputs
//  obj: original metro layout
//
//  Outputs
//  none
//
void Metro::cloneLayout( const Metro & obj )
{
    graph.clear();
    graph = obj.g();

    _nLines = obj._nLines;
    _nStations = obj._nStations;
    _nEdges = obj._nEdges;
    _nLabels = obj._nLabels;
    _distAlpha = obj._distAlpha;
    _distBeta = obj._distBeta;
    _meanVSize = obj._meanVSize;

    for ( unsigned int k = 0; k < _nLines; ++k ) {

        // copy the line color
        for ( unsigned i = 0; i < 3; ++i )
            _lineColor[ k ][ i ] = obj._lineColor[ k ][ i ];

        // copy the line name
        strcpy( _lineName[ k ], obj._lineName[ k ] );

        // copy line info
        vector< EdgeDescriptor > line;
        for ( unsigned i = 0; i < obj.line().size(); ++i ){
            line.push_back( obj.line()[ k ][ i ] );
        }
        _line.push_back( line );

        // copy line station info
        vector< VertexDescriptor > lineSta;
        for ( unsigned i = 0; i < obj.lineSta().size(); ++i ){
            lineSta.push_back( obj.lineSta()[ k ][ i ] );
        }
        _lineSta.push_back( lineSta );
    }
}

//
//  Metro::cloneLabel --    copy label node to the graph
//
//  Inputs
//  obj: original metro layout
//
//  Outputs
//  none
//
void Metro::cloneLabel( void )
{
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexHomeMap       vertexHome      = get( vertex_myhome, graph );
//    VertexGeoMap        vertexGeo       = get( vertex_mygeo, graph );
//    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, graph );
//    VertexExternalMap   vertexExternal  = get( vertex_myexternal, graph );
//    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    VertexNameMap       vertexName      = get( vertex_myname, graph );
//    VertexTexNameMap    vertexTexName   = get( vertex_mytexname, graph );
//    VertexScaleMap      vertexScale     = get( vertex_myscale, graph );
//    VertexSizeMap       vertexSize      = get( vertex_mysize, graph );
//    EdgeIndexMap        edgeIndex       = get( edge_index, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeGeoAngleMap     edgeGeoAngle    = get( edge_mygeoangle, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeLineIDMap       edgeLineID      = get( edge_mylineid, graph );
//    EdgeSelectShiftMap  edgeSelectShift = get( edge_myselectshift, graph );
//    EdgeSelectCtrlMap   edgeSelectCtrl  = get( edge_myselectctrl, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
//    EdgeIsLeaderMap     edgeIsLeader    = get( edge_myisleader, graph );
//
//    _nLabels = 0;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//
//        // add labels
//        if( ( vertexIsStation[ vd ] == true ) && ( vertexExtstate[ vd ] == true ) &&
//            ( vertexSelectMag[ vd ] == true ) ){
//
//            Coord2 geoSite = vertexExternal[ vd ].geoSite();
//            Coord2 curSite = vertexExternal[ vd ].curSite();
//            //cerr << "VID = " << vertexID[ vd ] << " eW = " << vertexExternal[ vd ].leaderWeight() << endl;
//            //cerr << "id = " << vertexID[ vd ] << " " << vertexName[ vd ] << endl;
//            //cerr << "coord = " << vertexCoord[ vd ];
//            //cerr << "geoSite = " << geoSite;
//            //cerr << "curSite = " << curSite;
//            //cerr << "VID = " << vertexID[ vd ] << " eW = " << vertexExternal[ vd ].leaderWeight() << endl;
//
//            // add label node
//            VertexDescriptor curVD =  add_vertex( graph );
//
//            vertexCoord[ curVD ]            = curSite;
//            vertexSmooth[ curVD ]           = curSite;
//            vertexGeo[ curVD ]              = geoSite;
//            vertexHome[ curVD ]             = geoSite;
//            vertexIndex[ curVD ]            = _nStations;
//            vertexID[ curVD ]               = _nStations;
//            vertexName[ curVD ]             = vertexName[ vd ];
//            vertexTexName[ curVD ]          = vertexTexName[ vd ];
//            vertexScale[ curVD ]            = vertexScale[ vd ];
//            vertexSelectMag[ curVD ]        = true;
//            vertexIsStation[ curVD ]        = false;
//            vertexSize[ curVD ]             = vertexSize[ vd ];
//            vertexExternal[ curVD ].leaderWeight() = vertexExternal[ vd ].leaderWeight();
//
//            _nStations++;
//
//            // add label edge
//            pair<EdgeDescriptor, unsigned int> foreE = add_edge( vd, curVD,
//                                                                 1.0, graph );
//            EdgeDescriptor foreED = foreE.first;
//            Coord2 coordO = vertexCoord[ vd ];
//            Coord2 coordD = vertexCoord[ curVD ];
//            double diffX = coordD.x() - coordO.x();
//            double diffY = coordD.y() - coordO.y();
//            double angle = atan2( diffY, diffX );
//
//            edgeIndex[ foreED ]             = _nEdges;
//            edgeID[ foreED ]                = _nEdges;
//            edgeGeoAngle[ foreED ]          = angle;
//            edgeSmoAngle[ foreED ]          = angle;
//            edgeCurAngle[ foreED ]          = angle;
//            edgeIsLine[ foreED ]            = false;
//            edgeIsLeader[ foreED ]          = true;
//            edgeWeight[ foreED ]            = vertexExternal[ vd ].leaderWeight();
//            // ( vertexCoord[ vd ] - vertexCoord[ curVD ] ).norm()/_distBeta;
//            _nEdges++;
//            _nLabels++;
//        }
//
//    }
//
//    _nSbeforeSim = _nStations;
//    _nEbeforeSim = _nEdges;
#ifdef  DEBUG
    cout << " numStations = " << _nStations << " num_vertices = " << num_vertices( graph ) << endl;
    cout << " numEdges = " << _nEdges << " num_edges = " << num_edges( graph ) << endl;
    cout << " nLabels = " << _nLabels << endl;
#endif  // DEBUG
}

//
//  Metro::cloneSubLabel --    copy label node to the graph
//
//  Inputs
//  obj: original metro layout
//
//  Outputs
//  none
//
void Metro::cloneSubLabel( void )
{
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexIDMap         vertexID        = get( vertex_myid, graph );
//    VertexHomeMap       vertexHome      = get( vertex_myhome, graph );
//    VertexGeoMap        vertexGeo       = get( vertex_mygeo, graph );
//    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, graph );
//    VertexExternalMap   vertexExternal  = get( vertex_myexternal, graph );
//    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//    VertexNameMap       vertexName      = get( vertex_myname, graph );
//    VertexTexNameMap    vertexTexName   = get( vertex_mytexname, graph );
//    VertexScaleMap      vertexScale     = get( vertex_myscale, graph );
//    EdgeIndexMap        edgeIndex       = get( edge_index, graph );
//    EdgeIDMap           edgeID          = get( edge_myid, graph );
//    EdgeWeightMap       edgeWeight      = get( edge_weight, graph );
//    EdgeGeoAngleMap     edgeGeoAngle    = get( edge_mygeoangle, graph );
//    EdgeSmoAngleMap     edgeSmoAngle    = get( edge_mysmoangle, graph );
//    EdgeCurAngleMap     edgeCurAngle    = get( edge_mycurangle, graph );
//    EdgeLineIDMap       edgeLineID      = get( edge_mylineid, graph );
//    EdgeSelectShiftMap  edgeSelectShift = get( edge_myselectshift, graph );
//    EdgeSelectCtrlMap   edgeSelectCtrl  = get( edge_myselectctrl, graph );
//    EdgeIsLineMap       edgeIsLine      = get( edge_myisline, graph );
//
//    // add sub labels
//    double w = 15.0;
//    vector< Coord2 > corner( 3 );
//
//    unsigned int count = 0;
//    unsigned int nL = _nLabels;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ) {
//
//        if( vertexIsStation[ vd ] == false && count < nL ) {
//
//            // cerr << " count = " << count << " nLabels = " << _nLabels << endl;
//            VertexDescriptor preVD = vd;
//
//            OutEdgeIterator e, e_end;
//            tie( e, e_end ) = out_edges( vd, graph );
//            VertexDescriptor staVD = target( *e, graph );
//            double vSign = 1.0, hSign = 1.0;
//            Coord2 vec = vertexCoord[ vd ] - vertexCoord[ staVD ];
//            cerr << vec;
//            if( vec.x() > 0.0 && vec.y() > 0.0 ) {
//                hSign = 1.0;
//                vSign = 1.0;
//            }
//            else if( vec.x() > 0.0 && vec.y() < 0.0 ) {
//                hSign = 1.0;
//                vSign = -1.0;
//            }
//            else if( vec.x() < 0.0 && vec.y() < 0.0 ) {
//                hSign = -1.0;
//                vSign = -1.0;
//            }
//            else {
//                hSign = -1.0;
//                vSign = 1.0;
//            }
//            corner[ 0 ].x() = 0.0;
//            corner[ 0 ].y() = vSign * w;
//            corner[ 1 ].x() = hSign * w;
//            corner[ 1 ].y() = vSign * w;
//            corner[ 2 ].x() = hSign * w;
//            corner[ 2 ].y() = 0.0;
//
//            Coord2 site = vertexCoord[ vd ];
//            unsigned int size = corner.size();
//            for( unsigned int i = 0; i < size; i++ ){
//
//                // add sub label node
//                VertexDescriptor lVD =  add_vertex( graph );
//
//                vertexCoord[ lVD ]            = site + corner[ i ];
//                vertexSmooth[ lVD ]           = site + corner[ i ];
//                vertexGeo[ lVD ]              = site + corner[ i ];
//                vertexHome[ lVD ]             = site + corner[ i ];
//                vertexIndex[ lVD ]            = _nSbeforeSim + size*count + i;
//                vertexID[ lVD ]               = _nSbeforeSim + size*count + i;
//                vertexName[ lVD ]             = vertexName[ vd ];
//                vertexTexName[ lVD ]          = vertexTexName[ vd ];
//                vertexScale[ lVD ]            = vertexScale[ vd ];
//                vertexSelectMag[ lVD ]        = false;
//                vertexIsStation[ lVD ]        = false;
//
//
//                _nStations++;
//                _nLabels++;
//
//                // add sub label edge
//                pair<EdgeDescriptor, unsigned int> lE = add_edge( preVD, lVD,
//                                                                  1.0, graph );
//                EdgeDescriptor lED = lE.first;
//                Coord2 coordA = vertexCoord[ preVD ];
//                Coord2 coordB = vertexCoord[ lVD ];
//                double diffX = coordB.x() - coordA.x();
//                double diffY = coordB.y() - coordA.y();
//                double angle = atan2( diffY, diffX );
//
//                edgeIndex[ lED ]             = _nEbeforeSim + (size+1)*count + i;
//                edgeID[ lED ]                = _nEbeforeSim + (size+1)*count + i;
//                edgeWeight[ lED ]            = vertexScale[ vd ];
//                edgeGeoAngle[ lED ]          = angle;
//                edgeSmoAngle[ lED ]          = angle;
//                edgeCurAngle[ lED ]          = angle;
//                edgeIsLine[ lED ]            = false;
//
//                preVD = lVD;
//                _nEdges++;
//            }
//
//            // add sub label edge
//            pair<EdgeDescriptor, unsigned int> lE = add_edge( preVD, vd,
//                                                              1.0, graph );
//            EdgeDescriptor lED = lE.first;
//            Coord2 coordA = vertexCoord[ preVD ];
//            Coord2 coordB = vertexCoord[ vd ];
//            double diffX = coordB.x() - coordA.x();
//            double diffY = coordB.y() - coordA.y();
//            double angle = atan2( diffY, diffX );
//
//            edgeIndex[ lED ]             = _nEbeforeSim + (size+1)*count + 4;
//            edgeID[ lED ]                = _nEbeforeSim + (size+1)*count + 4;
//            edgeWeight[ lED ]            = vertexScale[ vd ];
//            edgeGeoAngle[ lED ]          = angle;
//            edgeSmoAngle[ lED ]          = angle;
//            edgeCurAngle[ lED ]          = angle;
//            edgeIsLine[ lED ]            = false;
//
//            _nEdges++;
//            count++;
//        }
//    }
}

//
//  Metro::cloneSmooth --    copy all smooth vertex position to the input metro
//
//  Inputs
//  obj: original metro layout
//
//  Outputs
//  none
//
void Metro::cloneSmooth( const Metro & obj )
{
//    UndirectedGraph               sGraph          = obj.g();
//    VertexIndexMap      sVertexIndex    = get( vertex_index, sGraph );
//    VertexSmoothMap     sVertexSmooth   = get( vertex_mysmooth, sGraph );
//    VertexCoordMap      sVertexCoord    = get( vertex_mycoord, sGraph );
//    VertexIsStationMap  sVertexIsStation= get( vertex_myisstation, sGraph );
//
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexTempMap       vertexTemp      = get( vertex_mytemp, graph );
//    VertexExternalMap   vertexExternal  = get( vertex_myexternal, graph );
//    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, graph );
//    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//
//    // copy stations
//    BGL_FORALL_VERTICES( vdS, sGraph, UndirectedGraph ){
//        BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//            if ( vertexIndex[ vd ] == sVertexIndex[ vdS ] ) {
//                vertexSmooth[ vd ] = sVertexSmooth[ vdS ];
//                //vertexTemp[ vd ] = vertexCoord[ vdS ];
//                vertexCoord[ vd ] = sVertexCoord[ vdS ];
//            }
//        }
//    }
//
//    // copy labels
//    vector< VertexDescriptor > vdSVec;
//    BGL_FORALL_VERTICES( vdS, sGraph, UndirectedGraph ){
//
//        if( sVertexIsStation[ vdS ] == false && vdSVec.size() < obj.nLabels() ) {
//            vdSVec.push_back( vdS );
//        }
//    }
//    // cerr << "id = " << vdSVec.size() << " nLabels = " << obj.nLabels() << endl;
//    int count = 0;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//
//        if( vertexSelectMag[ vd ] == true && vertexExtstate[ vd ] == true && count < vdSVec.size() ) {
//
//            VertexDescriptor vdL = vdSVec[ count ];
//            Coord2 coord = sVertexCoord[ vdL ];
//#ifdef  SKIP
//            // compute min radius
//            double minR = INFINITY;
//            OutEdgeIterator e, e_end;
//            for ( tie( e, e_end ) = out_edges( vdL, sGraph ); e != e_end; ++e ) {
//                EdgeDescriptor ed = *e;
//                VertexDescriptor vS = source( ed, sGraph );
//                VertexDescriptor vT = target( ed, sGraph );
//                double radius = ( sVertexCoord[ vS ] - sVertexCoord[ vT ] ).norm();
//                if( vertexIsStation[ vS ] == false && vertexIsStation[ vT ] == false )
//                    radius /= 2.0;
//                if( radius < minR ) minR = radius;
//            }
//#endif  // SKIP
//
//            vertexExternal[ vd ].curSite() = coord;
//            //vertexExternal[ vd ].width() = minR * sqrt( 2.0 )/2.0;
//            //vertexExternal[ vd ].height() = minR * sqrt( 2.0 )/2.0;
//            vertexExternal[ vd ].width() = (vertexCoord[ vd ] - coord ).norm()/2.0;
//            vertexExternal[ vd ].height() = (vertexCoord[ vd ] - coord ).norm()/2.0;
//            vertexExternal[ vd ].leaderWeight() = (vertexCoord[ vd ] - coord ).norm()/_distBeta;
//            // cerr << "leaderW = " << vertexExternal[ vd ].leaderWeight() << endl;
//            // cerr << "minR = " << minR << endl;
//            // if( ( count + 4 ) < vdSVec.size() ) count = count + 4;
//            count = count + 1;
//        }
//    }
    //printGraph( graph );
}

//
//  Metro::cloneOctilinear --    copy all octilinear vertex position to the input metro
//
//  Inputs
//  obj: original metro layout
//
//  Outputs
//  none
//
void Metro::cloneOctilinear( const Metro & obj )
{
//    UndirectedGraph               sGraph          = obj.g();
//    VertexIndexMap      sVertexIndex    = get( vertex_index, sGraph );
//    VertexCoordMap      sVertexCoord    = get( vertex_mycoord, sGraph );
//    VertexIsStationMap  sVertexIsStation= get( vertex_myisstation, sGraph );
//
//    VertexIndexMap      vertexIndex     = get( vertex_index, graph );
//    VertexCoordMap      vertexCoord     = get( vertex_mycoord, graph );
//    VertexExternalMap   vertexExternal  = get( vertex_myexternal, graph );
//    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, graph );
//    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, graph );
//    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, graph );
//
//    // copy stations
//    BGL_FORALL_VERTICES( vdS, sGraph, UndirectedGraph ){
//        BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//            if ( vertexIndex[ vd ] == sVertexIndex[ vdS ] ) {
//                vertexCoord[ vd ] = sVertexCoord[ vdS ];
//            }
//        }
//    }
//
//    // copy labels
//    vector< VertexDescriptor > vdSVec;
//    BGL_FORALL_VERTICES( vdS, sGraph, UndirectedGraph ){
//
//        if( sVertexIsStation[ vdS ] == false && vdSVec.size() < obj.nLabels() ) {
//            vdSVec.push_back( vdS );
//        }
//    }
//    //cerr << "id = " << vdSVec.size() << " nLabels = " << obj.nLabels() << endl;
//    int count = 0;
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ){
//
//        if( vertexSelectMag[ vd ] == true && vertexExtstate[ vd ] == true && count < vdSVec.size() ) {
//
//            VertexDescriptor vdL = vdSVec[ count ];
//            Coord2 coord = sVertexCoord[ vdL ];
//#ifdef  SKIP
//            // compute min radius
//            double minR = INFINITY;
//            OutEdgeIterator e, e_end;
//            for ( tie( e, e_end ) = out_edges( vdL, sGraph ); e != e_end; ++e ) {
//                EdgeDescriptor ed = *e;
//                VertexDescriptor vS = source( ed, sGraph );
//                VertexDescriptor vT = target( ed, sGraph );
//                double radius = ( sVertexCoord[ vS ] - sVertexCoord[ vT ] ).norm();
//                if( vertexIsStation[ vS ] == false && vertexIsStation[ vT ] == false )
//                    radius /= 2.0;
//                if( radius < minR ) minR = radius;
//            }
//
//            // compute min distance to metro lines
//            BGL_FORALL_EDGES( ed, graph, UndirectedGraph ){
//
//                VertexDescriptor vdS = source( ed, graph );
//                VertexDescriptor vdT = target( ed, graph );
//                Coord2 & pointV = coord;
//                Coord2 & pointS = vertexCoord[ vdS ];
//                Coord2 & pointT = vertexCoord[ vdT ];
//
//                double m = ( pointS.y() - pointT.y() ) / ( pointS.x() - pointT.x() );
//                double k = pointS.y() - m * pointS.x();
//                double dist = fabs( m * pointV.x() - pointV.y() + k ) / sqrt( SQUARE( m ) + 1.0 );
//                Coord2 vec, p;
//                vec.x() = ( pointS - pointT ).y();
//                vec.y() = -( pointS - pointT ).x();
//                vec = vec / vec.norm();
//                p = pointV + dist * vec;
//                if( fabs( p.y() - m * p.x() - k ) > 1.0e-8 ){
//                    p = pointV - dist * vec;
//                }
//                double r = ( p - pointT ).x() / ( pointS - pointT ).x() ;
//                if( r > 0 && ( 1-r ) > 0 ) {
//                    //cerr << "r = " << r << ", " << vertexIndex[ vdS ] << ", " << vertexIndex[ vdT ] << endl;
//                    if( minR > dist ) minR = dist;
//                }
//            }
//#endif  // SKIP
//            vertexExternal[ vd ].curSite() = coord;
//            //vertexExternal[ vd ].width() = minR * sqrt( 2.0 )/2.0;
//            //vertexExternal[ vd ].height() = minR * sqrt( 2.0 )/2.0;
//            vertexExternal[ vd ].width() = (vertexCoord[ vd ] - coord ).norm()/2.0;
//            vertexExternal[ vd ].height() = (vertexCoord[ vd ] - coord ).norm()/2.0;
//            // cerr << "minR = " << minR << endl;
//            // if( ( count + 4 ) < vdSVec.size() ) count = count + 4;
//            count = count + 1;
//        }
//    }
}

void Metro::updateTempCoord( void )
{
//    VertexIDMap             vertexID            = get( vertex_myid, graph );
//    VertexCoordMap          vertexCoord         = get( vertex_mycoord, graph );
//    VertexTempMap           vertexTemp          = get( vertex_mytemp, graph );
//
//    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ) {
//        vertexTemp[ vd ] = vertexCoord[ vd ];
//    }
}

//
//  Metro::load --    load the metro information
//
//  Inputs
//  filename: file name of the network
//
//  Outputs
//  none
//
void Metro::load( const string filename )
{
    clear();

    ifstream    ifs( filename.c_str() );
    char        buf[ MAX_STR ];

    if ( !ifs ) {
        cerr << "Cannot open the target file : " << filename << endl;
        return;
    }

    _line.clear();
    _lineSta.clear();
    _nLines = 0;
    _nStations = 0;
    _nEdges = 0;

    while ( true ) {

        istringstream istr;
        char prompt[ MAX_STR ], argument[ MAX_STR ];

        //------------------------------------------------------------------------------
        //  Read line name
        //------------------------------------------------------------------------------

        ifs.getline( buf, sizeof( buf ) );
        if ( strcmp( "End", buf ) == 0 ) break;

        istr.clear();
        istr.str( buf );
        istr >> prompt >> argument;

        assert( strcmp( prompt, "Line:" ) == 0 );
        strcpy( _lineName[ _nLines ], argument );

        //------------------------------------------------------------------------------
        //  Read line color
        //------------------------------------------------------------------------------
        double r, g, b;

        ifs.getline( buf, sizeof( buf ) );
        istr.clear();
        istr.str( buf );
        istr >> prompt >> r >> g >> b;
        assert( strcmp( prompt, "Color:" ) == 0 );
        _lineColor[ _nLines ][ 0 ] = r/255.0;
        _lineColor[ _nLines ][ 1 ] = g/255.0;
        _lineColor[ _nLines ][ 2 ] = b/255.0;

        //------------------------------------------------------------------------------
        //  Read station data
        //------------------------------------------------------------------------------
        vector< VertexDescriptor > ptrSta;
        while ( true ) {

            double x, y, pri, w, h;

            ifs.getline( buf, sizeof( buf ) );
            // the end of station info.
            if ( buf[ 0 ] == '#' ) {
                break;
            }

            istr.clear();
            istr.str( buf );
            istr >> argument >> y >> x >> pri >> w >> h;
            //cout << x <<", " << y << "\n";


            // Check if the station exists or not
            VertexDescriptor curVD = NULL;
            BGL_FORALL_VERTICES( vd, graph, UndirectedGraph )
            {
                string name = * graph[ vd ].namePtr;
                if( strcmp( argument, name.c_str() ) == 0 ){


                    curVD = vd;
                    break;
                }
            }

            if ( curVD == NULL ){

                curVD = add_vertex( graph );
                vector< unsigned int > lineID = graph[ curVD ].lineID;
                lineID.push_back( _nLines );

	            graph[ curVD ].id               = _nStations;

	            graph[ curVD ].coordPtr         = new Coord2( x, y );
	            graph[ curVD ].geoPtr           = new Coord2( x, y );
	            graph[ curVD ].smoothPtr        = new Coord2( x, y );
	            graph[ curVD ].octilinearPtr    = new Coord2( x, y );
	            graph[ curVD ].namePtr          = new string( argument );
                graph[ curVD ].metroShape       = true;
                graph[ curVD ].inflectionPoint  = false;    // todo
                graph[ curVD ].smoothPath       = false;
                graph[ curVD ].interstate       = false;

				graph[ curVD ].weight           = 1.0;

				graph[ curVD ].lineID.push_back( _nLines );
                _nStations++;
            }
            else{
                graph[ curVD ].lineID.push_back( _nLines );
            }

#ifdef DEBUG

#endif  // DEBUG
            ptrSta.push_back( curVD );
        }
        _lineSta.push_back( ptrSta );

        //------------------------------------------------------------------------------
        //  Read line data
        //------------------------------------------------------------------------------
        vector< EdgeDescriptor > eachline;

        while ( true ) {
            int     origID, destID;
            double  weight;

            ifs.getline( buf, sizeof( buf ) );

            // the end of station info.
            if ( buf[ 0 ] == '#' ) {
                break;
            }

            istr.clear();
            istr.str( buf );
            istr >> origID >> destID >> weight;


            // search previous edge
            bool found = false;
            EdgeDescriptor oldED;
            VertexDescriptor oldVS = ptrSta[ origID ];
            VertexDescriptor oldVT = ptrSta[ destID ];
            unsigned int indexS = graph[ oldVS ].id;
            unsigned int indexT = graph[ oldVT ].id;

            tie( oldED, found ) = edge( oldVS, oldVT, graph );

            if ( found == true ) {

                graph[ oldED ].lineID.push_back( _nLines );
                eachline.push_back( oldED );
                bool test = false;
                tie( oldED, test ) = edge( oldVT, oldVS, graph );
            }
            else{
                VertexDescriptor src = vertex( indexS, graph );
                VertexDescriptor tar = vertex( indexT, graph );

                // handle fore edge
                pair<EdgeDescriptor, unsigned int> foreE = add_edge( src, tar, graph );
                EdgeDescriptor foreED = foreE.first;

                // calculate geographical angle
                Coord2 coordO;
                Coord2 coordD;
                if( graph[ src ].id < graph[ tar ].id ){
                    coordO = *graph[ src ].coordPtr;
                    coordD = *graph[ tar ].coordPtr;
                }
                else{
                    coordO = *graph[ tar ].coordPtr;
                    coordD = *graph[ src ].coordPtr;
                }
                double diffX = coordD.x() - coordO.x();
                double diffY = coordD.y() - coordO.y();
                double angle = atan2( diffY, diffX );

                graph[ foreED ].id          = _nEdges;
                graph[ foreED ].weight      = weight;
                graph[ foreED ].geoAngle    = angle;
                graph[ foreED ].smoAngle    = angle;
                graph[ foreED ].curAngle    = angle;
                graph[ foreED ].lineID.push_back( _nLines );
                graph[ foreED ].isCurved    = false;
                graph[ foreED ].isMetro     = true;
                graph[ foreED ].matchPath   = false;
                graph[ foreED ].isSmoothPath = false;

                eachline.push_back( foreED );
                _nEdges++;
            }

        }

        _line.push_back( eachline );
        _shapeLines.push_back(true); // init all with part of metro Shape
        _nLines++;

    }
    ifs.close();

    // Label instates
    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ) {
        DegreeSizeType degree = out_degree( vd, graph );
        if (degree > 2 )
            graph[vd].interstate = true;
        if (degree == 0)
            cerr << "degree 0 exists" << *graph[vd].namePtr << endl;
    }
    vector<VertexDescriptor> removeVertex;
    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph ) {
        DegreeSizeType degree = out_degree( vd, graph );
        if (degree < 1 ) {
            removeVertex.push_back(vd);
        }
    }
    for (VertexDescriptor vd : removeVertex) {
        // remove_vertex(vd, graph);
        // _nStations--;
        // cerr << "removing a circle" << endl;
    }

    return;
}

//
//  Metro::loadLabel --    load label information of the network
//
//  Inputs
//  filename: file name of the network label
//
//  Outputs
//  none
//
void Metro::loadLabel( const string filename ){}

//
//  Metro::exportData --  export metro data
//
//  Inputs
//      none
//
//  Outputs
//      none
//
void Metro::exportData( const string filename )
{
	ofstream            ofs( filename.c_str() );

	if ( !ofs ) {
		cerr << "Cannot open the target file : " << filename << endl;
		return;
	}
	cerr << " now writing the metro network and stations... " << endl;

	vector< unsigned int > conv( _nStations );

	for ( unsigned int nL = 0; nL < _nLines; ++nL ) {

		//------------------------------------------------------------------------------
		//      Header of each line
		//------------------------------------------------------------------------------
		ofs << "Line:" << "\t" << _lineName[ nL ] << endl;
		ofs << "Color:" << "\t"
		    << _lineColor[ nL ][ 0 ] * 255.0  << "\t"
		    << _lineColor[ nL ][ 1 ] * 255.0  << "\t"
		    << _lineColor[ nL ][ 2 ] * 255.0  << endl;

		//------------------------------------------------------------------------------
		//      Prepare convertion table for station IDs
		//------------------------------------------------------------------------------
		for ( unsigned int m = 0; m < _nStations; ++m )
			conv[ m ] = _nStations;

		VertexDescriptor vd = NULL;
		// for each station on the line
		for ( unsigned int k = 0; k < _lineSta[ nL ].size(); ++k ) {
			vd = _lineSta[ nL ][ k ];
			ofs << * graph[ vd ].namePtr << "\t"
			    << graph[ vd ].coordPtr->y() << "\t"
			    << graph[ vd ].coordPtr->x() << "\t"
			    << "1.0" << "\t"
			    << "1.0" << "\t"
			    << "1.0" << endl;
			// set the conversion table
			conv[ graph[ vd ].id ] = k;
		}
		ofs << "#" << endl;

		for ( unsigned int k = 0; k < _line[ nL ].size(); ++k ) {
			VertexDescriptor vdSrc = source( _line[ nL ][ k ], graph );
			VertexDescriptor vdTar = target( _line[ nL ][ k ], graph );
			ofs << conv[ graph[ vdSrc ].id ] << "\t" << conv[ graph[ vdTar ].id ] << "\t"
			    << graph[ _line[ nL ][ k ] ].weight << endl;
		}
		ofs << "#" << endl;
	}
	ofs << "End" << endl;

	ofs.close();
}

bool Metro::isOctolinearTurningPoint(VertexDescriptor vd) {
    bool turningpoint = true;
    DegreeSizeType degrees = out_degree( vd, graph );
    if (degrees == 2 and !graph[vd].smoothPath) {
        OutEdgeIterator e, e_end;
        double tarAngle = 3 * M_PI;
        for ( tie( e, e_end ) = out_edges( vd, graph ); e != e_end; ++e ) {
            EdgeDescriptor ed = *e;
            double tar = graph[ed].target;
            if (tarAngle > 2 * M_PI) {
                tarAngle = tar;
            } else if (tarAngle == tar or
            tarAngle == tar + M_PI or
            tarAngle == tar - M_PI) {
                turningpoint = false;
            }
        }
    }
    return turningpoint;
}

//
//  Metro::exportGraphML --  export metro data in graphml format
//
//  Inputs
//      none
//
//  Outputs
//      none
//
void Metro::exportGraphML( const string filename )
{
    ofstream            ofs( filename.c_str() );

    double ww = ZuKai::Base::Common::getMainwidgetWidth();
    double wh = ZuKai::Base::Common::getMainwidgetHeight();

    if ( !ofs ) {
        cerr << "Cannot open the target file : " << filename << endl;
        return;
    }
    cerr << " now writing the metro network and stations... " << endl;

    ofs << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
           "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
           "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
           "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n"
           "     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
           "  <key id=\"x\" for=\"node\" attr.name=\"x_coordinate\" attr.type=\"int\"/>\n"
           "  <key id=\"y\" for=\"node\" attr.name=\"y_coordinate\" attr.type=\"int\"/>\n"
           "  <key id=\"label\" for=\"node\" attr.name=\"station_name\" attr.type=\"string\"/>\n"
           "  <key id=\"line\" for=\"edge\" attr.name=\"line\" attr.type=\"string\"/>\n"
           "  <key id=\"smooth\" for=\"edge\" attr.name=\"smooth\" attr.type=\"boolean\"/>\n"
           "  <key id=\"turningpoint\" for=\"node\" attr.name=\"turningpoint\" attr.type=\"boolean\"/>\n"
           "  <key id=\"smoothNode\" for=\"node\" attr.name=\"smoothNode\" attr.type=\"boolean\"/>"
        << endl;

    // export line information
    for ( unsigned int nL = 0; nL < _nLines; ++nL ) {

        //------------------------------------------------------------------------------
        //      Header of each line
        //------------------------------------------------------------------------------
        ofs << "  <key id=\"l" << nL
            << "\" for=\"edge\" attr.name=\"" << _lineName[ nL ]
            << "\" attr.type=\"boolean\" color.r=\"" << _lineColor[ nL ][ 0 ] * 255
            << "\" color.g=\"" << _lineColor[ nL ][ 1 ] * 255
            << "\" color.b=\"" << _lineColor[ nL ][ 2 ] * 255
            << "\">\n"
            << "    <default>FALSE</default>\n"
            << "  </key>"
            << endl;
    }

    // export graph information
    ofs << "  <graph id=\"G\" edgedefault=\"undirected\">" << endl;

    // export vertex information
    BGL_FORALL_VERTICES( vd, graph, UndirectedGraph )
        {
            string turningPoint = "true";
            if (!isOctolinearTurningPoint(vd))
                turningPoint = "false";

            string smoothNode = "false";
            if (graph[vd].smoothPath)
                smoothNode = "true";
            
            string id = "n" + std::to_string( graph[vd].id );
            Coord2 & coord = *graph[ vd ].coordPtr;
            string & name = *graph[vd].namePtr;
            ofs << "    <node id=\"" << id << "\">" << endl;
            ofs << "      <data key=\"x\">" << (int)(coord.x() + ww) << "</data>" << endl;
            ofs << "      <data key=\"y\">" << (int)(coord.y() + wh) << "</data>" << endl;
            ofs << "      <data key=\"label\">" << name << "</data>" << endl;
            ofs << "      <data key=\"turningpoint\">" << turningPoint << "</data>" << endl;
            ofs << "      <data key=\"smoothNode\">" << smoothNode << "</data>" << endl;
            ofs << "    </node>" << endl;
        }
    ofs << endl;

    // export edge information

    // with a For all edges loop
    BGL_FORALL_EDGES( edge, graph, UndirectedGraph ){

        string lineIdString = "";
        auto lineIDs = graph[edge].lineID;
        for (auto lineID : lineIDs) {
            if (lineIdString.compare("") != 0)
                lineIdString += ",";
            lineIdString += "l" + to_string(lineID);
        }

        VertexDescriptor sVD = source(edge, graph);
        VertexDescriptor tVD = target(edge, graph);

        int id = graph[edge].id;
        int sid = graph[sVD].id;
        int tid = graph[tVD].id;

        string path = "false";
        if (graph[edge].isSmoothPath)
            path = "true";


        ofs << "    <edge id=\"e"  << id
        << "\" source=\"n" << sid
        << "\" target=\"n" << tid
        << "\">" << endl;
        ofs << "      <data key=\"line\">" << lineIdString << "</data>" << endl;
        ofs << "      <data key=\"smooth\">" << path << "</data>" << endl;
        ofs << "    </edge>" << endl;
    }

    // using the _nLines
    /*
    unsigned int counter = 0;
    for ( unsigned int nL = 0; nL < _nLines; ++nL ) {

        for ( unsigned int k = 1; k < _lineSta[ nL ].size(); ++k ) {

            string id = "e" + std::to_string( counter );
            string sid = "n" + std::to_string( graph[_lineSta[ nL ][ k-1 ]].id );
            string tid = "n" + std::to_string( graph[_lineSta[ nL ][ k ]].id );
            auto sidI = graph[_lineSta[ nL ][ k-1 ]].id;
            auto tidI = graph[_lineSta[ nL ][ k ]].id;
            VertexDescriptor sVD = NULL;
            VertexDescriptor tVD = NULL;
            BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) { // Sorry
                if (graph[vertex].id == sidI)
                    sVD = vertex;
                else if (graph[vertex].id == tidI)
                    tVD = vertex;
            }

            // checking if edge is smooth path
            string path = "false";
            if (sVD != NULL and tVD != NULL ) {
                auto edgePair = boost::edge(sVD, tVD, graph);

                if (edgePair.second) {
                    if (graph[edgePair.first].isSmoothPath) {
                        path = "true";
                    }
                }
            }


            // string path = std::to_string(graph[_lineSta[ nL ][ k ]].smoothPath);
            ofs << "    <edge id=\"" << id
                << "\" source=\"" << sid
                << "\" target=\"" << tid
                << "\">" << endl;
            // ofs << "      <data key=\"l" << nL
            //     << "\">true</data>" << endl;
            // ofs << "      <data key=\"smooth\">" << path << "</data>" << endl;

            ofs << "      <data key=\"line\">l" << nL << "</data>" << endl;
            ofs << "      <data key=\"smooth\">" << path << "</data>" << endl;
            ofs << "    </edge>" << endl;

            counter++;
        }
    }
    */
    ofs << "  </graph>" << endl;

    ofs << "</graphml>"
        << endl;

    cerr << "exporting " << filename << "..." << endl;

    ofs.close();
}
// end of header file
// Do not add any stuff under this line.

