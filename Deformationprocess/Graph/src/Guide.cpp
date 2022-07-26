/*
*	Guid.cpp
*   Based on FocusContext/Smooth.cpp
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   Guide Path / Guide Graph 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cmath>

using namespace std; 

#include "Guide.h"

//------------------------------------------------------------------------------
//	Private functions
//------------------------------------------------------------------------------

/*
* Adds new Vertex to Graph
* 
* Inputs
*   double x position.y
*   double y position.x
*   pointer for new Vertex
*   int id of new vertex
*
* Outputs
*   none
*/
VertexDescriptor Guide::setNewVertex(double x, double y, vector < VertexDescriptor > &ptrPoints, int _id) {
    VertexDescriptor curVD = add_vertex( graph );
    graph[ curVD ].geoPtr = new Coord2( x, y );
    graph[ curVD].coordPtr = new Coord2( x, y ); 
    graph[ curVD ].weight = 1.0;
    graph[ curVD ].id = _id;
    graph[ curVD ].inflectionPoint = false;
    ptrPoints.push_back(curVD);
    return curVD;
}

/*
* Calculates magnitude of a vector from [0,0] to v
*
* Inputs    Coordinates of point
* Outputs   Magnitude
*/
double Guide::magnitude( Coord2 v ) {
    double d = v.x() * v.x() + v.y() * v.y(); 
    return sqrt(d);
}


//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------
Guide::Guide( void ) {}

Guide::Guide( const Guide & obj ) {}

Guide::~Guide( void ) {}

/*
* Load new Graph from filen
* 
* Input 
*   string filename path to file 
* 
* Outputs
*   none
*/
void Guide::load( string filename) {
    ifstream(); 
    ifstream    ifs(filename.c_str());
    char        buf[ MAX_STR ];
    _nVertex = 0; 
    _nEdges = 0;
    _nLines = 0;

    if( !ifs ) {
        cerr << "Cannnot open the target file: " << filename << endl; 
        return; 
    } 

    istringstream istr; 
    vector <VertexDescriptor> ptrPoints;
    vector <VertexDescriptor> allPoints;
    string name, prompt, value;
    int _nGuids = 0; 

    while( true ) {
        ifs.getline(buf, sizeof(buf)); 
        istr.clear(); 
        istr.str(buf); 
        istr >> prompt >> value;
        cerr << prompt << " " << value << endl;     //TODO
        
        if (prompt.compare("End") == 0) break;
        else if (prompt.compare("Line:") == 0) {
            name = value;
            _lineName[_nLines] = value;
        }
        vector< VertexDescriptor > ptrSta;
        while ( true ) {    //Read in Vertex
            double x,y, pri;
            _nGuids++; 
            ifs.getline(buf, sizeof(buf));
            if (buf[ 0 ] == '#') break; 
            istr.clear(); 
            istr.str(buf); 
            istr >> y >> x >> pri;

            // check if allready one close enough
            VertexDescriptor toCloseVD = NULL;
            for (auto prevVD : allPoints ) {
                Coord2 prevCord = *graph[prevVD].coordPtr;
                Coord2 delta = Coord2(prevCord.x() - x , prevCord.y() - y);
                if (delta.norm() < MIN_NODE_DIST) {
                    toCloseVD = prevVD;
                    break;
                }
            }

            VertexDescriptor cVD;
            if (toCloseVD != NULL ) {
                ptrSta.push_back(toCloseVD);
                ptrPoints.push_back(toCloseVD);
                allPoints.push_back(toCloseVD);
                // _nVertex++;
            } else {
                cVD = setNewVertex(x, y, ptrPoints, _nVertex);
                ptrSta.push_back(cVD);
                allPoints.push_back(cVD);
                _nVertex++;
            }

        }
        _lineSta.push_back(ptrSta);

        while (true) {      //Read in Edges
            int sor, tar, w; 
            ifs.getline(buf, sizeof(buf)); 
            if (buf[0] == '#') break; 
            istr.clear(); 
            istr.str(buf); 
            istr >> sor >> tar >> w; 

            add_edge(ptrPoints[sor], ptrPoints[tar], graph); 
            _nEdges++; 
        }
        ptrPoints.clear();
        _nLines++;

    }

    // Check if complex Shape
    complexGuide = false;
    if (_nLines > 1 ) {
        complexGuide = true;
        cout << "has more than 1 line" << endl;
    }
    else {
        BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
        {
            auto degrees = out_degree( vertex, graph );
            if (degrees > 2 ) {
                complexGuide = true;
                cout << "has a vertex with degree > 2" << endl;
                break;
            }
        }
    }
    cout << "is a complex guide: " << complexGuide << endl;
    if (complexGuide) {
        // find most left vertex
        VertexDescriptor minXDS = NULL;
        double minX = 10e10;
        BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
        {
            Coord2 coord = *graph[vertex].coordPtr;
            if (minXDS == NULL || minX > coord.x()) {
                minXDS = vertex;
                minX = coord.x();
            }
        }
        outerPath.push_back(minXDS);

        Coord2 centerCoord = getCenterCoord();
        VertexDescriptor curVD = minXDS;
        double prevAngle = 0;
        while(true) {
            auto neighbours = adjacent_vertices(curVD, graph);
            vector<VertexDescriptor> possibleneighbours;
            for (auto vd : make_iterator_range(neighbours)) {

                // check if in outerpath
                if(std::find(outerPath.begin(), outerPath.end(), vd) != outerPath.end()) {
                    /* outerPath contains vd */
                } else {
                    // if not add to possible list
                    // std::cout << "curVD has adjacent vertex " << graph[vd].id << "\n";
                    possibleneighbours.push_back(vd);
                }
            }
            if (possibleneighbours.size() < 1 )
                break;
            else if (possibleneighbours.size() == 1 ) {
                auto vd = possibleneighbours[0];
                Coord2 nextEdge = (*graph[curVD].coordPtr - *graph[vd].coordPtr).normalize();
                double angle = atan2(nextEdge.y(), nextEdge.y());
                if (angle < 0 ) angle += 2 * M_PI;

                prevAngle = angle;
                curVD = vd;
                outerPath.push_back(curVD);
            }
            else { // turn as much left as possbile
                double largestAngle = 0;
                double largestDeltaAngle = 0;
                VertexDescriptor nextVd = NULL;
                for (auto vd : possibleneighbours) {
                    Coord2 delta = (centerCoord - *graph[vd].coordPtr).normalize();
                    Coord2 nextEdge = (*graph[curVD].coordPtr - *graph[vd].coordPtr).normalize();
                    // double angle = atan2(delta.y(), delta.x()) - atan2(nextEdge.y(), nextEdge.y());
                    double angle = atan2(nextEdge.y(), nextEdge.y());
                    if (angle < 0 ) angle += 2 * M_PI;
                    double deltaAngle = abs(prevAngle - angle);

                    if (deltaAngle > largestDeltaAngle || nextVd == NULL ) {
                        largestAngle = angle;
                        largestDeltaAngle = deltaAngle;
                        nextVd = vd;
                    }
                }

                prevAngle = largestAngle;
                curVD = nextVd;
                outerPath.push_back(curVD);
            }
        }
    }
    else {
        BGL_FORALL_VERTICES( vertex, graph, UndirectedGraph )
        {
            outerPath.push_back(vertex);
        }
    }
}

/**
* closest point to a vertex laying on a edge
* 
* Inputs
*    Coord2 source of the edge
*    Coord2 target of the edge
*    Coord2 vertex position 
* Outputs
*   Coord2 closest point on the edge
*/
Coord2 Guide::closestPointOnEdge( Coord2 source, Coord2 target, Coord2 vertex ) {
    Coord2 e = target - source;
    double eNorm = e.norm();
    if (fabs(eNorm) < 0.0001) return source;
    e = e / e.norm();
    Coord2 sv = vertex - source; 
    Coord2 closestPoint = (sv.x() * e.x() + sv.y() * e.y()) * e;    //sv.dot(e) * e;
    closestPoint = source + closestPoint;
    if (closestPoint[0] < min(source[0], target[0])) {
        if (source[0] < target[0]) closestPoint = source; 
        else closestPoint = target;
    } else if (closestPoint[0] > max(source[0], target[0])) {
        if (source[0] > target[0] ) closestPoint = source; 
        else closestPoint = target; 
    }
    if (isnan(closestPoint.x())) cerr << "Closest Point on X is nan" << endl;
    return closestPoint; 
}

/*
* closest point to a vertex laying on a edge
* 
* Inputs
*    Edgediscriptor of edge
*    Coord2 vertex position 
* Outputs
*   Coord2 closest point on the edge
*/
Coord2 Guide::closestPointOnEdge( EdgeDiscriptor edge, Coord2 vertex ) {
    VertexDescriptor svd = source(edge, graph); 
    VertexDescriptor tvd = target(edge, graph);
    Coord2 source = *graph[ svd ].coordPtr; //xy1
    Coord2 target = *graph[ tvd ].coordPtr;  //xy2

    Coord2 e = target - source;
    double eNorm = e.norm();
    if (fabs(eNorm) < 0.01 ) {
        return source;
    }
    e = e / eNorm; //e.normalize(); // e / e.norm();
    Coord2 sv = vertex - source;
    Coord2 closestPoint = (sv.x() * e.x() + sv.y() * e.y()) * e;    //sv.dot(e) * e;
    closestPoint = source + closestPoint;

    // if ( abs(source.x() - target.x() ) < 1. ) closestPoint.setX( vertex.x() );
    // if ( abs(source.y() - target.y() ) < 1. ) closestPoint.setY( vertex.y() );

    // checking x
    if (closestPoint.x() < min(source.x(), target.x())) {
        if (source.x() < target.x()) closestPoint = source;
        else closestPoint = target;
    } else if (closestPoint.x() > max(source.x(), target.x())) {
        if (source.x() > target.x() ) closestPoint = source;
        else closestPoint = target; 
    }
    // checking y
    else if ( closestPoint.y() < min(source.y(), target.y()) ) {
        if ( source.y() < target.y() ) closestPoint = source;
        else closestPoint = target;
    } else if ( closestPoint.y() > max(source.y(), target.y()) ) {
        if (source.y() > target.y() ) closestPoint = source;
        else closestPoint = target;
    }
    if (isnan(closestPoint.x())) cerr << "Closest Point on X is nan" << endl;
    return closestPoint; 
}

/*
* calculates the minimum distance to an edge
* 
* Inputs 
*   edgediscriptor of edge
*   Coord2 of vertex
* 
* Outputs
*   distance to edge
*/
double Guide::distanceToClosestPointOnEdge( EdgeDiscriptor edge, Coord2 vertex ) {
    Coord2 v = vertex - closestPointOnEdge( edge, vertex); 
    return v.norm();
}

/* 
* calculates the geographical center of a guide
* 
* Inputs none
* Outputs 
*   Coordinates of the center
*/
Coord2 Guide::getCenterCoord() {
    double maxX = -10e10;
    double maxY = -10e+10;
    double minX = 10e10;
    double minY = 10e+10;
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 cor = *graph[vertex].coordPtr;
        maxX = max(maxX, cor.x());
        minX = min(minX, cor.x());
        maxY = max(maxY, cor.y());
        minY = min(minY, cor.y());
    }
    Coord2 delta = Coord2( maxX - minX, maxY - minY );
    Coord2 center = Coord2( maxX - (delta.x() * 0.5), maxY - (delta.y() * .5));
    return center;
}

/**
 * Finds the shortes paths between two vertex (Smallest number of edges).
 * @param ui starting vertex
 * @param un end vertex
 * @return sequence of VertexDescriptor that belong to the shortest path.
 * Including ui, and un.
 */
vector<VertexDescriptor> Guide::shortesPath(VertexDescriptor ui, VertexDescriptor un) {
    vector<VertexDescriptor> result;
    if (ui == un) {
        result.push_back(un);
        result.push_back(ui);
        return result;
    }
    //create vertex set Q
    vector<VertexDescriptor> Q;
    map<VertexDescriptor, double> dist;
    map<VertexDescriptor, VertexDescriptor> prev;

    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        dist.insert(pair<VertexDescriptor, double>(vertex, 10e10));
        prev.insert(pair<VertexDescriptor, VertexDescriptor>(vertex, NULL));
        Q.push_back(vertex);
    }
    dist[ui] = 0.0;
    while(Q.size() > 0) {
        VertexDescriptor u = NULL;
        int eraseIndex = 0;
        for (int j = 0; j < Q.size(); j++) {
            auto vd = Q[j];
            if (dist[vd] < dist[u] || u == NULL ) {
                u = vd;
                eraseIndex = j;
            }
        }
        Q.erase(Q.begin() + eraseIndex);
        if (u == un)
            break;

        auto allNeighbours = adjacent_vertices(u, graph);
        vector<VertexDescriptor> neighbors;
        for (auto vd : make_iterator_range(allNeighbours)) {
            if(std::find(Q.begin(), Q.end(), vd) != Q.end()) {
                neighbors.push_back(vd);
            }
        }
        for (auto vd : neighbors) {
            double alt = dist[u] + 1; //dist[u] + length(u, v)
            if (alt < dist[vd]) {
                dist[vd] = alt;
                prev[vd] = u;
            }
        }
    }

    //Backtrack
    VertexDescriptor prevVD = un;
    while( prevVD != ui ) {
        result.push_back(prevVD);
        prevVD = prev[prevVD];
    }
    result.push_back(ui);
    return result;
}



/*
* Calculates the geo. closest edge to a vertex
* 
* Inputs
*   Coord2 of the vertex
* 
* Outputs
*   Edgdiscriptor of the closest edge to the vertex
*/
EdgeDiscriptor Guide::closestEdge(Coord2 vertex) {
    EdgeDiscriptor closestEdge;
    double minEdgeDistance = 10e+10;
    bool first = false;

    //find for each vertex the closest Edge 
    BGL_FORALL_EDGES(edge, graph, UndirectedGraph) {
        VertexDescriptor uS = source(edge, graph);
        VertexDescriptor uT = target(edge, graph);
        Coord2 closePoint = closestPointOnEdge(edge, vertex);
        // double distance = (closePoint - vertex).norm();
        Coord2 e = vertex - closePoint;
        double distance = e.norm();
            if ( (graph[uS].id == 7 && graph[uT].id == 8) ||
                 (graph[uT].id == 7 && graph[uS].id == 8) ) {
                    // cerr << "dostance: " << distance << endl;
            }
            
        if (distance < minEdgeDistance || first ) {
            minEdgeDistance = distance; 
            closestEdge = edge;
            first = false;
        }
    }
    return closestEdge; 
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
void Guide::exportGraphML( const string filename )
{

    ofstream            ofs( filename.c_str() );

    double ww = ZuKai::Base::Common::getMainwidgetWidth();
    double wh = ZuKai::Base::Common::getMainwidgetHeight();

    if ( !ofs ) {
        cerr << "Cannot open the target file : " << filename << endl;
        return;
    }
    cerr << " now writing the Guide network and stations... " << endl;

    ofs << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
           "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
           "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
           "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n"
           "     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
           "  <key id=\"x\" for=\"node\" attr.name=\"x_coordinate\" attr.type=\"int\"/>\n"
           "  <key id=\"y\" for=\"node\" attr.name=\"y_coordinate\" attr.type=\"int\"/>\n"
           "  <key id=\"label\" for=\"node\" attr.name=\"station_name\" attr.type=\"string\"/>\n"
           "  <key id=\"line\" for=\"edge\" attr.name=\"line\" attr.type=\"string\"/>"
        << endl;

    // export line information
    for ( unsigned int nL = 0; nL < _nLines; ++nL ) {

        //------------------------------------------------------------------------------
        //      Header of each line
        //------------------------------------------------------------------------------
        ofs << "  <key id=\"l" << nL
            << "\" for=\"edge\" attr.name=\"" << _lineName[ nL ]
            << "\" attr.type=\"boolean\" color.r=\"" << 30
            << "\" color.g=\"" << 30
            << "\" color.b=\"" << 30
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
            string id = "n" + std::to_string( graph[vd].id );
            Coord2 & coord = *graph[ vd ].coordPtr;
            string n = "guidevertex";
            string & name = n; // *graph[vd].namePtr;     //todo
            ofs << "    <node id=\"" << id << "\">" << endl;
            ofs << "      <data key=\"x\">" << (int)(coord.x() + ww) << "</data>" << endl;
            ofs << "      <data key=\"y\">" << (int)(coord.y() + wh) << "</data>" << endl;
            ofs << "      <data key=\"label\">" << name << "</data>" << endl;
            ofs << "    </node>" << endl;
        }
    ofs << endl;
    // export edge information
    unsigned int counter = 0;
    for ( unsigned int nL = 0; nL < _nLines; ++nL ) {

        for ( unsigned int k = 1; k < _lineSta[ nL ].size(); ++k ) {
            string id = "e" + std::to_string( counter );
            string sid = "n" + std::to_string( graph[_lineSta[ nL ][ k-1 ]].id );
            string tid = "n" + std::to_string( graph[_lineSta[ nL ][ k ]].id );
            ofs << "    <edge id=\"" << id
                << "\" source=\"" << sid
                << "\" target=\"" << tid
                << "\">" << endl;
            ofs << "      <data key=\"line\">" << "l" << nL
                << "</data>" << endl;
            ofs << "    </edge>" << endl;

            counter++;
        }
    }

    ofs << "  </graph>" << endl;

    ofs << "</graphml>"
        << endl;

    cerr << "exporting " << filename << "..." << endl;

    ofs.close();

}

double Guide::getHeight() {
    double minY = 10e10;
    double maxY = -10e10;
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 coord = *graph[vertex].coordPtr;
        minY = min(minY, coord.y());
        maxY = max(maxY, coord.y());
    }
    return maxY - minY;
}

double Guide::getWidth() {
    double minX = 10e10;
    double maxX = -10e10;
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 coord = *graph[vertex].coordPtr;
        minX = min(minX, coord.x());
        maxX = max(maxX, coord.x());
    }

    return maxX - minX;
}

void Guide::translate(Coord2 trans) {
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 v = *graph[vertex].coordPtr;
        graph[vertex].coordPtr->setX(v.x() + trans.x());
        graph[vertex].coordPtr->setY(v.y() + trans.y());
        graph[vertex].geoPtr->setX(v.x() + trans.x());
        graph[vertex].geoPtr->setY(v.y() + trans.y());
    }
}

void Guide::scale(double s) {
    BGL_FORALL_VERTICES(vertex, graph, UndirectedGraph) {
        Coord2 v = *graph[vertex].coordPtr;
        v *= s;
        graph[vertex].coordPtr->setX(v.x());
        graph[vertex].coordPtr->setY(v.y());
    }

}

// end of header file
// Do not add any stuff under this line.