/*
*	Smooth.cpp
*   Based on FocusContext/Smooth.cpp
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   Calculating Smooth integrating guide Layout 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cmath>

using namespace std;

#include "Smooth.h"

void Smooth::_init( Metro * __metro, Guide * __guide,  double __half_width, double __half_height )
{
    _metro                      = __metro;
    _guide                      = __guide;
    UndirectedGraph        & g            = _metro->g();

    // update smooth angle
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ){

        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );

        Coord2 coordO;
        Coord2 coordD;
        if( g[ vdS ].id < g[ vdT ].id ){
            coordO = *g[ vdS ].smoothPtr;
            coordD = *g[ vdT ].smoothPtr;
        }
        else{
            coordO = *g[ vdT ].smoothPtr;
            coordD = *g[ vdS ].smoothPtr;
        }
        double diffX = coordD.x() - coordO.x();
        double diffY = coordD.y() - coordO.y();
        double angle = atan2( diffY, diffX );

        g[ edge ].smoAngle = angle;
        g[ edge ].curAngle = angle;
    }

    // _metro->reorderID();

    unsigned int nVertices      = _metro->nStations();
    unsigned int nEdges         = _metro->nEdges();
    _d_Alpha                    = _metro->dAlpha();
    _targetLength                     = _metro->dBeta();

    cerr << "distance beta " << _targetLength << endl;

    // initialization
    _nVars = _nConstrs = 0;  
    _half_width = __half_width;
    _half_height = __half_height;

    // Total number of linear variables
    _nVars = 2 * nVertices;

    // Regular edge length 
    _nConstrs += 2 * nEdges;

    // Maximal angles of incident edges
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph )
    {
        DegreeSizeType degrees = out_degree( vertex, g );
        if( degrees > 1 ) _nConstrs += 2 * degrees;
    }

    // Positional constraints
    _nConstrs += 2 * nVertices;

    #ifdef ALONGEDGE
    _nConstrs += 2 * nVertices; 
    #endif

    #ifdef PARALLEL
    _nConstrs += 2 * nEdges; 
    #endif

    #ifdef OVERLAPPING
    _nConstrs += 2 * nVertices; 
    #endif

    // Decide if Path mode
    _pathMode = false;
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ) {
        if (g[edge].matchPath) {
            _pathMode = true;
            break;
        }
    }
    // _pathMode = true;
    _pathMode = false;


    // INIT CLOSEST POINTS
    _calculateAlongEdgeParameter();
    _handleHigherDegreeVertex();
    _removeSingleNonItersectingVertex();

    _initVars();
    _initCoefs();
    _initOutputs();
    _updateCoefs();
    _updateOutputs();
}

void Smooth::_calculateAlongEdgeParameter() {
    UndirectedGraph &g = _metro->g();
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
            Coord2 smoothCoord = *g[vertex].smoothPtr;
            EdgeDescriptor closestEdge = _guide->closestEdge(smoothCoord);
            Coord2 closestPoint = _guide->closestPointOnEdge(closestEdge, smoothCoord);

            int countIntersections = calculateIntersections(smoothCoord, closestPoint);
            g[vertex].intersectionsCount = countIntersections;
            g[vertex].smoothDistance = (closestPoint - smoothCoord).norm();
            g[vertex].initDistance = (closestPoint - smoothCoord).norm();
            g[vertex].closestPointOnEdge = closestPoint;

            if (countIntersections > 0) {
                g[vertex].intersection = true;
                if (_exclusivlyPathMode) {
                    g[vertex].smoothPath = g[vertex].autoPath;
                } else {
                    g[vertex].smoothPath = false;
                }

            } else {
                g[vertex].intersection = false;
                if (_exclusivlyPathMode) {
                    g[vertex].smoothPath = g[vertex].autoPath;
                }
                else {
                    g[vertex].smoothPath = true;
                }

            }
        }
}

void Smooth::_handleHigherDegreeVertex() {
    UndirectedGraph &g = _metro->g();
    // Handle nonIntersecting vertex with more than three neigbords
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
            DegreeSizeType degrees = out_degree( vertex, g );
            // look at all vertex with degree > 2 and non intersecting
            if (degrees > 2 && g[vertex].intersectionsCount < 1) {
                //cout neigborss of this vertex
                int countNonIntersectingNeighbors = 0;
                OutEdgeIterator e, e_end;

                for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                    EdgeDescriptor ed = *e;

                    VertexDescriptor vS = source(ed, g);
                    VertexDescriptor vT = target(ed, g);

                    if (g[vS].intersectionsCount < 1 && g[vT].intersectionsCount < 1) {
                        countNonIntersectingNeighbors++;
                    }
                }

                if (countNonIntersectingNeighbors > 2 ) {
                    // Look if there is a pair, remove rest
                    vector<unsigned int> lineIdCount;
                    int numbersOfEdgesWithSharedId = 0;
                    bool shareID = false;
                    unsigned int sharedID = 0;
                    for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                        EdgeDescriptor ed = *e;
                        VertexDescriptor vS = source(ed, g);
                        VertexDescriptor vT = target(ed, g);

                        if (g[vS].intersectionsCount < 1 && g[vT].intersectionsCount < 1) {
                            auto lineIDs = g[ed].lineID;
                            for (unsigned int line : lineIDs) {
                                if (find(lineIdCount.begin(), lineIdCount.end(), line) != lineIdCount.end()) {
                                    shareID = true;
                                    sharedID = line;
                                    numbersOfEdgesWithSharedId++;
                                } else {
                                    lineIdCount.push_back(line);
                                }
                            }
                        }
                    }
                    if (shareID && numbersOfEdgesWithSharedId <= 2 && false) { // there is a pair of edges that share the same line
                        int count = 0;
                        for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                            EdgeDescriptor ed = *e;
                            if (find(g[ed].lineID.begin(), g[ed].lineID.end(), sharedID) == g[ed].lineID.end() || count > 1) {
                                // remove vertex
                                VertexDescriptor vS = source(ed, g);
                                if (vS == vertex ) vS = target(ed, g);
                                g[vS].intersectionsCount = 1;
                                g[vS].intersection = true;
                                g[vS].smoothPath = false;
                            } else {
                                count++;
                            }

                        }
                    } else { // or remove all instead of the two best
                        map<VertexDescriptor, int> neighborDegree;
                        for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                            EdgeDescriptor ed = *e;

                            VertexDescriptor vS = source(ed, g);
                            VertexDescriptor vT = target(ed, g);

                            VertexDescriptor vj = vS;
                            if (vj == vertex) vj = vT;

                            OutEdgeIterator ej, ej_end;
                            int vj_count = 0;

                            // count neigbors of all neigbors
                            for ( tie( ej, ej_end ) = out_edges( vj, g ); ej != ej_end; ++ej ) {
                                EdgeDescriptor edj = *ej;
                                VertexDescriptor vjS = source(edj, g);
                                VertexDescriptor vjT = target(edj, g);
                                if (g[vjS].intersectionsCount < 1 && g[vjT].intersectionsCount < 1) {
                                    vj_count++;
                                }
                            }
                            neighborDegree.insert(pair<VertexDescriptor, int>(vj, vj_count));
                        }
                        // find the two neigbors with the two highest values
                        VertexDescriptor vd_1, vd_2;
                        int best = -1, secondBest = -1;

                        for (auto const&[vj, value] : neighborDegree) {
                            // cerr << g[vj].id << ", " << value << ", " << to_string(out_degree( vj, g )) << endl;
                            if (value > best ) {
                                vd_2 = vd_1;
                                secondBest = best;
                                vd_1 = vj;
                                best = value;
                            } else if (value > secondBest ) {
                                secondBest = value;
                                vd_2 = vj;
                            }
                        }
                        for (auto const&[vj, value] : neighborDegree) {
                            if (vj != vd_1 && vj != vd_2 ) {
                                g[vj].intersection = true;
                                g[vj].smoothPath = false;
                                g[vj].intersectionsCount = g[vj].intersectionsCount + 1;
                            }
                        }
                    }
                }
            }
        }
}

void Smooth::_removeSingleNonItersectingVertex() {
    UndirectedGraph &g = _metro->g();
    if (!_pathMode) {
        BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
            if (g[vertex].intersectionsCount == 0 ) {
                OutEdgeIterator e, e_end;
                bool remove = true;
                for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end && remove; ++e ) {
                    EdgeDescriptor ed = *e;
                    VertexDescriptor vS = source(ed, g);
                    VertexDescriptor vT = target(ed, g);
                    if (g[vS].intersectionsCount == 0 && g[vT].intersectionsCount == 0)
                        remove = false;
                }
                if (remove) {
                    g[vertex].intersectionsCount = g[vertex].intersectionsCount + 1;
                    g[vertex].intersection = true;
                    g[vertex].smoothPath = false;
                }
            }
        }
    }
}

/*
*   initialize the coefficient
*
*   Inputs  none
*   Outputs none
*/ 
void Smooth::_initCoefs( void )
{
    UndirectedGraph        & g            = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();
    unsigned int nVertices      = _metro->nStations();
    
    // initialization
    unsigned int nRows = 0;
    _coef.resize( _nConstrs, _nVars ); 
    _coef << Eigen::MatrixXd::Zero( _nConstrs, _nVars );
    
    // Regular edge length
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        unsigned int idS = g[ vdS ].id;
        unsigned int idT = g[ vdT ].id;

        double wL = _w_contextlength;

        _coef( nRows, idS ) = wL; 
        _coef( nRows, idT ) = -wL;
        nRows++;        
        // y
        _coef( nRows, idS + nVertices ) = wL;
        _coef( nRows, idT + nVertices ) = -wL;
        nRows++;
    }

    // Maximal angles of incident edges
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        DegreeSizeType degrees = out_degree( vertex, g );
        if( degrees > 1 ){

            // sort the embedding
            map< double, EdgeDescriptor > circM;
            OutEdgeIterator e, e_end;
            for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                EdgeDescriptor ed = *e;
                VertexDescriptor vS = source( ed, g );
                VertexDescriptor vT = target( ed, g );
                double angle = g[ ed ].geoAngle;

                if ( g[ vS ].id > g[ vT ].id ) {
                    if ( angle > 0 ) angle = -M_PI + g[ ed ].geoAngle;
                    else angle = M_PI + g[ ed ].geoAngle;
                }
                circM.insert( pair< double, EdgeDescriptor >( angle, ed ) );
            }

            // set coefficient
            double tanTheta = tan( (double)( degrees - 2 ) * M_PI / 2.0 / (double) degrees );
        
            map< double, EdgeDescriptor >::iterator itN = circM.begin();
            itN++;
            for ( map< double, EdgeDescriptor >::iterator it = circM.begin();
                  it != circM.end(); it++ ) {

                EdgeDescriptor edC = it->second;
                EdgeDescriptor edN = itN->second;
                VertexDescriptor vS = source( edC, g );
                VertexDescriptor vTC = target( edC, g );
                VertexDescriptor vTN = target( edN, g );
                unsigned int idS = g[ vS ].id;
                unsigned int idTC = g[ vTC ].id;
                unsigned int idTN = g[ vTN ].id;

                // x
                _coef( nRows, idS ) = _w_angle;   
                _coef( nRows, idTC ) = -0.5 * _w_angle;   
                _coef( nRows, idTN ) = -0.5 * _w_angle;   
                _coef( nRows, idTC + nVertices ) = -0.5 * _w_angle * tanTheta;   
                _coef( nRows, idTN + nVertices ) =  0.5 * _w_angle * tanTheta;   
                nRows++;

                // y
                _coef( nRows, idS + nVertices ) = _w_angle;   
                _coef( nRows, idTC + nVertices ) = -0.5 * _w_angle;   
                _coef( nRows, idTN + nVertices ) = -0.5 * _w_angle;   
                _coef( nRows, idTC ) =  0.5 * _w_angle * tanTheta;   
                _coef( nRows, idTN ) = -0.5 * _w_angle * tanTheta;   
                nRows++;

                itN++;
                if( itN == circM.end() ) {
                    itN = circM.begin();
                }
            }
        }   
    }
     
    // Positional constraints
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        unsigned int id = g[ vertex ].id;
        // x
        _coef( nRows, id ) = _w_position;
        nRows++;
        // y
        _coef( nRows, id + nVertices ) = _w_position;
        nRows++;
    }

   #ifdef ALONGEDGE 
   BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {

        double w = _w_alongEdge;
        DegreeSizeType degrees = out_degree( vertex, g );
        w *= degrees;
        // w *= pow(degrees * 0.5, 2);
        // w *= min(1.0, 1.0 / (minEdgeDistance * _constSmooth));
       if (_pathMode) {
           if (g[vertex].autoPath)// && g[vertex].intersectionsCount == 0)
               w *= 1;
           else
               w = 0;
       } else {
           if (g[vertex].smoothPath)
               w *= min(1.0, 1.0 / (g[vertex].initDistance * _constSmooth));
           else
               w = 0;
       }



        unsigned int id = g[ vertex ].id;
        _coef( nRows, id ) = w;
        nRows++;
        _coef( nRows, id + nVertices ) = w;
        nRows++;
    }
   #endif

   #ifdef PARALLEL
    BGL_FORALL_EDGES(edge, g, UndirectedGraph) {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        unsigned int idS = g[ vdS ].id;
        unsigned int idT = g[ vdT ].id;

        Coord2 vi = *g[ vdS ].geoPtr; 
        Coord2 vj = *g[ vdT ].geoPtr;
        Coord2 vji = vj - vi;
        Coord2 cji = (vj + vi); //Center of the edge 
        cji *= 0.5; 

        EdgeDescriptor closestEdge = _guide->closestEdge(cji); 
        Coord2 closestPoint = _guide->closestPointOnEdge(closestEdge, cji); 
        double dist = magnitude((closestPoint - cji)); 

        double wL = _w_parallel;
        wL *= min(1.0, 1.0 / (dist * _constSmooth));
        // x
        _coef( nRows, idS ) = wL; 
        _coef( nRows, idT ) = -wL;
        nRows++;        
        
        // y
        _coef( nRows, idS + nVertices ) = wL;
        _coef( nRows, idT + nVertices ) = -wL;
        nRows++;
    }
   #endif

    #ifdef OVERLAPPING 
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
        unsigned int id = g[vertex].id; 
        Coord2 vi = *g[vertex].geoPtr; 
        double w = 0; 
        Coord2 gi = _metro->closestPointOnGraph(vertex); 
        Coord2 d = vi - gi; 
        if( magnitude(d) < _gamma && false) {
            w = _w_overlapping; 
        }
        _coef( nRows, id ) = w; 
        nRows++; 
        _coef( nRows, id + nVertices ) = w; 
        nRows++; 
    }
    #endif
}

/*
* Resets smooth layout back to geo. layout
*
* Inputs    none
* Outputs   none
*/
void Smooth::_resetSmoothPtr() 
{
    UndirectedGraph        & g            = _metro->g();
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        Coord2 cord = *g[vertex].geoPtr;
        *g[vertex].smoothPtr = cord; 
    }

}

/*
* initialize the variables
*
* Inputs    none
* Outputs   none
*/
void Smooth::_initVars( void )
{
    UndirectedGraph        & g            = _metro->g();
    unsigned int nVertices      = _metro->nStations();

    _var.resize( _nVars );

    unsigned int nRows = 0;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        _var( nRows, 0 ) = g[ vertex ].smoothPtr->x();
        _var( nRows + nVertices, 0 ) = g[ vertex ].smoothPtr->y();
        nRows++;
    }
    assert( _nVars == 2*nRows );
}

/*
* initialize the output
*
* Inputs  none
* Outputs none
*/
void Smooth::_initOutputs( void )
{
    UndirectedGraph        & g            = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();

    // initialization
    unsigned int nRows = 0;
    _output.resize( _nConstrs );
    _output << Eigen::VectorXd::Zero( _nConstrs );
 
    // Regular edge length
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi = *g[ vdS ].geoPtr;
        Coord2 vj = *g[ vdT ].geoPtr;
        Coord2 vji = vi - vj;
        double targetLength = _targetLength;
        if (!g[vdS].isStation or !g[vdT].isStation)
            targetLength *= 0.5;
        double s = targetLength * g[ edge ].weight / vji.norm();
        double cosTheta = 1.0, sinTheta = 0.0;

        double wL = _w_contextlength;
        // if (!g[ vdS ].isStation or !g[vdT].isStation)
        //     wL *= 0.5;
        // x
        _output( nRows, 0 ) = wL * s * ( cosTheta * vji.x() - sinTheta * vji.y() );
        nRows++;
        // y
        _output( nRows, 0 ) = wL * s * ( sinTheta * vji.x() + cosTheta * vji.y() );
        nRows++;
    }

    // Maximal angles of incident edges
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){

        DegreeSizeType degrees = out_degree( vertex, g );

        if( degrees > 1 ){
            OutEdgeIterator eo_cur, eo_nxt, eo_end;
            for( tie( eo_cur, eo_end ) = out_edges( vertex, g ); eo_cur != eo_end; ++eo_cur ){

                // x
                _output( nRows, 0 ) = 0.0;
                nRows++;

                // y
                _output( nRows, 0 ) = 0.0;
                nRows++;
            }
        }
    }

    // Positional constraints
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        _output( nRows, 0 ) = _w_position * g[ vertex ].geoPtr->x();
        nRows++;
        _output( nRows, 0 ) = _w_position * g[ vertex ].geoPtr->y();
        nRows++;
    }

    #ifdef ALONGEDGE    
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
        Coord2 cord = *g[vertex].geoPtr; 
        EdgeDescriptor closestEdge = EdgeDescriptor();
        double minEdgeDistance = 10e+10; 
        bool first = true; 

        //find for each vertex the closest Edge 
        BGL_FORALL_EDGES(edge, gGuide, UndirectedGraph) {
            Coord2 sourceCord = *gGuide[ source(edge, gGuide) ].coordPtr; //xy1
            Coord2 targetCord = *gGuide[ target(edge, gGuide) ].coordPtr;  //xy2
            
            Coord2 closePoint = closestPointOnEdge(sourceCord, targetCord, cord); 
            double distance = magnitude(closePoint - cord); 
             
            if (distance < minEdgeDistance || first) {
                minEdgeDistance = distance; 
                closestEdge = edge; 
                first = false; 
            }
        }
        Coord2 gs = *gGuide[ source(closestEdge, gGuide) ].coordPtr;
        Coord2 gt = *gGuide[ target(closestEdge, gGuide) ].coordPtr; 
        
        Coord2 v_init; 
        if ( g[vertex].metroShape ) v_init = closestPointOnEdge(gs, gt, cord); 
        else v_init = cord; 
        double w = _w_alongEdge; 
        DegreeSizeType degrees = out_degree( vertex, g );
        w *= degrees;
        // w *= pow(degrees * 0.5, 2);
        // w *= min(1.0, 1.0 / (minEdgeDistance * _constSmooth));
        if (_pathMode) {
            if (g[vertex].autoPath )//&& g[vertex].intersectionsCount == 0)
                w *= 1;
            else
                w = 0;
        } else {
            if (g[vertex].smoothPath) // TODO metro shape adden
                w *= min(1.0, 1.0 / (g[vertex].initDistance * _constSmooth));
            else
                w = 0;
        }


        _output(nRows, 0) =  v_init[0] * w; 
        nRows++; 
        _output(nRows, 0) =  v_init[1] * w;  
        nRows++;
    }
    #endif

    #ifdef PARALLEL
    BGL_FORALL_EDGES(edge, g, UndirectedGraph) {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi = *g[ vdS ].geoPtr;
        Coord2 vj = *g[ vdT ].geoPtr;
        Coord2 vji = vj - vi;
        Coord2 cji = (vj + vi); //Center of the edge 
        cji *= 0.5; 

        EdgeDescriptor closeEdgeDiscriptor = _guide->closestEdge(cji);
        Coord2 guideSource = *gGuide[ source( closeEdgeDiscriptor, gGuide) ].coordPtr; 
        Coord2 guideTarget = *gGuide[ target( closeEdgeDiscriptor, gGuide) ].coordPtr; 

        Coord2 guide = guideTarget - guideSource; 
        double guideTheta = atan2(guide.y(), guide.x()); 
        double edgeTheta = atan2(vji.y(), vji.x()); 
        double theta = guideTheta - edgeTheta; 
        
        #ifdef SHORTESTANGLE
        if (theta > M_PI /2 && theta < M_PI ) {
            theta = ( M_PI - theta ) * -1; 
        } else if (theta < -M_PI / 2 && theta > - M_PI ) {
            theta = (- M_PI - theta ) * -1; 
        }
        #endif

        vi = vi - cji; 
        vj = vj - cji; 
        if (g[vdS].metroShape || g[vdT].metroShape ) {
            vi.setX( vi.x() * cos(theta) - vi.y() * sin(theta)); 
            vi.setY( vi.x() * sin(theta) + vi.y() * cos(theta));   
            vj.setX( vj.x() * cos(theta) - vj.y() * sin(theta)); 
            vj.setY( vj.x() * sin(theta) + vi.y() * cos(theta)); 
        }

        vi = vi + cji; 
        vj = vj + cji; 


        Coord2 closestPoint = _guide->closestPointOnEdge(closeEdgeDiscriptor, cji); 
        double dist = magnitude((closestPoint - cji)); 
        double wL = _w_parallel;
        wL *= min(1.0, 1.0 / (dist * _constSmooth)); 

        _output( nRows, 0 ) = wL * (vi.x() - vj.x());
        nRows++;
        // y
        _output( nRows, 0 ) = wL * (vi.y() - vj.y());
        nRows++;
    }
    #endif

    #ifdef OVERLAPPING
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
        Coord2 vi = *g[vertex].geoPtr;
        double w = 0; 
        Coord2 dir = Coord2(0,0); 
        Coord2 gi = _metro->closestPointOnGraph(vertex); 
        Coord2 d = vi - gi; 
        if (magnitude(d) < _gamma && false) {
            w = _w_overlapping; 
            dir = d; 
        } 

        _output(nRows, 0) = w * (vi.x() + dir.x()); 
        nRows++; 
        _output(nRows, 0) = w * (vi.y() + dir.y()); 
        nRows++; 
    }

    #endif // Overlapping 
}

/*
* update the coefs
*
* Inputs    none
* Outputs   none
*/
void Smooth::_updateCoefs( void )
{
    UndirectedGraph               & g             = _metro->g();
    unsigned int        nVertices       = _metro->nStations();
    unsigned int        nEdges          = _metro->nEdges();
    unsigned int        nVE             = 0;
    unsigned int        nB              = 0;
    vector< double >    ratioR          = _metro->ratioR();
    unsigned int nRows                  = 0; 
    Eigen::MatrixXd     oldCoef;
    
    oldCoef = _coef.block( 0, 0, _nConstrs, _nVars );

    #ifdef ALONGEDGE
    nRows = _nConstrs - _nVars - (2 * nEdges) - _nVars;
        BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
                double w = _w_alongEdge;
                DegreeSizeType degrees = out_degree( vertex, g );
                w *= degrees;
                if (_pathMode) {
                    if (g[vertex].autoPath)// && g[vertex].intersectionsCount == 0)
                        w *= 1;
                    else
                        w = 0;
                } else {
                    if (g[vertex].smoothPath) // TODO metro shape adden
                        w *= min(1.0, 1.0 / (g[vertex].initDistance * _constSmooth));
                    else
                        w = 0;
                }

                unsigned int id = g[ vertex ].id;
                oldCoef( nRows, id ) = w;
                nRows++;
                oldCoef( nRows, id + nVertices ) = w;
                nRows++;
            }
    #endif

    #ifdef OVERLAPPING
    nRows = _nConstrs - _nVars; 
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
        unsigned int id = g[ vertex ].id; 
        Coord2 vi = *g[ vertex ].smoothPtr;
        double w = 0; 
        Coord2 gi = _metro->closestPointOnGraph( vertex ); 
        Coord2 d = vi - gi; 
        if ( magnitude(d) < _gamma && false) {
            w = _w_overlapping; 
        } 
        oldCoef( nRows, id ) = w; 
        nRows++; 
        oldCoef ( nRows, id + nVertices ) = w; 
        nRows++; 
    }
    #endif // OVERLAPPINg

#ifdef  SMOOTH_CONFLICT
    nVE = _metro->VEconflict().size();
#endif  // SMOOTH_CONFLICT
#ifdef  SMOOTH_BOUNDARY
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph )
    {
        double minD = _targetLength / 2.0;

        if ( g[ vd ].smoothPtr->x() <= -( _half_width - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->x() >= ( _half_width - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->y() <= -( _half_height - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->y() >= ( _half_height - minD ) ) nB++;
    }
#endif  // SMOOTH_BOUNDARY
    _coef.resize( _nConstrs + nB + 2*nVE, _nVars );
    // copy old coefficient
    _coef << oldCoef, Eigen::MatrixXd::Zero( nB + 2*nVE, _nVars );
    nRows = _nConstrs;

#ifdef  SMOOTH_BOUNDARY
    // add boundary coefficient
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph ){

        unsigned int id = g[ vd ].id;
        double minD = _targetLength / 2.0;

        if ( g[ vd ].smoothPtr->x() <= -( _half_width - minD ) ) {
            _coef( nRows, id ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].smoothPtr->x() >= ( _half_width - minD ) ) {
            _coef( nRows, id ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].smoothPtr->y() <= -( _half_height - minD ) ) {
            _coef( nRows, id + nVertices ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].smoothPtr->y() >= ( _half_height - minD ) ) {
            _coef( nRows, id + nVertices ) = _w_boundary;
            nRows++;
        }
    }
#endif  // SMOOTH_BOUNDARY


#ifdef  SMOOTH_CONFLICT
    // add conflict coefficient
    unsigned int countVE = 0;
    for ( VEMap::iterator it = _metro->VEconflict().begin();
          it != _metro->VEconflict().end(); ++it ) {
        VertexDescriptor vdV = it->second.first;
        EdgeDescriptor ed = it->second.second;
        VertexDescriptor vdS = source( ed, g );
        VertexDescriptor vdT = target( ed, g );
        unsigned int idV = g[ vdV ].id;
        unsigned int idS = g[ vdS ].id;
        unsigned int idT = g[ vdT ].id;
        double r = ratioR[ countVE ];

        // x
        _coef( nRows, idV ) = 1.0 * _w_crossing;
        _coef( nRows, idS ) = -r * _w_crossing;
        _coef( nRows, idT ) = ( r - 1.0 ) * _w_crossing;
        nRows++;

        // y
        _coef( nRows, idV + nVertices ) = 1.0 * _w_crossing;
        _coef( nRows, idS + nVertices ) = -r * _w_crossing;
        _coef( nRows, idT + nVertices ) = ( r - 1.0 ) * _w_crossing;
        nRows++;
        
        countVE++;
    } 
#endif  // SMOOTH_CONFLICT

}

/*
* update the output
*
* Inputs    none
* Outputs   none
*/
void Smooth::_updateOutputs( void ) // TODO Check all divide zeros
{
    UndirectedGraph               & g             = _metro->g();
    UndirectedGraph               & gGuide        = _guide->g();
    unsigned int        nVertices       = _metro->nStations();
    unsigned int        nVE             = 0;
    unsigned int        nB              = 0;
    vector< double >    ratioR          = _metro->ratioR();
    vector< double >    ratioGeo        = _metro->ratioGeo();

    unsigned int nRows = 0;
    Eigen::VectorXd     oldOutput;
    oldOutput = _output;
#ifdef  SMOOTH_CONFLICT
    nVE = _metro->VEconflict().size();
#endif  // SMOOTH_CONFLICT
#ifdef  SMOOTH_BOUNDARY
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph )
    {
        double minD = _targetLength / 2.0;
        if ( g[ vd ].smoothPtr->x() <= -( _half_width - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->x() >= ( _half_width - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->y() <= -( _half_height - minD ) ) nB++;
        if ( g[ vd ].smoothPtr->y() >= ( _half_height - minD ) ) nB++;
    }
#endif  // SMOOTH_BOUNDARY
    _output.resize( _nConstrs + nB + 2*nVE );
    _output << Eigen::VectorXd::Zero( _nConstrs + nB + 2*nVE );

    // Regular edge length
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi = *g[ vdS ].geoPtr;
        Coord2 vj = *g[ vdT ].geoPtr;
        Coord2 vji = vi - vj;
        Coord2 pi = *g[ vdS ].smoothPtr;
        Coord2 pj = *g[ vdT ].smoothPtr;
        Coord2 pji = pi - pj;
        // double vjiNorm = vji.norm();
        double targetLength = _targetLength;
        if (!g[vdS].isStation or !g[vdT].isStation)
            targetLength *= 0.5;

        double s = targetLength * g[ edge ].weight / vji.norm();
        double angleV = atan2( vji.y(), vji.x() );
        double angleP = atan2( pji.y(), pji.x() );
        double diffTheta = angleP - angleV;
        double cosTheta = cos( diffTheta ); 
        double sinTheta = sin( diffTheta );

        double wL = _w_contextlength;
        // if (!g[ vdS ].isStation or !g[vdT].isStation)
        //     wL *= 0.5;

        // x
        _output( nRows, 0 ) = wL * s * ( cosTheta * vji.x() - sinTheta * vji.y() );
        nRows++;

        // y
        _output( nRows, 0 ) = wL * s * ( sinTheta * vji.x() + cosTheta * vji.y() );
        nRows++;
    }


    // Maximal angles of incident edges
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){

        DegreeSizeType degrees = out_degree( vertex, g );

        if( degrees > 1 ){

            double tanTheta = tan( (double)( degrees - 2 ) * M_PI / 2.0 / (double) degrees );
            OutEdgeIterator eo_cur, eo_nxt, eo_end;
            for( tie( eo_cur, eo_end ) = out_edges( vertex, g ); eo_cur != eo_end; ++eo_cur ){

                // x
                _output( nRows, 0 ) = oldOutput( nRows, 0 );
                nRows++;

                // y
                _output( nRows, 0 ) = oldOutput( nRows, 0 );
                nRows++;
            }
        }
    }

    // Positional constraints
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){

        // x
        _output( nRows, 0 ) = oldOutput( nRows, 0 );
        nRows++;

        // y
        _output( nRows, 0 ) = oldOutput( nRows, 0 );
        nRows++;
    }

    #ifdef ALONGEDGE
    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
        Coord2 cord = *g[vertex].geoPtr;
        Coord2 scord = *g[vertex].smoothPtr; 
        // EdgeDescriptor closestEdge = EdgeDescriptor();
        double minEdgeDistance;

        EdgeDescriptor closestEdge = _guide->closestEdge(scord);
        Coord2 closestPoint = g[vertex].closestPointOnEdge;
        Coord2 closestPointOrig = _guide->closestPointOnEdge(_guide->closestEdge(cord), cord);
        int currentIntersections = calculateIntersections(scord, closestPoint);
        Coord2 v_updated = closestPoint;

        int countIntersections = 0;

        if ( g[vertex].smoothPath ) {
            countIntersections = g[vertex].intersectionsCount;
        }
        else {
            v_updated = scord;
        }
        if (_pathMode) {
            if (!g[vertex].autoPath ) v_updated = scord;
            if (currentIntersections > 0 ) v_updated = scord; // TODO not sure if remove
        } else {
            if (currentIntersections > 0) {
                v_updated = scord;
            }
            if (countIntersections > 0 ) {
                v_updated = scord; // TOD0
            }
        }


        double w = _w_alongEdge; 
        DegreeSizeType degrees = out_degree( vertex, g );
        w *= degrees;
        // w *= min(1.0, 1.0 / (minEdgeDistance * _constSmooth)); // TODO other side of the equation has to be updated as well
        if (_pathMode) {
            if (g[vertex].autoPath)
                w *= 1;
            else
                w = 0;
        } else {
            if (g[vertex].smoothPath) { // TODO metro shape adden
                w *= min(1.0, 1.0 / (g[vertex].initDistance * _constSmooth));
            } else
                w = 0;
        }

        _output(nRows, 0) =  v_updated.getX() * w;
        nRows++; 
        _output(nRows, 0) =  v_updated.getY() * w; 
        nRows++;  
    }
    #endif

    #ifdef PARALLEL
    BGL_FORALL_EDGES(edge, g, UndirectedGraph) {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi = *g[ vdS ].smoothPtr;
        Coord2 vj = *g[ vdT ].smoothPtr;
        Coord2 vji = vj - vi;
        Coord2 cji = vji;
        cji /= 2;  
        cji += vi; 

        Coord2 vi_geo = *g[vdS].geoPtr; 
        Coord2 vj_geo = *g[vdT].geoPtr; 
        Coord2 vji_geo = vj_geo - vi_geo; 
 
        Coord2 cji_geo = vji_geo; 
        cji_geo /= 2; 
        cji_geo += vi_geo; 

        EdgeDescriptor closeEdgeDiscriptor = _guide->closestEdge(cji);
        Coord2 guideSource = *gGuide[ source( closeEdgeDiscriptor, gGuide) ].coordPtr; 
        Coord2 guideTarget = *gGuide[ target( closeEdgeDiscriptor, gGuide) ].coordPtr; 

        Coord2 guide = guideTarget - guideSource; 
        double guideAngle = atan2(guide.y(), guide.x()); 
        double curAngle = atan2(vji.y(), vji.x()); 
        
        if ( guideAngle < 0.0 && curAngle >= 0.0 ) guideAngle += M_PI; 
        else if ( guideAngle > 0.0 && curAngle < 0.0 ) guideAngle -= M_PI; 

        double theta = guideAngle - curAngle; 
        if (fabs(theta) > M_PI_4 && fabs(theta) < M_PI * 0.75) {
            if (theta < 0 ) theta = - M_PI_2; 
            else theta = M_PI_2; 
        } 

        vi = vi - cji; 
        vj = vj - cji; 

        if (g[vdS].metroShape || g[vdT].metroShape ) {
            vi.setX( vi.x() * cos(theta) - vi.y() * sin(theta)); 
            vi.setY( vi.x() * sin(theta) + vi.y() * cos(theta));   
            vj.setX( vj.x() * cos(theta) - vj.y() * sin(theta)); 
            vj.setY( vj.x() * sin(theta) + vi.y() * cos(theta)); 
        }

        vi = vi + cji; 
        vj = vj + cji; 

        Coord2 closestPoint = _guide->closestPointOnEdge(closeEdgeDiscriptor, cji_geo); 
        double dist = (closestPoint - cji_geo).norm();
        double wL = _w_parallel;
        wL *= min(1.0, 1.0 / (dist * _constSmooth));

        _output( nRows, 0 ) = wL * (vi.x() - vj.x());
        nRows++;
        // y
        _output( nRows, 0 ) = wL * (vi.y() - vj.y());
        nRows++;
    }
    #endif

    #ifdef OVERLAPPING 
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ) { 
        Coord2 vi = *g[vertex].smoothPtr; 
        double w = 0; 
        Coord2 dir = Coord2(0,0); 
        Coord2 gi = _metro->closestPointOnGraph(vertex); 
        Coord2 d = vi - gi; 
        if (magnitude(d) < _gamma && false ) {
            w = _w_overlapping; 
            dir = d; 
        }
        double magDir = magnitude(dir); 
        if (magDir > 0) {
            double r = _gamma / magDir; 
            dir *= r; 
        }
        _output(nRows, 0) = w * (vi.x() + dir.x());
        nRows++; 
        _output(nRows, 0) = w * (vi.y() + dir.y()); 
        nRows++; 
    }
    #endif // OVERLAPPING

#ifdef  SMOOTH_BOUNDARY
    // boundary constraints
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph ){

        double minD = _targetLength / 2.0;
        if ( g[ vd ].smoothPtr->x() <= -( _half_width - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * -( _half_width - minD );
            nRows++;
        }
        if ( g[ vd ].smoothPtr->x() >= ( _half_width - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * ( _half_width - minD );
            nRows++;
        }
        if ( g[ vd ].smoothPtr->y() <= -( _half_height - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * -( _half_height - minD );
            nRows++;
        }
        if ( g[ vd ].smoothPtr->y() >= ( _half_height - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * ( _half_height - minD );
            nRows++;
        }
    }
#endif  // SMOOTH_BOUNDARY

#ifdef  SMOOTH_CONFLICT
    unsigned int countVE = 0;
    for ( VEMap::iterator it = _metro->VEconflict().begin();
          it != _metro->VEconflict().end(); ++it ) {
        VertexDescriptor vdV = it->second.first;
        EdgeDescriptor ed = it->second.second;
        VertexDescriptor vdS = source( ed, g );
        VertexDescriptor vdT = target( ed, g );
        double r = ratioR[ countVE ];
        double rGeo = ratioGeo[ countVE ];
        r = rGeo; // zero issue but i think the right one.

        Coord2 vCur = *g[ vdV ].geoPtr; // geo ptr?
        Coord2 pCur = r * *g[ vdS ].geoPtr + ( 1.0-r ) * *g[ vdT ].geoPtr;
        Coord2 v = *g[vdV].geoPtr;
        Coord2 p = r * *g[ vdS ].geoPtr + ( 1.0-r ) * *g[ vdT ].geoPtr;

        double minD = _targetLength;
        if (!g[vdV].isStation or !g[vdV].isStation) {
            minD *= 0.5;
        }

        double delta;
        if (( vCur - pCur ).norm() != 0.0) delta = minD / ( vCur - pCur ).norm(); // multiply to get target distance
        else {  // v lies on edge –> v = p; move v in the direction of the longest adeazent edge
            delta = minD;
            v = v + (*g[vdV].smoothPtr - v).normalize();
        }
        // x
        _output( nRows, 0 ) = _w_crossing * delta * ( v - p ).x();
        nRows++;

        // y
        _output( nRows, 0 ) = _w_crossing * delta * ( v - p ).y();
        nRows++;

        countVE++;
    }
#endif  // SMOOTH_CONFLICT
}

int Smooth::calculateIntersections( Coord2 vertexCord, Coord2 closestPointOnGuide ) {
    UndirectedGraph &gMetro = _metro->g();
    int countIntersections = 0;

    double x1 = vertexCord.x();
    double y1 = vertexCord.y();
    double x2 = closestPointOnGuide.x();
    double y2 = closestPointOnGuide.y();

    double A1 = y2 - y1;
    double B1 = x1 - x2;
    double C1 = A1 * x1 +B1 * y1;

    // reflected
    double x1_ = closestPointOnGuide.x();
    double y1_ = closestPointOnGuide.y();
    double x2_ = closestPointOnGuide.x() + (closestPointOnGuide.x() - vertexCord.x());
    double y2_ = closestPointOnGuide.y() + (closestPointOnGuide.y() - vertexCord.y());

    double A1_ = y2_ - y1_;
    double B1_ = x1_ - x2_;
    double C1_ = A1_ * x1_ +B1_ * y1_;

    BGL_FORALL_EDGES(edge, gMetro, UndirectedGraph) {
        Coord2 edgeSource = *gMetro[source(edge, gMetro)].smoothPtr;
        Coord2 edgeTarget = *gMetro[target(edge, gMetro)].smoothPtr;

        if (edgeSource.x() != vertexCord.x() && edgeTarget.x() != vertexCord.x()) {
            double x3 = edgeSource.x();
            double y3 = edgeSource.y();
            double x4 = edgeTarget.x();
            double y4 = edgeTarget.y();

            double A2 = y4 - y3;
            double B2 = x3 - x4;
            double C2 = A2 * x3 + B2 * y3;

            double x,y;
            double x_, y_;
            double det = A1 * B2 - A2 * B1;
            double det_ = A1_ * B2 - A2 * B1_;
            if ( fabs(det) < 0.001 || fabs(det_) < 0.001 ) {
                //Lines are parallel
                x = y = 0;
                x_ = y_ = 0;
            } else {
                x = (B2 * C1 - B1 * C2) / det;
                y = (A1 * C2 - A2 * C1) / det;
                x_ = (B2 * C1_ - B1_ * C2) / det_;
                y_ = (A1_ * C2 - A2 * C1_) / det_;
                // check if it lays on ray
                double tol = 0.01;
                tol = MIN_DISTANCE;
                if (    (min(x1,x2) - tol <= x && x <= max(x1, x2) + tol ) &&
                        (min(x3,x4) - tol <= x && x <= max(x3, x4) + tol ) &&
                        (min(y1,y2) - tol <= y && y <= max(y1, y2) + tol ) &&
                        (min(y3,y4) - tol <= y && y <= max(y3, y4) + tol)
                       ) {
                    // "is on edge");
                    countIntersections++;
                }
                else if (   (min(x1_,x2_) - tol <= x_ && x_ <= max(x1_, x2_) + tol ) &&      // check reflected ray
                            (min(x3 ,x4 ) - tol <= x_ && x_ <= max(x3 , x4 ) + tol ) &&
                            (min(y1_,y2_) - tol <= y_ && y_ <= max(y1_, y2_) + tol ) &&
                            (min(y3 ,y4 ) - tol <= y_ && y_ <= max(y3 , y4 ) + tol )
                        )  {
                    countIntersections++;
                } else {
                    // println("is not on edge");
                }
            }
        }
    }
    return countIntersections;
}

/*
* Calculates magnitude of a vector from [0,0] to v
*
* Inputs    Coordines of point
* Outputs   Magnitude
*/
double Smooth::magnitude( Coord2 v ) {
    double d = v.x() * v.x() + v.y() * v.y(); 
    return sqrt(d);
}

/*
* Calculates the closest point of a point laying on a edge
* Inputs    source Coord. of edge
*           target Coord. of edge
* 
* Outputs   Coord. of point 
*/ 
Coord2 Smooth::closestPointOnEdge( Coord2 source, Coord2 target, Coord2 vertex ) {
    Eigen::RowVectorXd se(2); 
    Eigen::RowVectorXd te(2); 
    Eigen::RowVectorXd v(2);
    Eigen::RowVectorXd e(2); 
    Eigen::RowVectorXd sv(2);
    Eigen::RowVectorXd closestPoint(2); 
    se << source.x(), source.y(); 
    te << target.x(), target.y(); 
    v << vertex.x(), vertex.y();
    e = te - se; 
    e.normalize();   
    sv = v - se; 
    closestPoint = sv.dot(e) * e;
    closestPoint = se + closestPoint;  

    if (closestPoint[0] < min(se[0], te[0])) {
        if (se[0] < te[0]) {
            closestPoint = se; 
        } else {
            closestPoint = te; 
        }
    } else if (closestPoint[0] > max(se[0], te[0])) {
        if (se[0] > te[0] ) closestPoint = se; 
        else closestPoint = te; 
    }
    
    Coord2 result; 
    result.set(closestPoint[0], closestPoint[1]);  
    return result; 
}

/*
* LeastSquare optimization
*
* Inputs    none
* Outputs   none
*/
double Smooth::LeastSquare( unsigned int iter )
{
    double mse = 0.0;
    for( unsigned int i = 0; i < iter; i++ ) {

        Eigen::VectorXd last_var = _var;

        // ### optimization ###
        _var = ( ( _coef.transpose() * _coef ) ).inverse() * _coef.transpose() * _output;

        // ### retrieve the result ###
        retrieve();

        // ### update coefficient matrix ###
        _updateCoefs();

        // ### update target values ###
        _updateOutputs();

        // node movement 
        Eigen::VectorXd err = last_var - _var;
        mse = err.adjoint() * err;

        if( ( mse ) < 1.0e-1 ) break;
    }
    return mse;
}

/*
* ConjugateGradient optimization
*
* Inputs    none
* Outputs   none
*/
double Smooth::ConjugateGradient( unsigned int iter )
{
    // initialization, prepare the square matrix
    Eigen::MatrixXd A;
    Eigen::VectorXd b, Ap;
    A = _coef.transpose() * _coef;
    b = _coef.transpose() * _output;

    // initialization
    Eigen::VectorXd err = b - A * _var;    
    Eigen::VectorXd p = err;  

    // main algorithm
    for( int i = 0; i < (int) (iter * 1.0); i++ ) {  // TODO reset to (int) iter
        A = _coef.transpose() * _coef;
        b = _coef.transpose() * _output;
        Ap = A * p;

        double alpha = (double)( p.transpose() * err ) / (double)( p.transpose() * Ap );
        _var = _var + alpha * p; 
        err = b - A * _var;

        if ( sqrt( err.adjoint() * err ) < 1e-10 ) {   
            cerr << "sqrterror(" << i << ") = " << sqrt( err.adjoint() * err ) << endl;
            break;
        }
        else {
            double beta = -1.0 * (double)( err.transpose() * Ap ) / (double)( p.transpose() * Ap );
            p = err + beta * p;
        }
        // show progress in console
        if (i % 25 == 0) {
            cout << "Smooth ... " << i << " / " << iter << " ... " << endl;
        }
        _calculateAlongEdgeParameter();
        if (i == int(iter * 0.25) && !_exclusivlyPathMode && false) {
            _pathMode = false;
        } else if (i > iter * 0.25 and true ) {
            if (!_guide->isComplex())
               _handleHigherDegreeVertex();
            _removeSingleNonItersectingVertex();
        }
        retrieve();
        _updateCoefs();
        _updateOutputs();
    }
    _calculateAlongEdgeParameter();
    if (!_guide->isComplex())
       _handleHigherDegreeVertex();


    _removeSingleNonItersectingVertex();
    _asignSmoothEdge();
    return sqrt( err.adjoint() * err );
}

void Smooth::_asignSmoothEdge() {
    UndirectedGraph        & graph = _metro->g();
    BGL_FORALL_EDGES(edge, graph,UndirectedGraph) {
        VertexDescriptor vdS = source(edge, graph);
        VertexDescriptor vdT = target(edge, graph);

        if (graph[vdS].smoothPath && graph[vdT].smoothPath) {
            graph[edge].isSmoothPath = true;
        } else {
            graph[edge].isSmoothPath = false;
        }
    }
}

/*
* retrieve the result
*
* Inputs    none
* Outputs   none
*/
void Smooth::retrieve( void )
{
    UndirectedGraph        & g            = _metro->g();
    unsigned int nVertices      = _metro->nStations();

    vector< vector< VertexDescriptor > > vdMatrix;
    // find the vertex that is too close to an edge
    for ( VEMap::iterator it = _metro->VEconflict().begin();
          it != _metro->VEconflict().end(); ++it ) {

        vector< VertexDescriptor > vdVec;
        VertexDescriptor vdV = it->second.first;
        EdgeDescriptor ed = it->second.second;
        VertexDescriptor vdS = source( ed, g );
        VertexDescriptor vdT = target( ed, g );
        vdVec.push_back( vdV );
        vdVec.push_back( vdS );
        vdVec.push_back( vdT );

        vdMatrix.push_back( vdVec );
    }
    unsigned int nRows = 0;
    double scale = 1.0;

    nRows = 0;

    // vector<Coord2> oldPositions;
    map<VertexDescriptor, Coord2> oldPositions;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        Coord2 curCoord = *g[ vertex ].smoothPtr;
        Coord2 newCoord = Coord2( _var( nRows, 0 ), _var( nRows + nVertices, 0 ) );
        if (isnan(newCoord.y())) {
            // newCoord = curCoord;
        }
        // oldPositions.push_back( *g[ vertex ].coordPtr );
        oldPositions[vertex] = *g[ vertex ].coordPtr;
        Coord2 newClosPoint = _metro->closestPointOnGraph(vertex, newCoord);
        double distance = (newCoord - newClosPoint).norm();
        if (distance < MIN_DISTANCE) {
            g[ vertex ].coordPtr->x() = g[ vertex ].smoothPtr->x() = curCoord.x() * .5 + _var( nRows, 0 ) * .5;
            g[ vertex ].coordPtr->y() = g[ vertex ].smoothPtr->y() = curCoord.y() * .5 + _var( nRows + nVertices, 0 ) * .5;
            g[vertex].collision = true;
        } else {
            g[ vertex ].coordPtr->x() = g[ vertex ].smoothPtr->x() = newCoord.x();
            g[ vertex ].coordPtr->y() = g[ vertex ].smoothPtr->y() = newCoord.y();
            g[vertex].collision = false;
            assert( vertexID[ vertex ] == nRows );
        }
        nRows++;
    }

    // Check if violates Planarity
    int count = 0;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        bool planarViolation = false;
        OutEdgeIterator e, e_end;
        for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
            EdgeDescriptor ed = *e;
            auto intersectingEdges = _metro->intersectingEdge(ed);
            if (intersectingEdges.size() > 0){
                planarViolation = true;
                for (auto edInt : intersectingEdges) {
                    auto sD = source(edInt.first, g);
                    auto tD = target(edInt.first, g);
                    g[ sD ].coordPtr->x() = g[ sD ].smoothPtr->x() = oldPositions[sD].x();
                    g[ tD ].coordPtr->x() = g[ tD ].smoothPtr->x() = oldPositions[tD].x();
                }
            }
        }
        if (planarViolation) {
            // g[ vertex ].coordPtr->x() = g[ vertex ].smoothPtr->x() = g[ vertex ].coordPtr->x() * 0.5 + oldPositions[count].x() * 0.5;
            g[ vertex ].coordPtr->x() = g[ vertex ].smoothPtr->x() = oldPositions[vertex].x();
            // g[ vertex ].coordPtr->y() = g[ vertex ].smoothPtr->y() = g[ vertex ].coordPtr->y() * 0.5 + oldPositions[count].y() * 0.5;
            g[ vertex ].coordPtr->y() = g[ vertex ].smoothPtr->y() = oldPositions[vertex].y();
        }



        count++;
    }

    // update smooth angle
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ){

        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );

        Coord2 coordO;
        Coord2 coordD;
        if( g[ vdS ].id < g[ vdT ].id ){
            coordO = *g[ vdS ].smoothPtr;
            coordD = *g[ vdT ].smoothPtr;
        }
        else{
            coordO = *g[ vdT ].smoothPtr;
            coordD = *g[ vdS ].smoothPtr;
        }
        double diffX = coordD.x() - coordO.x();
        double diffY = coordD.y() - coordO.y();
        double angle = atan2( diffY, diffX );

        g[ edge ].smoAngle = angle;
        g[ edge ].curAngle = angle;
    }

    // check possible conflict
    _metro->checkVEConflicts();
    // UPDATE CLOSEST POINTS
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        Coord2 smoothCoord = *g[vertex].smoothPtr; // why sometimes nan?
        EdgeDescriptor closestEdge = _guide->closestEdge(smoothCoord);
        Coord2 closestPoint = _guide->closestPointOnEdge(closestEdge, smoothCoord);
        g[vertex].smoothDistance = (closestPoint - smoothCoord).norm();
        g[vertex].closestPointOnEdge = closestPoint;
    }
}


/* 
* memory management
*
* Inputs    none
* Outputs   none
*/
void Smooth::clear( void ) {}

//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------

/*
* default constructor
*
* Inputs  none
* Outputs none
*/
Smooth::Smooth( void )
{   
}

/*
* copy constructor
*
* Inputs obj : object of this class
* Outputs none
*/
Smooth::Smooth( const Smooth & obj ) {}

/*
* destructor
*
* Inputs    none
* Outputs   none
*/
Smooth::~Smooth( void ) {}

