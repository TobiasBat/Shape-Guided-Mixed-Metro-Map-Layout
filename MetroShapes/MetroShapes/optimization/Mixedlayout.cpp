/*
*	Mixedlayout.cpp
*   Based on FocusContext/Octilinear.cpp
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   To calculate the Mixed layout
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstring>
#include <cmath>

using namespace std;

#include "Mixedlayout.h"

//------------------------------------------------------------------------------
//	Protected functions
//------------------------------------------------------------------------------

/*
* initialize the constrained optimization problem
*
* Inputs
*     __metro : pointer to metro
*
* Outputs
*     none
*/
void Mixedlayout::_init( Metro * __metro, Guide * __guid ,double __half_width, double __half_height )
{
    _metro                      = __metro;
    // _metro->re
    unsigned int nVertices      = _metro->nStations();
    unsigned int nEdges         = _metro->nEdges();
    _d_Alpha                    = _metro->dAlpha();
    _d_Beta                     = _metro->dBeta();
    _guide                       = __guid; 

    _updateEdgeCurAngle();
    _metro->checkVEConflicts();
    _metro->reorderID();
    _setTargetAngle2();

#ifdef DEGREE_TWO_HEURISTIC
    _metro->simplifyMixedLayout();
    _metro->reorderID();
    _updateEdgeCurAngle();
    _setTargetAngle2();

    _metro->reorderID();
    _metro->checkVEConflicts();
#endif

    nVertices      = _metro->nStations();
    nEdges         = _metro->nEdges();

    _nVars = _nConstrs = 0;
    _half_width = __half_width;
    _half_height = __half_height;

    //  Total number of linear variables
    _nVars = 2 * _metro->nStations();

    // Regular edge length 
    _nConstrs += 2 * nEdges;

    // Positional constraints
    _nConstrs += 2 * nVertices;

    #ifdef OVERLAPPIN_CONST
    // Overlapping Constraints 
    _nConstrs += 2 * nVertices;  
    #endif

    _initCoefs();
    _initVars();
    _initOutputs();
    _updateCoefs();
    _updateOutputs();
}

void Mixedlayout::removeStation(VertexDescriptor vd)
{
    UndirectedGraph        & g            = _metro->g();
    OutEdgeIterator e, e_end;
    vector<EdgeDiscriptor> edgesRemove;
    vector<VertexDescriptor> newEdgesVertex;

    for ( tie( e, e_end ) = out_edges( vd, g ); e != e_end; ++e ) {
        EdgeDescriptor ed = *e;
        VertexDescriptor vS = source( ed, g );
        VertexDescriptor vT = target( ed, g );
        if (vS != vd)
            newEdgesVertex.push_back(vS);
        else if (vT !=  vd)
            newEdgesVertex.push_back(vT);
        edgesRemove.push_back(ed);
    }
    for (auto ed : edgesRemove)
        remove_edge(ed, g);
    remove_vertex(vd, g);


}


/*
* initialize the coefficient
*
*  Inputs
*      none
*
*  Outputs
*      none
*/
void Mixedlayout::_initCoefs( void )
{
    UndirectedGraph        & g            = _metro->g();
    unsigned int nVertices      = _metro->nStations();
    
    // initialization
    unsigned int nRows = 0;
    _coef.resize( _nConstrs, _nVars );
    _coef << Eigen::MatrixXd::Zero( _nConstrs, _nVars );

    // Regular edge
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        unsigned int idS = MIN2( g[ vdS ].id, g[ vdT ].id );
        unsigned int idT = MAX2( g[ vdS ].id, g[ vdT ].id );

        auto w = _w_mixed * g[edge].weight;
        if (g[vdS].smoothPath and g[vdT].smoothPath)
            w = 0;
        // x
        _coef( nRows, idS ) = w;
        _coef( nRows, idT ) = -w;
        nRows++;

        // y
        _coef( nRows, idS + nVertices ) = w;
        _coef( nRows, idT + nVertices ) = -w;
        nRows++;
    }

    // Positional constraints
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        unsigned int id = g[ vertex ].id;
        if (g[vertex].smoothPath) {
            _coef( nRows, id ) = _w_position_path;  // x
            nRows++;
            _coef( nRows, id + nVertices ) = _w_position_path;   // y
            nRows++;
        } else {
            _coef( nRows, id ) = _w_position;
            nRows++;
            _coef( nRows, id + nVertices ) = _w_position;
            nRows++;
        }
    }

    #ifdef OVERLAPPIN_CONST
    // Node overlapping Constratins 
    BGL_FORALL_VERTICES( vertexI, g, UndirectedGraph ) {
        
        unsigned int id = g[ vertexI ].id;
        Coord2 vi = *g[ vertexI ].coordPtr; 
        double w = 0; 
        
        Coord2 vj = _metro->closestPointOnGraph( vertexI );
        Coord2 d = vi - vj; 
        if ( magnitude(d) < _gama ) {
            w = _w_overlap; 
            g[ vertexI ].validPos = false; 
        } else g[ vertexI ].validPos = true;           
                 
        _coef( nRows, id ) = w;
        nRows++; 
        _coef( nRows, id + nVertices ) = w;
        nRows++;  

    }
    #endif //Overlapping
}


/*
* itialize the variables
*
*  Inputs
*      none
*
*  Outputs
*      none
*/
void Mixedlayout::_initVars( void )
{
    UndirectedGraph        & g            = _metro->g();
    unsigned int nVertices      = _metro->nStations();

    // initialization
    _var.resize( _nVars );

    unsigned int nRows = 0;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){

        _var( nRows, 0 ) = g[ vertex ].coordPtr->x();
        _var( nRows + nVertices, 0 ) = g[ vertex ].coordPtr->y();
        nRows++;
    }
}

/*
* updates curent Angle of Edges
* 
* Inputs
*   none
* Outputs
*   none
*/
void Mixedlayout::_updateEdgeCurAngle( void )
{
    UndirectedGraph & g = _metro->g();

    // initialization
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ){

        VertexDescriptor vS = source( edge, g );
        VertexDescriptor vT = target( edge, g );
        Coord2 vi, vj;
        if( g[ vS ].id < g[ vT ].id ){
            vi = *g[ vS ].coordPtr;
            vj = *g[ vT ].coordPtr;
        }
        else{
            vi = *g[ vT ].coordPtr;
            vj = *g[ vS ].coordPtr;
        }
        Coord2 vji = vj - vi;

        double angle = atan2( vji.y(), vji.x() );
        g[ edge ].curAngle = angle;
    }
}
void Mixedlayout::_setTargetAngle2() {
    UndirectedGraph & g = _metro->g();

    double sector[ 8 ] = { -M_PI, -3.0*M_PI/4.0, -M_PI/2.0, -M_PI/4.0, 0.0,
                           M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0};
    vector< vector< VertexDescriptor > > vdVec( 8 );
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        DegreeSizeType degrees = out_degree( vertex, g );
        if (degrees < 8 and degrees != 0)
            vdVec[ 8 - degrees ].push_back( vertex );
    }

    // initialization
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ){
        g[ edge ].target = 2.0*M_PI;
        g[ edge ].isCurved = false;
        g[ edge ].isSmoothPath = false;
        g[ edge ].targetSet = false;
    }

    vector<EdgeDescriptor> edgesToSubdivide;

    // set target angles
    for( int i = 0; i < (int) vdVec.size(); i++ ){
        for( int j = 0; j < (int) vdVec[ i ].size(); j++ ){

            map< double, EdgeDescriptor > circM;
            // sort the angle
            OutEdgeIterator e, e_end;
            for ( tie( e, e_end ) = out_edges( vdVec[i][j], g ); e != e_end; ++e ) {
                EdgeDescriptor ed = *e;
                VertexDescriptor vS = source( ed, g );
                VertexDescriptor vT = target( ed, g );
                double angle = g[ ed ].curAngle;

                if ( g[ vS ].id > g[ vT ].id ) {
                    if ( angle > 0 ) angle = -M_PI + g[ ed ].curAngle;
                    else angle = M_PI + g[ ed ].curAngle;
                }
                circM.insert( pair< double, EdgeDescriptor >( angle, ed ) );
            }

            vector<EdgeDescriptor> circMVector;

            for (auto it : circM) {
                EdgeDescriptor ed = it.second;
                circMVector.push_back(ed);
            }

            // vector<EdgeDescriptor> assignments[8];
            vector< vector<EdgeDescriptor>> assignments = { {}, {}, {}, {}, {}, {}, {}, {}};

            // Add each edge to closes sector
            for (auto ed : circMVector) {
                int closestSector = 0;
                double minDist = 10e6;
                for (int k = 0; k < 8; k++) {
                    double dist = getOffsetToSector(ed, k);
                    if (dist <= minDist) {
                        minDist = dist;
                        closestSector = k;
                    }
                }
                if (fabs(minDist) >= M_PI * 0.2)
                    cout << "min Dist: " << minDist << endl;

                assignments[closestSector].push_back(ed);
            }

            int iterations = 0;
            vector<int> hasBeenMoved;

            while(continioueAssigning(assignments) and iterations < 100) {
                cout << "Assigment is not done " << endl;
                iterations++;
                for (int k = 0; k < assignments.size(); k++) { // check that can not move back
                    if (assignments[k].size() > 1 ) {
                        vector<double> minOffset;
                        vector<int> shift;
                        // Find for each edge the minimal left or right turn
                        for (auto ed : assignments[k]) {
                            if (g[ed].targetSet) {
                                minOffset.push_back(10e6);
                                shift.push_back(0);
                            } else {
                                VertexDescriptor sVD = source(ed, g);
                                VertexDescriptor tVD = target(ed, g);

                                int minSector = (k - 1) % 8;
                                if (minSector < 0 ) minSector = 7;
                                int maxSector = (k + 1) % 8;
                                if (maxSector > 7) maxSector = 0;


                                double deltaMinus = fabs(getOffsetToSector(ed, minSector));
                                if (containsASetEdge(assignments[minSector]) or
                                count(hasBeenMoved.begin(), hasBeenMoved.end(), g[ed].id) or
                                (g[sVD].smoothPath and g[tVD].smoothPath) )// or (assignments[(k - 1) % 8].size() > 1 and iterations > 2))
                                    deltaMinus = 10e6;
                                double deltaPlus = fabs(getOffsetToSector(ed, (k + 1) % 8));
                                if (containsASetEdge(assignments[(k + 1) % 8]) or
                                count(hasBeenMoved.begin(), hasBeenMoved.end(), g[ed].id) or
                                (g[sVD].smoothPath and g[tVD].smoothPath) ) // or ( assignments[(k + 1) % 8].size() > 1 and iterations > 2))
                                    deltaPlus = 10e6;
                                if (deltaPlus < deltaMinus) {
                                    minOffset.push_back(deltaPlus);
                                    shift.push_back(1);
                                } else {
                                    minOffset.push_back(deltaMinus);
                                    shift.push_back(-1);
                                }

                            }
                        }

                        // Find the edge with min Cost to rotate
                        int m = 0;
                        int bestIndex = 0;
                        double minCost = 10e10;
                        for (auto offset : minOffset) {
                            if (offset <= minCost and shift[m] != 0) {
                                minCost = offset;
                                bestIndex = m;
                            }
                            m++;
                        }

                        // move best Edge left or right
                        if (minCost == 10e10) {
                            // cout << "    conflict but all set with id: " << g[vdVec[i][j]].id << endl;
                            // try to shift
                            // or split it and run it again
                            for (auto subEdge : assignments[k]) {
                                if (!count(edgesToSubdivide.begin(), edgesToSubdivide.end(), subEdge))
                                    edgesToSubdivide.push_back(subEdge);
                            }
                        } else {
                            EdgeDescriptor edge = assignments[k][bestIndex];
                            int id = g[edge].id;
                            hasBeenMoved.push_back(id);
                            int s = shift[bestIndex];
                            assignments[k].erase(assignments[k].begin() + bestIndex);
                            int newIndex = (k + s);
                            if (newIndex < 0) newIndex = 7;
                            if (newIndex > 7 ) newIndex = 0;
                            assignments[newIndex].push_back(edge);
                            // cout << "assigned it from sector " << k << " to sector " << newIndex << endl;
                        }
                    }

                }
            }

            if (continioueAssigning(assignments)) {
                cout << "could not find a valid assignment solution" << endl;
            }



            // Apply the assignment
            for (int m = 0; m < 8; m++) {
                auto vec = assignments[m];
                for (auto ed : vec) {

                    VertexDescriptor vS = source( ed, g );
                    VertexDescriptor vT = target( ed, g );

                    Coord2 cordVS =  *g[ vS ].smoothPtr;
                    Coord2 cordVT = *g[ vT ].smoothPtr;

                    if (!g[ed].targetSet and g[vS].smoothPath and g[vT].smoothPath) {
                        g[ ed ].isCurved = true;
                        g[ ed ].isSmoothPath = true;

                        EdgeDescriptor closestGuidEdgeSource = _guide->closestEdge( cordVS  );
                        EdgeDescriptor closestGuidEdgeTarget = _guide->closestEdge( cordVT );

                        double dSource = _guide->distanceToClosestPointOnEdge(closestGuidEdgeSource, cordVS);
                        double dTarget = _guide->distanceToClosestPointOnEdge(closestGuidEdgeTarget, cordVT);

                        Coord2 guidSCord = _guide->closestPointOnEdge(closestGuidEdgeSource, cordVS);
                        Coord2 guidTCord = _guide->closestPointOnEdge(closestGuidEdgeTarget, cordVT);
                        Coord2 gVec = guidTCord - guidSCord;
                        double guidAngle = atan2(gVec.y(), gVec.x());

                        double curAngle = g[ed].curAngle;

                        if (abs( guidAngle - curAngle) > M_PI_2 ) {
                            if (abs(guidAngle) > M_PI_2 && abs(curAngle) < M_PI_2 ) {
                                if( guidAngle > 0 ) guidAngle -= M_PI;
                                else guidAngle += M_PI;
                            }
                            else if (abs(guidAngle) < M_PI_2 && abs(curAngle) > M_PI_2) {
                                if( guidAngle > 0 ) guidAngle -= M_PI;
                                else guidAngle += M_PI;

                            }
                            if( abs(guidAngle - curAngle) > M_PI_4 ) {
                                guidAngle = curAngle;
                                g[ ed ].isCurved = false;
                            }
                        }
                    }
                    else if (!g[ed].targetSet) {
                        double targetAngle = sector[m];
                        if ( g[ vS ].id > g[ vT ].id ) {
                            if ( targetAngle > 0 ) {
                                g[ ed ].target = -M_PI + targetAngle;
                            }
                            else {
                                g[ ed ].target = M_PI + targetAngle;
                            }
                        }
                        else {
                            g[ ed ].target = targetAngle;
                        }
                        g[ed].targetSet = true;
                    }
                }
            }
        }
    }
}



bool Mixedlayout::containsASetEdge(vector<EdgeDescriptor> edges ) {
    UndirectedGraph & g = _metro->g();
    bool result = false;
    for (auto ed : edges ) {
        if (g[ed].targetSet)
            result = true;
    }
    return result;
}

bool Mixedlayout::continioueAssigning(vector<vector<EdgeDescriptor>> assignments) {
    bool cont = false;
    for (int i = 0; i < 8; i++) {
        if (assignments[i].size() > 1) {
            cont = true;
        }
    }
    return cont;
}

double Mixedlayout::getOffsetToSector( EdgeDescriptor ed, int sector) {
    UndirectedGraph & g = _metro->g();
    double sectors[ 8 ] = { -M_PI, -3.0*M_PI/4.0, -M_PI/2.0, -M_PI/4.0, 0.0,
                           M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0};

    VertexDescriptor vS = source( ed, g );
    Coord2 cordVS =  *g[ vS ].smoothPtr;
    VertexDescriptor vT = target( ed, g );
    Coord2 cordVT = *g[ vT ].smoothPtr;

    double angle = g[ ed ].curAngle;
    if ( g[ vS ].id > g[ vT ].id ) {
        if ( angle > 0 ) angle = -M_PI + g[ ed ].curAngle;
        else angle = M_PI + g[ ed ].curAngle;
    }

    double dist = fabs( sectors[ sector ] - angle );
    if (sector == 0) {
        // double dist2 = fabs( M_PI - angle );
        // if (dist2 < dist)
        //     dist = dist2;
        dist = min(dist, fabs( M_PI - angle ));
    }
    return dist;
}

/*
* Set the target angle of an edge, octilienar or smooth edge
* 
* Inputs
*   none
* Outputs
*   none
*/
void Mixedlayout::_setTargetAngle( void )   
{
    UndirectedGraph & g = _metro->g();

    double sector[ 9 ] = { -M_PI, -3.0*M_PI/4.0, -M_PI/2.0, -M_PI/4.0, 0.0,
                           M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI };
    vector< vector< VertexDescriptor > > vdVec( OCTILINEAR_SECTOR );
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        DegreeSizeType degrees = out_degree( vertex, g );
        if (degrees < 9 and degrees != 0) {
            vdVec[ OCTILINEAR_SECTOR - degrees ].push_back( vertex );
        } else {
            cerr << "this is a problem" << endl;
        }
    }

    // initialization
    BGL_FORALL_EDGES( edge, g, UndirectedGraph ){
        g[ edge ].target = 2.0*M_PI;
        g[ edge ].isCurved = false;
        g[ edge ].isSmoothPath = false;
        g[ edge ].targetSet = false;
    }

    // set target angles
    for( int i = 0; i < (int) vdVec.size(); i++ ){
        for( int j = 0; j < (int) vdVec[ i ].size(); j++ ){

            map< double, EdgeDescriptor > circM;
            // sort the angle
            OutEdgeIterator e, e_end;
            for ( tie( e, e_end ) = out_edges( vdVec[i][j], g ); e != e_end; ++e ) {
                EdgeDescriptor ed = *e;
                VertexDescriptor vS = source( ed, g );
                VertexDescriptor vT = target( ed, g );
                double angle = g[ ed ].curAngle;

                if ( g[ vS ].id > g[ vT ].id ) {
                    if ( angle > 0 ) angle = -M_PI + g[ ed ].curAngle;
                    else angle = M_PI + g[ ed ].curAngle;
                }

                circM.insert( pair< double, EdgeDescriptor >( angle, ed ) );
            }

            vector<EdgeDescriptor> circMVector;

            // for ( map< double, EdgeDescriptor >::iterator it = circM.begin();
            // it != circM.end(); it++ ) {
            for (auto it : circM) {
                EdgeDescriptor ed = it.second;
                circMVector.push_back(ed);
            }
            // assign the sector
            int index = 0;
            int firstIndex = 0;
            bool firstIndexSet = false;
            bool validAssignment = false;
            int rotations = 0;
            // for ( map< double, EdgeDescriptor >::iterator it = circM.begin();
            //      it != circM.end(); it++ ) {
            //    EdgeDescriptor ed = it->second;

            vector<double> assignments;
            double minCosts = 10e10;

            // Try all rotations
            while (!validAssignment and rotations < circMVector.size()) {
                validAssignment = true;
                for (EdgeDescriptor ed : circMVector) {
                    VertexDescriptor vS = source( ed, g );
                    Coord2 cordVS =  *g[ vS ].smoothPtr;
                    VertexDescriptor vT = target( ed, g );
                    Coord2 cordVT = *g[ vT ].smoothPtr;

                    double curAngle = g[ed].curAngle;

                    if (!g[ed].targetSet) {
                        double angle = g[ ed ].curAngle;
                        if ( g[ vS ].id > g[ vT ].id ) {
                            if ( angle > 0 ) angle = -M_PI + g[ ed ].curAngle;
                            else angle = M_PI + g[ ed ].curAngle;
                        }
                        double target = 0, minDist = M_PI * 4.0;


                        /**
                         * index should always be below 9 as long as degree of the vertex is larger than 8
                         * which never happens in one of the test example
                         */
                        if (index < 9 + firstIndex) {
                            /**
                             * loop will not run when index is larger than 8.
                             * Than the target angle of the edge will be the default as defined above. –> double target = 0
                             * So independet of the angle of the edge.
                             *
                             * Maybe such an scenario is the problem:
                             *      the first angle is close (lets say index 0)
                             *      but we try all other sections as well.
                             *      And because of some error, the last possible segment (index 8)
                             *      is just a bit closer.
                             *      And therefore all the next edges whit a greater angle
                             *      have no "free segement" left
                             */
                            for ( int k = index; k < 9 + firstIndex; k++ ) {
                                double dist = fabs( sector[ k % 9 ] - angle );
                                if( minDist > dist) {
                                    minDist = dist;
                                    target = sector[ k % 9 ];
                                    index = k + 1;
                                }
                            }
                        } else {
                            validAssignment = false;
                            /**
                             * this section of the code should never be executed in our case.
                             * Because we sorted the edges first.
                             *
                             * Quick fix I have added: in case the index > 9 –> so according to the programm
                             * there is no "free angle availabe anymore" we just take the
                             * direction of the 9 possible that is the closest to the current angle.
                             * Has the problem that we can assign the same direction to incident edges.
                             */
                            cout << "!!!!!!!!!! index >= 9 " << angle << ", " << index << " first index: " << firstIndex << endl;
                            double minDelta = 10e10;
                            for (auto sec : sector) {
                                double delta = abs(angle - sec);
                                if (delta < minDelta) {
                                    minDelta = delta;
                                    target = sec;
                                }
                            }
                            cout << "chosen index: " << target << endl;
                        }


                        if ( g[ vS ].id > g[ vT ].id ) {
                            if ( target > 0 ) {
                                // g[ ed ].target = -M_PI + target;
                                assignments.push_back(-M_PI + target);
                            }
                            else {
                                // g[ ed ].target = M_PI + target;
                                assignments.push_back(M_PI + target);
                            }
                        }
                        else {
                            // g[ ed ].target = target;
                            assignments.push_back(target);
                        }
                    }
                    else {
                        double target = g[ ed ].target;
                        assignments.push_back(target);

                        if ( g[ vS ].id > g[ vT ].id ) {
                            if ( target > 0 ) target = -M_PI + g[ ed ].target;
                            else target = M_PI + g[ ed ].target;
                        }
                        for( int k = index; k < 9; k++ ){
                            if( target == sector[ k ] ) {
                                index = k+1;
                            }
                        }
                    }

                    if (!firstIndexSet) {
                        firstIndex = index - 1;
                        firstIndexSet = true;
                        if (firstIndex > 8)
                            cout << "first Index is > 8 " << firstIndex << endl; // DEBUG
                    }

                    if (g[vS].smoothPath and g[vT].smoothPath) {
                        g[ ed ].isCurved = true;
                        g[ ed ].isSmoothPath = true;

                        EdgeDescriptor closestGuidEdgeSource = _guide->closestEdge( cordVS  );
                        EdgeDescriptor closestGuidEdgeTarget = _guide->closestEdge( cordVT );

                        double dSource = _guide->distanceToClosestPointOnEdge(closestGuidEdgeSource, cordVS);
                        double dTarget = _guide->distanceToClosestPointOnEdge(closestGuidEdgeTarget, cordVT);

                        Coord2 guidSCord = _guide->closestPointOnEdge(closestGuidEdgeSource, cordVS);
                        Coord2 guidTCord = _guide->closestPointOnEdge(closestGuidEdgeTarget, cordVT);
                        Coord2 gVec = guidTCord - guidSCord;
                        double guidAngle = atan2(gVec.y(), gVec.x());

                        // double guidAngle = atan2(gVec.y(), gVec.x());
                        double curAngle = g[ed].curAngle;

                        if (abs( guidAngle - curAngle) > M_PI_2 ) {
                            if (abs(guidAngle) > M_PI_2 && abs(curAngle) < M_PI_2 ) {
                                if( guidAngle > 0 ) guidAngle -= M_PI;
                                else guidAngle += M_PI;
                            }
                            else if (abs(guidAngle) < M_PI_2 && abs(curAngle) > M_PI_2) {
                                if( guidAngle > 0 ) guidAngle -= M_PI;
                                else guidAngle += M_PI;

                            }
                            if( abs(guidAngle - curAngle) > M_PI_4 ) {
                                guidAngle = curAngle;
                                g[ ed ].isCurved = false;
                            }
                        }
                        // g[ ed ].target = guidAngle;     // set the target angle to guid angle
                        assignments.back() = guidAngle;

                    }
                }

                // Todo calc cost of assignment
                // rotate vector to try another sequence
                int counter_i = 0;
                double sumCosts = 0;
                for (auto tar : assignments) {
                    auto ed = circMVector[counter_i];
                    double curAngle = g[ed].curAngle;
                    double cost = fabs(tar - curAngle);
                    sumCosts += cost;
                    counter_i++;
                }
                // High Cost for not valid assignemnt
                if (!validAssignment) {
                    sumCosts += 1000.0;
                }
                // Apply the target if better than the previouse asignment
                if (sumCosts < minCosts ) {
                    minCosts = sumCosts;
                    cout << "min Costs " << minCosts << endl;
                    counter_i = 0;
                    for (auto tar : assignments) {
                        auto ed = circMVector[counter_i];
                        auto sVD = source(ed, g);
                        auto tVD = target(ed, g);
                        if (!g[sVD].smoothPath and !g[tVD].smoothPath) {
                            g[ed].target = tar;
                        }
                        counter_i++;
                    }
                }
                assignments.clear();

                rotations += 1;
                if (!validAssignment or true) {
                    index = 0;
                    firstIndex = 0;
                    firstIndexSet = false;

                    rotate(circMVector.begin(), circMVector.begin()+1, circMVector.end());
                    if (rotations >= circMVector.size() - 1 and minCosts > 100)
                        cout << "could not find a valid rotation: " << rotations << " id: " << g[vdVec[i][j]].id << endl;

                }
            }

            // set target as set
            for (auto ed : circMVector) {
                g[ed].targetSet = true;
            }
        }
    }
}



/*
* initialize the output
*
*  Inputs
*      none
*  Outputs
*      none
*/
void Mixedlayout::_initOutputs( void )
{
    UndirectedGraph        & g            = _metro->g();

    _updateEdgeCurAngle();

    // initialization
    unsigned int nRows = 0;
    _output.resize( _nConstrs );
    _output << Eigen::VectorXd::Zero( _nConstrs );


    // Regular edge octilinearty
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi, vj;
        if( g[ vdS ].id < g[ vdT ].id ){
            vi = *g[ vdS ].smoothPtr;
            vj = *g[ vdT ].smoothPtr;
        }
        else{
            vi = *g[ vdT ].smoothPtr;
            vj = *g[ vdS ].smoothPtr;
        }
        Coord2 vji = vi - vj;

        double angle = g[ edge ].curAngle;
        double theta = g[ edge ].target - angle;

        double cosTheta = cos( theta ), sinTheta = sin( theta );
       // double beta = _targetLength;
       // if (!g[ vdS ].isStation or !g[vdT].isStation) {
       //     beta *= 0.5; // TODO
       // }
       // beta *= ( g[edge].stationsOnEdge + 1.0 );

       // double s = beta * g[ edge ].weight / vji.norm();
        double s = 1.0;
        auto w = _w_mixed * g[edge].weight;
        if (g[vdS].smoothPath and g[vdT].smoothPath)
            w = 0;
        // x
        _output( nRows, 0 ) = w * s * ( cosTheta * vji.x() - sinTheta * vji.y() );
        nRows++;

        // y
        _output( nRows, 0 ) = w * s * ( sinTheta * vji.x() + cosTheta * vji.y() );
        nRows++;
    }


    // Positional constraints
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        double w = _w_position;
        Coord2 pos = *g[ vertex ].smoothPtr;
        if (g[vertex].smoothPath) {
            w = _w_position_path;
            pos = g[vertex].closestPointOnEdge;
        }

        _output( nRows, 0 ) = w * pos.x();
        nRows++;
        _output( nRows, 0 ) = w * pos.y();
        nRows++;
    }


    #ifdef OVERLAPPIN_CONST
    // Node overlapping Constratins 
    BGL_FORALL_VERTICES( vertexI, g, UndirectedGraph ) {
        Coord2 vi = *g[ vertexI ].smoothPtr; 
        double w = 0; 
        Coord2 direction = Coord2(0,0); 

        Coord2 vj = _metro->closestPointOnGraph( vertexI );
        Coord2 d = vi - vj; 
        if ( magnitude(d) < _gama ) {
            w = _w_overlap; 
            direction = d; 
        }       
        
        double magDir = magnitude(direction); 
        if (magDir > 0) {
            double r = _gama * 0.5 / magDir; 
            direction *= r; 
        }

        _output(nRows,0) = w * ( g[ vertexI ].smoothPtr->x() + direction.x()); 
        nRows++; 
        _output( nRows, 0 ) = w * (g[ vertexI ].smoothPtr->y() + direction.y());   
        nRows++; 
    }
    #endif //Overlapping
}


/*
* update the coefs
*
*  Inputs
*      none
*
*  Outputs
*      none
*/
void Mixedlayout::_updateCoefs( void )
{
    UndirectedGraph               & g             = _metro->g();
    unsigned int        nVertices       = _metro->nStations();
    unsigned int        nVE             = 0;
    unsigned int        nB              = 0;
    vector< double >    ratioR          = _metro->ratioR();

    Eigen::MatrixXd     oldCoef;
    oldCoef = _coef.block( 0, 0, _nConstrs, _nVars );
    unsigned int nRows = 0; 

    #ifdef OVERLAPPIN_CONST
    nRows =  _nConstrs - _nVars; 
    BGL_FORALL_VERTICES( vertexI, g, UndirectedGraph ) {
    
        unsigned int id = g[ vertexI ].id;
        Coord2 vi = *g[ vertexI ].coordPtr; 
        double w = 0;       
        
        Coord2 vj = _metro->closestPointOnGraph( vertexI );
        Coord2 d = vi - vj; 
        if ( magnitude(d) < _gama) {
            w = _w_overlap; 
            g[ vertexI ].validPos = false; 
        } else g[ vertexI ].validPos = true;

        oldCoef( nRows, id ) = w;
        nRows++; 
        oldCoef( nRows, id + nVertices) = w;
        nRows++;  

    }
    #endif // OVERLAPPIN_CONST

    #ifdef OCTILINEAR_CONFLICT
        nVE = _metro->VEconflict().size();
    #endif // OCTILINEAR_CONFLICT
    #ifdef OCTILINEAR_BOUNDARY
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph )
    {
        double minD = _d_Beta/2.0;
        if ( g[ vd ].coordPtr->x() <= -( _half_width - minD ) ) nB++;
        if ( g[ vd ].coordPtr->x() >= ( _half_width - minD ) ) nB++;
        if ( g[ vd ].coordPtr->y() <= -( _half_height - minD ) ) nB++;
        if ( g[ vd ].coordPtr->y() >= ( _half_height - minD ) ) nB++;
    }
    #endif // OCTILINEAR_BOUNDARY
    _coef.resize( _nConstrs + nB + 2*nVE, _nVars );
    // copy old coefficient
    _coef << oldCoef, Eigen::MatrixXd::Zero( nB + 2*nVE, _nVars );
    nRows = _nConstrs;

    #ifdef  OCTILINEAR_BOUNDARY
    // add boundary coefficient
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph ){

        unsigned int id = g[ vd ].id;
        double minD = _d_Beta/2.0;

        if ( g[ vd ].coordPtr->x() <= -( _half_width - minD ) ) {
            _coef( nRows, id ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].coordPtr->x() >= ( _half_width - minD ) ) {
            _coef( nRows, id ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].coordPtr->y() <= -( _half_height - minD ) ) {
            _coef( nRows, id + nVertices ) = _w_boundary;
            nRows++;
        }
        if ( g[ vd ].coordPtr->y() >= ( _half_height - minD ) ) {
            _coef( nRows, id + nVertices ) = _w_boundary;
            nRows++;
        }
    }
    #endif  // OCTILINEAR_BOUNDARY
    #ifdef  OCTILINEAR_CONFLICT
    // copy conflict coefficient
    unsigned int countVE = 0;
    for (auto it : _metro->VEconflict()) {
            VertexDescriptor vdV = it.second.first;
            EdgeDescriptor ed = it.second.second;
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
    #endif  // OCTILINEAR_CONFLICT
}

/* 
* update the output
* 
*  Inputs
*      none
*
*  Outputs
*      none
*/
void Mixedlayout::_updateOutputs( void )
{
    UndirectedGraph     & g             = _metro->g();
    unsigned int        nVE             = 0;
    unsigned int        nB              = 0;
    vector< double >    ratioR          = _metro->ratioR();
    vector< double >    ratioGeo        = _metro->ratioGeo();

    unsigned int nRows = 0;
    Eigen::VectorXd     oldOutput;
    oldOutput = _output;
    #ifdef  OCTILINEAR_CONFLICT
    nVE = _metro->VEconflict().size();
    #endif  // OCTILINEAR_CONFLICT
    #ifdef  OCTILINEAR_BOUNDARY
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph )
    {
        double minD = _d_Beta/2.0;
        if ( g[ vd ].coordPtr->x() <= -( _half_width - minD ) ) nB++;
        if ( g[ vd ].coordPtr->x() >= ( _half_width - minD ) ) nB++;
        if ( g[ vd ].coordPtr->y() <= -( _half_height - minD ) ) nB++;
        if ( g[ vd ].coordPtr->y() >= ( _half_height - minD ) ) nB++;
    }
    #endif  // OCTILINEAR_BOUNDARY
    _output.resize( _nConstrs + nB + 2*nVE );
    _output << Eigen::VectorXd::Zero( _nConstrs + nB + 2*nVE );

    _updateEdgeCurAngle();

    // Regular edge length
    BGL_FORALL_EDGES( edge, g, UndirectedGraph )
    {
        VertexDescriptor vdS = source( edge, g );
        VertexDescriptor vdT = target( edge, g );
        Coord2 vi, vj;
        if( g[ vdS ].id < g[ vdT ].id ){
            vi = *g[ vdS ].coordPtr;
            vj = *g[ vdT ].coordPtr;
        }
        else{
            vi = *g[ vdT ].coordPtr;
            vj = *g[ vdS ].coordPtr;
        }
        Coord2 vji = vi - vj;
        double angle = g[ edge ].curAngle;
        double theta = g[ edge ].target - angle;
        double cosTheta = cos( theta ), sinTheta = sin( theta );
        double s = 1.0;
        auto w = _w_mixed * g[edge].weight;
        if (g[vdS].smoothPath and g[vdT].smoothPath)
            w = 0;
        // x
        _output( nRows, 0 ) = w * s * ( cosTheta * vji.x() - sinTheta * vji.y() );
        nRows++;
        // y
        _output( nRows, 0 ) = w * s * ( sinTheta * vji.x() + cosTheta * vji.y() );
        nRows++;
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

    #ifdef OVERLAPPIN_CONST
    // Node overlapping Constratins 
    BGL_FORALL_VERTICES( vertexI, g, UndirectedGraph ) {
        Coord2 vi = *g[ vertexI ].coordPtr; 
        double w = 0; 
        Coord2 direction = Coord2(0,0); 

        Coord2 di = _metro->closestPointOnGraph( vertexI );
        Coord2 d = vi - di; 
       
        if ( magnitude(d) < _gama) {
            w = _w_overlap; 
            direction = d; 
            double magDir = direction.norm();
            if (direction.norm() < 0.1 ) {
                // cout << "could be a problem " << d << endl;
            } else {
                    double r = _gama / magDir; 
                    direction *= r; 
            }
        }          
        _output(nRows,0) = w * ( g[ vertexI ].coordPtr->x() + direction.x() ); 
        nRows++; 
        _output( nRows, 0 ) = w * ( g[ vertexI ].coordPtr->y() + direction.y() );   
        nRows++; 
    }
    #endif //Overlapping

    #ifdef  OCTILINEAR_BOUNDARY
    // boundary constraints
    BGL_FORALL_VERTICES( vd, g, UndirectedGraph ){

        double minD = _d_Beta/2.0;

        if ( g[ vd ].coordPtr->x() <= -( _half_width - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * -( _half_width - minD );
            nRows++;
        }
        if ( g[ vd ].coordPtr->x() >= ( _half_width - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * ( _half_width - minD );
            nRows++;
        }
        if ( g[ vd ].coordPtr->y() <= -( _half_height - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * -( _half_height - minD );
            nRows++;
        }
        if ( g[ vd ].coordPtr->y() >= ( _half_height - minD ) ) {
            _output( nRows, 0 ) = _w_boundary * ( _half_height - minD );
            nRows++;
        }
    }
#endif  // OCTILINEAR_BOUNDARY

#ifdef  OCTILINEAR_CONFLICT
    // copy conflict coefficient
    unsigned int countVE = 0;
    // for ( VEMap::iterator it = _metro->VEconflict().begin();
    //       it != _metro->VEconflict().end(); ++it ) {
    for (auto it : _metro->VEconflict()) {
        VertexDescriptor vdV = it.second.first;
        g[vdV].collision = true;
        EdgeDescriptor ed = it.second.second;
        VertexDescriptor vdS = source( ed, g );
        VertexDescriptor vdT = target( ed, g );
        double r = ratioR[ countVE ];
        // r = ratioGeo[ countVE ];

        Coord2 v = *g[ vdV ].coordPtr; // why is geoPtr here
        Coord2 p = r * *g[ vdS ].coordPtr + ( 1.0-r ) * *g[ vdT ].coordPtr; // why is geo ptr here
        double minD = _d_Beta/2.0;
        if (!g[vdS].isStation or !g[vdT].isStation)
            minD *= 0.5;
        double delta;
        if ( ( v - p ).norm() > 0.00001 ) {
            delta = minD / ( v - p ).norm();
        } else {
            // cerr << "lies directly on edge " << endl;
            delta = minD;
            Coord2 v1 = *g[vdS].coordPtr - *g[vdT].coordPtr;
            Coord2 normal = Coord2(0,0);
            double angle = M_PI/2;
            normal.setX(v1.x() * cos(angle) - v1.y() * sin(angle));
            normal.setY(v1.x() * sin(angle) + v1.y() * cos(angle));
            normal = normal.normalize();
            // cerr << "normal: " << normal << endl;
            v = v + normal;
            // move in direction of edge
            // v = v + (*g[vdV].smoothPtr - v).normalize(); // looks like this is the issue – normalize when both equal is devide by zero

        }


        // x
        _output( nRows, 0 ) = _w_crossing * delta * ( v - p ).x();
        nRows++;
        // y
        _output( nRows, 0 ) = _w_crossing * delta * ( v - p ).y();
        nRows++;

        countVE++;
    }
#endif  // OCTILINEAR_CONFLICT
}

/*
* optimization
*
*  Inputs
*      none
*  Outputs
*      err:
*/
double Mixedlayout::LeastSquare( unsigned int iter )
{
    double mse = 0.0;
    for( int i = 0; i < (int) iter; i++ ) {

        Eigen::VectorXd last_var = _var;
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

        if( (i+1) % 100 == 0 )
            cerr << setprecision( 10 ) << "Loop(" << i << ") mse = " << mse << endl;
        if( ( mse ) < 5.0e-4 ) break;
    }
    return mse;
}

/*
* ConjugateGradient optimization
*
*  Inputs
*      none
*  Outputs
*      err
*/
double Mixedlayout::ConjugateGradient( unsigned int iter )
{
    // initialization, prepare the square matrix
    // initialization, prepare the square matrix
    Eigen::MatrixXd A;
    Eigen::VectorXd b, Ap;
    A = _coef.transpose() * _coef;
    b = _coef.transpose() * _output;

    // initialization
    Eigen::VectorXd err = b - A * _var;
    Eigen::VectorXd p = err;
    double rsold = err.adjoint() * err;

    // main algorithm
    for( int i = 0; i < (int) (iter * 1.); i++ ) {

        // prepare the square matrix
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
            cout << "Mixed ... " << i << " / " << iter << " ... " << endl; 
        }

        // update
        retrieve();
        _updateCoefs();
        _updateOutputs();
    }

    return sqrt( err.adjoint() * err );
}

/**
 * Postprocessing of the mixed layout of Guide Shapes that have
 * vertex with degree > 2. Or more than one poly-line.
 */
void Mixedlayout::complexPostProcessing() {
    UndirectedGraph &g = _metro->g();
    UndirectedGraph &guideG = _guide->g();
    // finding out all in between vertex
    map<EdgeDescriptor, vector<VertexDescriptor>> betweenGVert;
    BGL_FORALL_EDGES(edge, g, UndirectedGraph ) {
        if (g[edge].isSmoothPath && g[edge].originalMetroEdge) {
            EdgeDescriptor edI = _guide->closestEdge(*g[source(edge, g)].coordPtr);
            EdgeDescriptor edJ = _guide->closestEdge(*g[target(edge, g)].coordPtr);
            VertexDescriptor uiS = source(edI, guideG);
            VertexDescriptor uiT = target(edI, guideG);
            VertexDescriptor ujS = source(edJ, guideG);
            VertexDescriptor ujT = target(edJ, guideG);

            // Find the longest Shortest Path. To decide if Source or Target of edge
            vector<VertexDescriptor> betweenVertex = _guide->shortesPath(uiS, ujS);
            vector<VertexDescriptor> betweenVertexST = _guide->shortesPath(uiS, ujT);
            vector<VertexDescriptor> betweenVertexTS = _guide->shortesPath(uiT, ujS);
            vector<VertexDescriptor> betweenVertexTT = _guide->shortesPath(uiT, ujT);
            // Select the longest of the four
            if (betweenVertexST.size() > betweenVertex.size() )
                betweenVertex = betweenVertexST;
            if (betweenVertexTS.size() > betweenVertex.size() )
                betweenVertex = betweenVertexTS;
            if (betweenVertexTT.size() > betweenVertex.size() )
                betweenVertex = betweenVertexTT;
            // Remove the two vertex that are not in the path between the two point ons guide edge
            if (betweenVertex.size() > 0) betweenVertex.erase(betweenVertex.end() - 1);
            if (betweenVertex.size() > 0) betweenVertex.erase(betweenVertex.begin());
            if (betweenVertex.size() > 0) {
                betweenGVert.insert(pair<EdgeDescriptor, vector<VertexDescriptor>>(edge, betweenVertex));
            }
        }
    }

    // Remove when two metro edges share a guide edge
    BGL_FORALL_EDGES(e_i, g, UndirectedGraph ) {
        if (g[e_i].isSmoothPath) {
            BGL_FORALL_EDGES(e_j, g, UndirectedGraph ) {
                if (g[e_i].isSmoothPath && g[e_j].isSmoothPath) {
                    // check if they share guide vertex in there paths
                    if (find_first_of (betweenGVert[e_i].begin(), betweenGVert[e_i].end(),
                                       betweenGVert[e_j].begin(), betweenGVert[e_j].end()) != betweenGVert[e_i].end()) {
                        if (betweenGVert[e_i].size() > betweenGVert[e_j].size()) {
                            g[e_i].isSmoothPath = false;
                        } else {

                        }
                    }
                }
            }
        }
    }

    // ADDing Additional Vertex
    int addVertex = 0;
    int addEdges = 0;
    for (auto item : betweenGVert) {
        EdgeDescriptor edge = item.first;
        auto uDescriptorInbetween = item.second;
        VertexDescriptor prevVert = NULL;
        VertexDescriptor sVD = source(edge, g);
        VertexDescriptor tVD = target(edge, g);
        VertexDescriptor first = NULL;
        VertexDescriptor last = NULL;

        if (g[edge].isSmoothPath && g[edge].originalMetroEdge && uDescriptorInbetween.size() > 0) {
            g[edge].isVisible = false;
            double count = 0;
            for (VertexDescriptor guideVD : uDescriptorInbetween) {
                // Add intermediat vertices
                double x = guideG[guideVD].coordPtr->x();
                double y = guideG[guideVD].coordPtr->y();
                VertexDescriptor vij = add_vertex(g);
                g[vij].coordPtr = new Coord2(x, y);
                g[vij].smoothPtr = new Coord2(x, y);
                g[vij].octilinearPtr = new Coord2(x, y);
                g[vij].intersection = false;
                g[vij].metroShape = true;
                g[vij].isStation = false;
                g[vij].inflectionPoint = false;
                g[vij].closestPointOnEdge = Coord2(1000, 1000);
                g[vij].id = _metro->getNumStations() + addVertex;
                addVertex++;

                // Find the first Vertex of the edge
                if (prevVert == NULL ) {
                    double distSource = (*guideG[sVD].coordPtr - *guideG[vij].coordPtr).norm();
                    double distTarget = (*guideG[tVD].coordPtr - *guideG[vij].coordPtr).norm();
                    if (distSource < distTarget) {
                        first = sVD;
                        last = tVD;
                    } else {
                        first = tVD;
                        last = sVD;
                    }
                    prevVert = first;
                }

                // Move vertex with offset
                double k = count / (uDescriptorInbetween.size() + 2.0);
                Coord2 distFirst = *g[first].coordPtr - g[first].closestPointOnEdge;
                Coord2 distLast = *g[last].coordPtr - g[last].closestPointOnEdge;
                // cout << "delta " << distFirst << " , " << distLast << endl;
                double offX = (1.0 - k) * distFirst.x() + k * distLast.x();
                double offY = (1.0 - k) * distFirst.y() + k * distLast.y();
                g[vij].coordPtr->setX(x + offX);
                g[vij].coordPtr->setY(y + offY);

                // Add new Edge
                EdgeDescriptor e1; bool b1;
                tie(e1,b1) = add_edge(prevVert, vij, g);
                g[e1].isSmoothPath = true;
                g[e1].isMetro = true;
                g[e1].isCurved = false;
                g[e1].matchPath = false;
                g[e1].lineID = g[edge].lineID;
                g[e1].originalMetroEdge = false;
                addEdges++;
                g[e1].id = _metro->getNumEdges() + addEdges;
                prevVert = vij;

                count++;
            }
            // ADD last intermediate edge
            if (prevVert != last) {
                EdgeDescriptor e1; bool b1;
                tie(e1,b1) = add_edge(prevVert, last, g);
                g[e1].isSmoothPath = true;
                g[e1].isMetro = true;
                g[e1].isCurved = false;
                g[e1].matchPath = false;
                g[e1].lineID = g[edge].lineID;
                g[e1].originalMetroEdge = false;
                addEdges++;
                g[e1].id = _metro->getNumEdges() + addEdges;
            }
        }
    }
}

void Mixedlayout::preparePostProcessing() {
    cout << "preparing post processing" << endl;
    UndirectedGraph &g = _metro->g();
    UndirectedGraph &guideG = _guide->g();
    // finding out all inbetween vertex
    map<EdgeDescriptor, vector<unsigned int>> betweenGVert;
    BGL_FORALL_EDGES(edge, g, UndirectedGraph ) {
        if (g[edge].isSmoothPath && g[edge].originalMetroEdge) {
            EdgeDescriptor edS = _guide->closestEdge(*g[source(edge, g)].coordPtr);
            EdgeDescriptor edT = _guide->closestEdge(*g[target(edge, g)].coordPtr);
            int lowGuideIndex = min(
                    min(guideG[source(edS, guideG)].id, guideG[target(edS, guideG)].id),
                    min(guideG[source(edT, guideG)].id, guideG[target(edT, guideG)].id)
                    );

            int highGuideIndex = max(
                    max(guideG[source(edS, guideG)].id, guideG[target(edS, guideG)].id),
                    max(guideG[source(edT, guideG)].id, guideG[target(edT, guideG)].id)
                    );
            VertexDescriptor gi, gn;
            BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                if (guideG[vertex].id == lowGuideIndex ) gi = vertex;
                else if (guideG[vertex].id == highGuideIndex ) gn = vertex;
            }

            VertexDescriptor vi = source(edge, g);
            VertexDescriptor vj = target(edge, g);

            vector<unsigned int> betweenVertex;

            if (highGuideIndex - lowGuideIndex < _guide->getNumEdges() - highGuideIndex + lowGuideIndex ) { // zerro is not inbetween
                BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                    if (guideG[vertex].id > lowGuideIndex && guideG[vertex].id < highGuideIndex ) {
                        betweenVertex.push_back(guideG[vertex].id);
                    }
                }
            } else {
                BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                    if (guideG[vertex].id < lowGuideIndex || guideG[vertex].id > highGuideIndex ) {
                        betweenVertex.push_back(guideG[vertex].id);
                    }
                }
            }
            betweenGVert.insert(pair<EdgeDescriptor, vector<unsigned int>>(edge, betweenVertex));
        }
    }
    // Remove when share guide Vertex
    BGL_FORALL_EDGES(e_i, g, UndirectedGraph ) {
        if (g[e_i].isSmoothPath) {
            BGL_FORALL_EDGES(e_j, g, UndirectedGraph ) {
                if (g[e_i].id < g[e_j].id && g[e_j].isSmoothPath && g[e_i].isSmoothPath) {
                // if ( g[e_j].isSmoothPath && g[e_i].isSmoothPath && g[e_i].id != g[e_j].id) {
                    auto e_iBetween = betweenGVert[e_i];
                    auto e_jBetween = betweenGVert[e_j];
                    // share element
                    if (find_first_of (e_iBetween.begin(), e_iBetween.end(),
                                       e_jBetween.begin(), e_jBetween.end()) != e_iBetween.end()) {
                        if (e_iBetween.size() > e_jBetween.size() ) {
                            g[e_i].isSmoothPath = false;
                        } else {
                            cout << "size: " << e_iBetween.size() << ", " << e_jBetween.size() << endl;
                            g[e_j].isSmoothPath = false;
                        }
                        cout << "removed edge with id: " << g[e_j].id << endl;
                    } else {

                    }
                }
            }
        }
    }
    cout << "ended preparing post processing" << endl;
}

void Mixedlayout::postProcessing() {
    cerr << "start post processing" << endl;
    UndirectedGraph &g = _metro->g();
    UndirectedGraph &guideG = _guide->g();

    int addVertex = 0;
    int addEdges = 0;
    BGL_FORALL_EDGES(edge, g, UndirectedGraph ) {
        if (g[edge].isSmoothPath && g[edge].originalMetroEdge) {
            EdgeDescriptor edS = _guide->closestEdge(*g[source(edge, g)].coordPtr);
            EdgeDescriptor edT = _guide->closestEdge(*g[target(edge, g)].coordPtr);

            int lowGuideIndex = min(
                    min(guideG[source(edS, guideG)].id, guideG[target(edS, guideG)].id),
                    min(guideG[source(edT, guideG)].id, guideG[target(edT, guideG)].id)
                    );

            int highGuideIndex = max(
                    max(guideG[source(edS, guideG)].id, guideG[target(edS, guideG)].id),
                    max(guideG[source(edT, guideG)].id, guideG[target(edT, guideG)].id)
            );

            VertexDescriptor gi, gn;
            BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                if (guideG[vertex].id == lowGuideIndex ) gi = vertex;
                else if (guideG[vertex].id == highGuideIndex ) gn = vertex;
            }

            VertexDescriptor vi = source(edge, g);
            VertexDescriptor vj = target(edge, g);

            vector<VertexDescriptor> betweenVertex;
            bool normalDir = false;

            if ((*guideG[gi].coordPtr - *g[vi].coordPtr).norm() <
                (*guideG[gn].coordPtr - *g[vi].coordPtr).norm())
            {
                normalDir = false;
            } else normalDir = true;

            if (highGuideIndex - lowGuideIndex < _guide->getNumEdges() - highGuideIndex + lowGuideIndex ) { // zerro is not inbetween
                BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                    if (guideG[vertex].id > lowGuideIndex && guideG[vertex].id < highGuideIndex ) {
                        betweenVertex.push_back(vertex);
                    }
                }
            } else {
                int countBetween0andLow = 0;
                int countBevore0 = 0;
                normalDir = false;
                cerr << "low guide index " << lowGuideIndex << endl;
                BGL_FORALL_VERTICES(vertex, guideG, UndirectedGraph ) {
                    if (guideG[vertex].id <= lowGuideIndex ) {
                        // auto pos = betweenVertex.end() - countBetween0andLow
                        betweenVertex.insert(betweenVertex.begin(), vertex);
                        countBetween0andLow++;
                    } else if ( guideG[vertex].id > highGuideIndex ) {
                        auto pos = betweenVertex.begin() + countBevore0;
                        betweenVertex.insert(pos, vertex);
                        countBevore0++;
                    }
                }
            }


            VertexDescriptor prevVert = vi;
            if (normalDir) {
                vi = vj;
                vj = prevVert;
                prevVert = vi;
            }

            Coord2 vi_dist = (g[vi].closestPointOnEdge - *g[vi].coordPtr);
            Coord2 vj_dist = (g[vj].closestPointOnEdge - *g[vj].coordPtr);
            g[vi].coordPtr->setX(g[vi].closestPointOnEdge.x());
            g[vi].coordPtr->setY(g[vi].closestPointOnEdge.y());
            g[vj].coordPtr->setX(g[vj].closestPointOnEdge.x());
            g[vj].coordPtr->setY(g[vj].closestPointOnEdge.y());


            int count = 1;
           for (auto it = betweenVertex.begin(); it != betweenVertex.end(); it++ ) {

                VertexDescriptor uid = *it;
                VertexDescriptor vij = add_vertex(g);
                double x = guideG[uid].coordPtr->x();
                double y = guideG[uid].coordPtr->y();
                double k;
                if (betweenVertex.size() != 0) k = count / betweenVertex.size();
                else k = 0;
                // x -= ( 1. - k ) * vi_dist.x() + ( k ) * vj_dist.x();
                // y -= ( 1. - k ) * vi_dist.y() + ( k ) * vj_dist.y();

                g[vij].coordPtr = new Coord2(x, y);
                g[vij].smoothPtr = new Coord2(x, y);
                g[vij].octilinearPtr = new Coord2(x, y);
                g[vij].intersection = false;
                g[vij].metroShape = true;
                g[vij].isStation = false;
                g[vij].inflectionPoint = false;
                g[vij].closestPointOnEdge = Coord2(1000, 1000);
                g[vij].id = _metro->getNumStations() + addVertex;
                addVertex++;

                EdgeDescriptor e1; bool b1;
                tie(e1,b1) = add_edge(prevVert, vij, g);
                g[e1].isSmoothPath = true;
                g[e1].isMetro = true;
                g[e1].isCurved = false;
                g[e1].matchPath = false;
                g[e1].lineID = g[edge].lineID;
                g[e1].originalMetroEdge = false;
                addEdges++;
                g[e1].id = _metro->getNumEdges() + addEdges;
                prevVert = vij;

                count++;

           }
            EdgeDescriptor e2; bool b2;
            tie(e2, b2) = add_edge(prevVert, vj, g);
            g[e2].isSmoothPath = true;
            g[e2].isMetro = true;
            g[e2].isCurved = false;
            g[e2].matchPath = false;
            g[e2].lineID = g[edge].lineID;
            g[e2].originalMetroEdge = false;
            addEdges++;
            g[e2].id = _metro->getNumEdges() + addEdges;

            // remove_edge(edge, g); // probably that was the problem
            g[edge].isVisible = false;
        }
        else {}
    }
    cerr << "end post processing" << endl;
}

/*
* retrieve the result
*
*  Inputs
*      none
*  Outputs
*      none
*/
void Mixedlayout::retrieve( void )
{
    UndirectedGraph        & g            = _metro->g();
    unsigned int nVertices      = _metro->nStations();

    // find the vertex that is too close to an edge
    vector< VertexDescriptor > vdVec;
    for ( VEMap::iterator it = _metro->VEconflict().begin();
          it != _metro->VEconflict().end(); ++it ) {
        VertexDescriptor vdV = it->second.first;
        vdVec.push_back( vdV );
    }


    unsigned int nRows = 0;
    /*
    // update coordinates
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        bool doClose = false;
        for( unsigned int i = 0; i < (int) vdVec.size(); i++ ){
            if( vertex == vdVec[i] ) doClose = true;
        }
        // if (doClose ) g[vertex].collision = true;
        if( doClose){
        // if( _metro->VEconflict().size() > 0 ){
            // Coord2 downscale;
            // downscale.x() = ( _var( nRows, 0 ) - g[ vertex ].coordPtr->x() )/2.0 + g[ vertex ].coordPtr->x();
            // downscale.y() = ( _var( nRows, 0 ) - g[ vertex ].coordPtr->x() )/2.0 + g[ vertex ].coordPtr->x();;
            // g[ vertex ].coordPtr->x() = downscale.x();
            // g[ vertex ].coordPtr->y() = downscale.y();

            g[ vertex ].coordPtr->x() = _var( nRows, 0 );
            g[ vertex ].coordPtr->y() = _var( nRows + nVertices, 0 );
            g[vertex].collision = true;
        }
        else{
            g[ vertex ].coordPtr->x() = _var( nRows, 0 );
            g[ vertex ].coordPtr->y() = _var( nRows + nVertices, 0 );
            g[vertex].collision = false;
        }
        nRows++;
    } */

    // check if two vertex are to close
    nRows = 0;
    map<VertexDescriptor, Coord2> oldPositions;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        Coord2 nextPos = Coord2(_var( nRows, 0 ), _var( nRows + nVertices, 0 ) );
        Coord2 newClosPoint = _metro->closestPointOnGraph(vertex, nextPos);

        oldPositions[vertex] = *g[vertex].coordPtr;
        if ((nextPos - newClosPoint).norm() < MIN_DISTANCE) {
            g[ vertex ].coordPtr->x() = g[ vertex ].coordPtr->x() * .5 + _var( nRows, 0 ) * .5;
            g[ vertex ].coordPtr->y() = g[ vertex ].coordPtr->y() * .5 + _var( nRows + nVertices, 0 ) * .5;
            g[vertex].collision = true;
        } else {
            g[ vertex ].coordPtr->x() = _var( nRows, 0 );
            g[ vertex ].coordPtr->y() = _var( nRows + nVertices, 0 );
            g[vertex].collision = false;
        }
        nRows++;
    }
    int count = 0;
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph ){
        bool planarViolation = false;
        OutEdgeIterator e, e_end;
        for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
            EdgeDescriptor ed = *e;
            auto intersectingEdges = _metro->intersectingEdge(ed);
            if (!intersectingEdges.empty()){
                planarViolation = true;
                // cout << "   not planar anymore: " << intersectingEdges.size() << endl;
            }
            for (const auto& edInt : intersectingEdges) {
                auto sD = source(edInt.first, g);
                auto tD = target(edInt.first, g);
                g[ sD ].coordPtr->x() = oldPositions[sD].x();
                g[ tD ].coordPtr->x() = oldPositions[tD].x();
            }
        }
        if (planarViolation) {
            g[ vertex ].coordPtr->x() = oldPositions[vertex].x();
            g[ vertex ].coordPtr->y() = oldPositions[vertex].y();
        }

        count++;
    }
    _metro->checkVEConflicts();
}

void Mixedlayout::_checkSameAnge() {
    /*
    cout << "checking if some vertex twice the angle " << endl;
    UndirectedGraph & g = _metro->g();
    // Check if for a vertex there is a outgoing edge that has the same target angle
    BGL_FORALL_VERTICES( vertex, g, UndirectedGraph) {
        DegreeSizeType degrees = out_degree( vertex, g );
        vector<double> targetAngles;
        bool conflict = false;
        if (degrees > 1 && !g[vertex].smoothPath ) {
            OutEdgeIterator e, e_end;
            for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                EdgeDescriptor ed = *e;
                double ta = g[ed].target;

                VertexDescriptor vj = target(ed, g);
                if (g[vj].id == g[vertex].id)
                    vj = source(ed,g);

                if (g[vertex].id > g[vj].id) {
                    ta += M_PI;

                    if (ta > M_PI)
                        ta = - 2.0 * M_PI + ta; // Bring to a negative sector

                }


                if (find(targetAngles.begin(), targetAngles.end(), ta) != targetAngles.end()) { // not sure if optimal
                    conflict = true;
                }
                targetAngles.push_back(ta);
            }
        }

        // Do something In case there is a conflict
        if (conflict && false) {
            cout << "   There is a conflict around vertex: " << g[vertex].id << ", " << *g[vertex].namePtr << endl;
            vector<EdgeDescriptor> edges;
            vector<VertexDescriptor> neighbors;
            vector<double> curAngles;

            Coord2 vi_cor = *g[vertex].coordPtr;

            OutEdgeIterator e, e_end;
            for ( tie( e, e_end ) = out_edges( vertex, g ); e != e_end; ++e ) {
                EdgeDescriptor edge = *e;
                edges.push_back(edge);
                VertexDescriptor neighbor = source(edge, g);
                if ( neighbor == vertex )
                    neighbor = target(edge, g);
                neighbors.push_back(neighbor);

                Coord2 vj_cor = *g[neighbor].coordPtr;
                Coord2 vij_cor = Coord2(vj_cor.x() - vi_cor.x(),
                                        vj_cor.y() - vi_cor.y());
                double angle = atan2(vij_cor.y(), vij_cor.x());
               curAngles.push_back(angle);
               // curAngles.push_back(g[edge].curAngle);
            }

            // Solve Issue:
            auto edgesSectors = _solveTargetAnglesConflict(edges, curAngles);
            cout << "    sectors around" << *g[vertex].namePtr << endl;
            for (auto p : edgesSectors) {

                EdgeDescriptor ed = p.first;
                auto sector = p.second;

                VertexDescriptor vj = target(ed, g);
                if (g[vj].id == g[vertex].id)
                    vj = source(ed,g);

                if (g[vertex].id > g[vj].id) {
                    sector -= M_PI;

                    if (sector > M_PI) {
                        sector = - 2.0 * M_PI + sector; // Bring to a negative sector
                    }
                    if (sector < -M_PI)
                        sector = 2.0 * M_PI + sector;

                }

                // if (g[vertex].id > g[vj].id)
                //  sector -= M_PI;

                // if (sector >= M_PI)
                //     sector = 2.0 * M_PI - sector;
                // if (sector <= -M_PI)
                //     sector = 2.0 * M_PI + sector;

                cout << "      " << *g[vj].namePtr << ": " <<  sector << endl;

                g[ed].target = sector;
            }
        }
    }
     */
}


vector<pair<EdgeDescriptor, double>> Mixedlayout::_solveTargetAnglesConflict(vector<EdgeDescriptor> edges, vector<double> curAngles) {
    vector<pair<EdgeDescriptor, double>> result;
    /*
    vector<double> sectors { // -M_PI,
        -3.0*M_PI/4.0, -M_PI/2.0, -M_PI/4.0, 0.0,
        M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI
    };

    UndirectedGraph & g = _metro->g();
    cout << "solve asignment for angles: " << endl;
    for (auto a : curAngles)
        cout << "   " << a << endl;
    // create all Permuations
    vector<vector<pair<EdgeDescriptor, double>>> permutations;
    vector<vector<double>> permut_costs;
    for (int i = 0; i < edges.size(); i++ ) {
        auto ed = edges[i];
        auto edge_angle = curAngles[i];
        // if (edge_angle < 0.0)
        //     edge_angle += 2.0 * M_PI;

        vector<pair<EdgeDescriptor, double>> allForEd;
        vector<double> costs;

        for (auto a : sectors) {
            double sector_angle = a;

            double sec_a = a;
            double edg_a = edge_angle;

            double cost = min(abs(sec_a - edg_a),
                              abs((M_PI - sec_a) - (M_PI - edg_a)));
            // cost *= degrees_source * degrees_target;
            allForEd.push_back(pair<EdgeDescriptor, double>(ed, a));
            costs.push_back(cost);
        }
        permutations.push_back(allForEd);
        permut_costs.push_back(costs);
    }

    // Select all valid rotations from Permutations
    vector<vector<pair<EdgeDescriptor, double>>> rotations;
    vector<vector<double>> rotations_costs;
    for (int ja = 0; ja < sectors.size(); ja++ ) {

        vector<pair<EdgeDescriptor, double>> rot;
        vector<double> costs;

        for (int i = 0; i < edges.size(); i++ ) {
            int index = ( i + ja ) % sectors.size();
            rot.push_back(pair<EdgeDescriptor, double>(permutations[i][index]));
            costs.push_back(permut_costs[i][index]);
        }
        rotations.push_back(rot);
        rotations_costs.push_back(costs);
    }

    // Select the best rotation
    double minCosts = 1000000.0;
    for (int i = 0; i < rotations.size(); i++) {
        auto rotation = rotations[i];
        auto costs = rotations_costs[i];
        double sumCosts = 0.0;
        cout << "rotation: \n" << "   ";
        for (int j = 0; j < rotation.size(); j++) {
            sumCosts += costs[j];
            cout << rotation[j].second << "/" << costs[j] << ", ";
        }
        cout << "; with costs " << sumCosts << endl;

        if (minCosts > sumCosts) {
            cout << "new min costs " << sumCosts << endl;
            minCosts = sumCosts;
            result = rotation;
        }
    }
    cout << "    changed to: " << endl;
    for (auto r : result)
        cout << "      " << r.second << endl;
    cout << "    with costs: " << minCosts << endl;
    */
     return result;

}

/*
* Calculates magnitude of a vector from [0,0] to v
*
* Inputs    Coordinates of point
* Outputs   Magnitude
*/
double Mixedlayout::magnitude( Coord2 v ) {
    double d = v.x() * v.x() + v.y() * v.y(); 
    return sqrt(d);
}

/*
* memory management
*
*  Inputs
*      none
*
*  Outputs
*      none
*/
void Mixedlayout::clear( void ){}

//------------------------------------------------------------------------------
//	Public functions
//------------------------------------------------------------------------------

/*
* default constructor
*
*  Inputs
*     none
*  Outputs
*     none
*/
Mixedlayout::Mixedlayout( void ) {}

/*
* copy constructor
*
*  Inputs
*    obj : object of this class
*  Outputs
*     none
*/
Mixedlayout::Mixedlayout( const Mixedlayout & obj )
{
}

/*
* destructor
*
*  Inputs
*     none
*  Outputs
*     none
*/
Mixedlayout::~Mixedlayout( void ) {}

// end of header file
// Do not add any stuff under this line.

