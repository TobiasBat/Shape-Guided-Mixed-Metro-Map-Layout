/*
*	Guid.h
*   Based on FocusContext/Smooth.cpp
*	@author Tobias Batik 
*	@version 1.0  14/03/2021
*
*   Guide Path / Guide Graph 
*/

#ifndef _Guid_H        // begining of header file
#define _Guid_H        // notifying that this file is included

#include <iostream>
#include <iomanip>


using namespace std;

#include "BaseUndirectedGraph.h"
#include "Grid2.h"

typedef Graph::BaseUndirectedGraph                      UndirectedGraph;
typedef Graph::BaseUndirectedGraph::vertex_descriptor   VertexDescriptor;
typedef Graph::BaseUndirectedGraph::edge_descriptor     EdgeDiscriptor;

#define MAX_LINES		(32)
#define MIN_NODE_DIST  (2.0)

class Guide {
    private: 
        UndirectedGraph     graph;

        int curVert; 
        double magnitude( Coord2 v );
        void printGraph( void );
        VertexDescriptor setNewVertex(double x, double y, vector < VertexDescriptor > &ptrPoints, int _id);
        vector<VertexDescriptor> outerPath;
        bool complexGuide;

    protected: 
        int _nVertex; 
        int _nEdges;
        unsigned int _nLines;   // TODO
        vector< vector< VertexDescriptor > >	_lineSta;   // TODO
        string					    _lineName   [ MAX_LINES ];   // TODO //was char array


    public: 
        Guide();
        Guide( const Guide & obj ); 
        virtual ~Guide();  

        void                    load( string filePath ); 
        EdgeDiscriptor          closestEdge(Coord2 vertex);
        UndirectedGraph&        g( void ) { return graph; } 
        Coord2                  closestPointOnEdge( Coord2 source, Coord2 target, Coord2 vertex ); 
        Coord2                  closestPointOnEdge( EdgeDiscriptor edge, Coord2 vertex ); 
        Coord2                  getCenterCoord(); 
        double                  distanceToClosestPointOnEdge( EdgeDiscriptor edge, Coord2 vertex );     
        int                     getNumVertex( void ) { return _nVertex; }
        int                     getNumEdges( void ) { return _nEdges; }
        void exportGraphML( const string filename );
        void                    translate(Coord2 trans);
        void                    scale(double s);
        vector<VertexDescriptor>    getOutherPath() { return outerPath; }
        bool                    isComplex() { return complexGuide; }
        vector<VertexDescriptor> shortesPath(VertexDescriptor ui, VertexDescriptor uj);
        double                  getHeight();
        double                  getWidth();
};

#endif