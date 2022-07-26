//******************************************************************************
// BaseUndirectedGraph.h
//	: header file for base undirected graph
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Dec 27 23:16:12 2018
//
//******************************************************************************

#ifndef _Graph_BaseUndirectedGraph_H
#define _Graph_BaseUndirectedGraph_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <cstdlib>

using namespace std;

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// force-directed layout
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>

using namespace boost;

#include "Coord2.h"
#include "BaseGraphProperty.h"
#include "BaseVertexProperty.h"
#include "BaseEdgeProperty.h"

//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace Graph {

    typedef adjacency_list< listS, listS, undirectedS,
            BaseVertexProperty, BaseEdgeProperty,
            BaseGraphProperty >  BaseUndirectedGraph;

    //------------------------------------------------------------------------------
    //	Special functions
    //------------------------------------------------------------------------------
    void printGraph( const BaseUndirectedGraph & g );
    void clearGraph( BaseUndirectedGraph & g );

} // namespace Graph

#endif  // _Graph_BaseUndirectedGraph_H
