//******************************************************************************
// BaseUndirectedGraph.cpp
//      : program file for graph function
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Mon Mar 14 02:16:23 2019
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "BaseUndirectedGraph.h"

namespace Graph {

    //------------------------------------------------------------------------------
    //	Customized Vertex Functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Customized Edge Functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Customized Layout Functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Special functions
    //------------------------------------------------------------------------------
    //
    //  BaseUndirectedGraph::printGraph -- print the graph.
    //
    //  Inputs
    //  g   : object of Grpah
    //
    //  Outputs
    //  none
    //
    void printGraph( const BaseUndirectedGraph & graph )
    {
        cerr << "num_vertices = " << num_vertices( graph ) << endl;
        cerr << "num_edges = " << num_edges( graph ) << endl;

    //#ifdef  DEBUG
        // print vertex information
        BGL_FORALL_VERTICES( vd, graph, BaseUndirectedGraph ) {

            //BaseUndirectedGraph::degree_size_type      degrees         = out_degree( vd, graph );
            cerr << " id = " << graph[vd].id
                 << " coord = " << *graph[ vd ].coordPtr;
        }
    //#endif  // DEBUG

    #ifdef  DEBUG
        // print edge information
        BGL_FORALL_EDGES( ed, graph, BaseUndirectedGraph ) {

            BaseUndirectedGraph::vertex_descriptor vdS = source( ed, graph );
            BaseUndirectedGraph::vertex_descriptor vdT = target( ed, graph );

            cerr << "eid = " << graph[ ed ].id << " ( " << graph[ vdS ].id << " == " << graph[ vdT ].id << " ) "
                 << " w = " << graph[ ed ].weight << endl;
        }
    #endif  // DEBUG
    }

    //
    //  BaseUndirectedGraph::clearGraph -- clear the graph.
    //
    //  Inputs
    //  g   : object of Graph
    //
    //  Outputs
    //  none
    //
    void clearGraph( BaseUndirectedGraph & graph )
    {
        // clear edges
        BaseUndirectedGraph::edge_iterator ei, ei_end, e_next;
        tie( ei, ei_end ) = edges( graph );
        for ( e_next = ei; ei != ei_end; ei = e_next ) {
            e_next++;
            remove_edge( *ei, graph );
        }

    #ifdef  SKIP
        BGL_FORALL_EDGES( edge, graph, BaseGraph )
        {
            remove_edge( edge, graph );
        }
    #endif  // SKIP

        // clear vertices
        pair< BaseUndirectedGraph::vertex_iterator, BaseUndirectedGraph::vertex_iterator > vp;
        for ( vp = vertices( graph ); vp.first != vp.second;  ) {
            BaseUndirectedGraph::vertex_descriptor vd = (*vp.first);
            ++vp.first;
            clear_vertex( vd, graph );
            remove_vertex( vd, graph );
        }
    }

} // namespace Graph