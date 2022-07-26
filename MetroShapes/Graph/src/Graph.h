
#ifndef _Graph_H
#define _Graph_H


#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>
#include <ctime>
#include <cstdlib>

using namespace std;

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>

using namespace boost;

#include "Common.h"
#include "Coord2.h"
#include "Label.h"

//----------------------------------------------------------------------
//  Defining macros
//----------------------------------------------------------------------

#define NO_VERTEX_INDEX         (9999)


//----------------------------------------------------------------------
//  Defining BGL data structures
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//  Boost customization
//----------------------------------------------------------------------
//
//// Airport attributes
//enum vertex_myid_t			{ vertex_myid };		// Ordered Vertex ID without leaking number in graph
//enum vertex_mycoord_t                   { vertex_mycoord };             // current 2D coordinates
//enum vertex_mysmooth_t                  { vertex_mysmooth };            // smooth 2D coordinates
//enum vertex_mygeo_t                     { vertex_mygeo };               // geographical 2D coordinates
//enum vertex_mytemp_t                    { vertex_mytemp };              // temp position for animation
//enum vertex_myhome_t                    { vertex_myhome };              // home 2D coordinates
//enum vertex_mylineid_t                  { vertex_mylineid };            // line ID
//enum vertex_myname_t                    { vertex_myname };              // vertex name
//enum vertex_mysize_t                    { vertex_mysize };              // size
//enum vertex_myscale_t                   { vertex_myscale };             // scale
//enum vertex_myselectmag_t               { vertex_myselectmag };         // true if the vertex is selected as magnification lens area
//enum vertex_myselecttop_t               { vertex_myselecttop };         // true if the vertex is selected as topology lens area
//enum vertex_mygeodesic_t                { vertex_mygeodesic };          // record the geodesic distance from the selected node
//enum vertex_myzone_t                    { vertex_myzone };              // record the zone from the selected node
//enum vertex_mytextureid_t               { vertex_mytextureid };         // label texture id
//enum vertex_mytexname_t                 { vertex_mytexname };           // link name to the label
//enum vertex_myexternal_t                { vertex_myexternal };          // label
//enum vertex_myextstate_t                { vertex_myextstate };          // label state
//enum vertex_myisstation_t               { vertex_myisstation };         // set true if it is a station node
//enum vertex_myremove_t                  { vertex_myremove };            // removed vertex
//enum vertex_weight_t                    { vertex_weight };              // vertex weight
//
//// Number of flights
//enum edge_myid_t			{ edge_myid };
//enum edge_mygeoangle_t                  { edge_mygeoangle };            // geographical angle
//enum edge_mysmoangle_t                  { edge_mysmoangle };            // smooth angle
//enum edge_mycurangle_t                  { edge_mycurangle };            // current angle
//enum edge_mytarget_t                    { edge_mytarget };              // target octilinear angle
//enum edge_mylineid_t			{ edge_mylineid };
//enum edge_myselectshift_t		{ edge_myselectshift };         // true if the edge is selected with shift key
//enum edge_myselectctrl_t		{ edge_myselectctrl };          // true if the edge is selected with ctrl key
//enum edge_myisline_t		        { edge_myisline };              // true if the edge is metro lines
//enum edge_myisleader_t		        { edge_myisleader };            // true if the edge is metro leader
//
//namespace boost {
//    // vertex properties
//    BOOST_INSTALL_PROPERTY( vertex, myid );
//    BOOST_INSTALL_PROPERTY( vertex, mycoord );
//    BOOST_INSTALL_PROPERTY( vertex, mysmooth );
//    BOOST_INSTALL_PROPERTY( vertex, mygeo );
//    BOOST_INSTALL_PROPERTY( vertex, mytemp );
//    BOOST_INSTALL_PROPERTY( vertex, myhome );
//    BOOST_INSTALL_PROPERTY( vertex, mylineid );
//    BOOST_INSTALL_PROPERTY( vertex, myname );
//    BOOST_INSTALL_PROPERTY( vertex, mysize );
//    BOOST_INSTALL_PROPERTY( vertex, myscale );
//    BOOST_INSTALL_PROPERTY( vertex, myselectmag );
//    BOOST_INSTALL_PROPERTY( vertex, myselecttop );
//    BOOST_INSTALL_PROPERTY( vertex, mygeodesic );
//    BOOST_INSTALL_PROPERTY( vertex, myzone );
//    BOOST_INSTALL_PROPERTY( vertex, mytextureid );
//    BOOST_INSTALL_PROPERTY( vertex, mytexname );
//    BOOST_INSTALL_PROPERTY( vertex, myexternal );
//    BOOST_INSTALL_PROPERTY( vertex, myextstate );
//    BOOST_INSTALL_PROPERTY( vertex, myisstation );
//    BOOST_INSTALL_PROPERTY( vertex, myremove );
//    BOOST_INSTALL_PROPERTY( vertex, weight );
//
//    // edge properties
//    BOOST_INSTALL_PROPERTY( edge,   myid );
//    BOOST_INSTALL_PROPERTY( edge,   mygeoangle );
//    BOOST_INSTALL_PROPERTY( edge,   mysmoangle );
//    BOOST_INSTALL_PROPERTY( edge,   mycurangle );
//    BOOST_INSTALL_PROPERTY( edge,   mytarget );
//    BOOST_INSTALL_PROPERTY( edge,   mylineid );
//    BOOST_INSTALL_PROPERTY( edge,   myselectshift );
//    BOOST_INSTALL_PROPERTY( edge,   myselectctrl );
//    BOOST_INSTALL_PROPERTY( edge,   myisline );
//    BOOST_INSTALL_PROPERTY( edge,   myisleader );
//}
//
////------------------------------------------------------------------------------
////  Customizing vertex properties
////------------------------------------------------------------------------------
//typedef property< vertex_myid_t, unsigned int >					MyVID;
//typedef property< vertex_mycoord_t, Coord2, MyVID >			        MyVCoord;
//typedef property< vertex_mysmooth_t, Coord2, MyVCoord >				MyVSmooth;
//typedef property< vertex_mygeo_t, Coord2, MyVSmooth >				MyVGeo;
//typedef property< vertex_mytemp_t, Coord2, MyVGeo >				MyVTemp;
//typedef property< vertex_myhome_t, Coord2, MyVTemp >				MyVHome;
//typedef property< vertex_mylineid_t, vector< unsigned int >, MyVHome >		MyVLineID;
//typedef property< vertex_myname_t, string, MyVLineID >			        MyVName;
//typedef property< vertex_mysize_t, double, MyVName >			        MyVSize;
//typedef property< vertex_myscale_t, double, MyVSize >			        MyVScale;
//typedef property< vertex_myselectmag_t, bool, MyVScale >			MyVSelectMag;
//typedef property< vertex_myselecttop_t, bool, MyVSelectMag >			MyVSelectTop;
//typedef property< vertex_mygeodesic_t, double, MyVSelectTop >			MyVGeodesic;
//typedef property< vertex_myzone_t, int, MyVGeodesic >			        MyVZone;
//typedef property< vertex_mytextureid_t, unsigned int, MyVZone >		        MyVTextureID;
//typedef property< vertex_mytexname_t, string, MyVTextureID >	                MyVTexName;
//typedef property< vertex_myexternal_t, Label, MyVTexName >			MyVExternal;
//typedef property< vertex_myextstate_t, bool, MyVExternal >			MyVExtstate;
//typedef property< vertex_myisstation_t, bool, MyVExtstate >			MyVIsStation;
//typedef property< vertex_myremove_t, bool, MyVIsStation >			MyVRemove;
//typedef property< vertex_weight_t, double, MyVRemove >			        MyVWeight;
//typedef property< vertex_index_t, unsigned int, MyVWeight >			MyVertexProperty;
//
////------------------------------------------------------------------------------
////  Customizing edge properties
////------------------------------------------------------------------------------
//typedef property< edge_myid_t, unsigned int >					MyEID;
//typedef property< edge_mygeoangle_t, double, MyEID >                            MyEGeoAngle;
//typedef property< edge_mysmoangle_t, double, MyEGeoAngle >                      MyESmoAngle;
//typedef property< edge_mycurangle_t, double, MyESmoAngle >                      MyECurAngle;
//typedef property< edge_mytarget_t, double, MyECurAngle >                        MyETarget;
//typedef property< edge_mylineid_t, vector< unsigned int >, MyETarget >          MyELineID;
//typedef property< edge_myselectshift_t, bool, MyELineID >                       MyESelectShift;
//typedef property< edge_myselectctrl_t, bool, MyESelectShift >                   MyESelectCtrl;
//typedef property< edge_myisline_t, bool, MyESelectCtrl >                        MyEIsLine;
//typedef property< edge_myisleader_t, bool, MyEIsLine >                          MyEIsLeader;
//typedef property< edge_weight_t, double, MyEIsLeader >                          MyEWeight;
//typedef property< edge_index_t, unsigned int, MyEWeight >			MyEdgeProperty;
//
//
//typedef adjacency_list< vecS, listS, undirectedS,
//                        MyVertexProperty, MyEdgeProperty >  Graph;
//typedef graph_traits< Graph >::degree_size_type             degree_size_type;
//typedef graph_traits< Graph >::adjacency_iterator           AdjacencyIterator;
//
//typedef graph_traits< Graph >                               GraphTraits;
//typedef graph_traits< Graph >::vertex_descriptor            VertexDescriptor;
//typedef graph_traits< Graph >::edge_descriptor              EdgeDescriptor;
//typedef pair< VertexDescriptor, VertexDescriptor >          VVPair;
//
//typedef graph_traits< Graph >::vertex_iterator              VertexIterator;
//
//typedef graph_traits< Graph >::edge_iterator                EdgeIterator;
//typedef graph_traits< Graph >::out_edge_iterator            OutEdgeIterator;
//typedef graph_traits< Graph >::in_edge_iterator             InEdgeIterator;
//typedef graph_traits< Graph >::edges_size_type		    EdgesSizeType;
//typedef graph_traits< Graph >::degree_size_type		    DegreeSizeType;
//
//
////------------------------------------------------------------------------------
////  Customized vertex maps
////------------------------------------------------------------------------------
//// Ordered Vertex ID without leaking number in graph
//typedef property_map< Graph, vertex_index_t >::type		VertexIndexMap;
//typedef property_map< Graph, vertex_myid_t >::type		VertexIDMap;
//// current 2D coordinates
//typedef property_map< Graph, vertex_mycoord_t >::type           VertexCoordMap;
//// smooth 2D coordinates
//typedef property_map< Graph, vertex_mysmooth_t >::type          VertexSmoothMap;
//// geographical 2D coordinates
//typedef property_map< Graph, vertex_mygeo_t >::type             VertexGeoMap;
//// temp 2D coordinates
//typedef property_map< Graph, vertex_mytemp_t >::type            VertexTempMap;
//// home 2D coordinates
//typedef property_map< Graph, vertex_myhome_t >::type            VertexHomeMap;
//// array of metro line IDs
//typedef property_map< Graph, vertex_mylineid_t >::type          VertexLineIDMap;
//typedef property_map< Graph, vertex_myname_t >::type            VertexNameMap;
//typedef property_map< Graph, vertex_mysize_t >::type            VertexSizeMap;
//typedef property_map< Graph, vertex_myscale_t >::type           VertexScaleMap;
//typedef property_map< Graph, vertex_myselectmag_t >::type       VertexSelectMagMap;
//typedef property_map< Graph, vertex_myselecttop_t >::type       VertexSelectTopMap;
//typedef property_map< Graph, vertex_mygeodesic_t >::type        VertexGeodesicMap;
//typedef property_map< Graph, vertex_myzone_t >::type            VertexZoneMap;
//typedef property_map< Graph, vertex_mytextureid_t >::type       VertexTextureIDMap;
//typedef property_map< Graph, vertex_mytexname_t >::type         VertexTexNameMap;
//typedef property_map< Graph, vertex_myexternal_t >::type        VertexExternalMap;
//typedef property_map< Graph, vertex_myextstate_t >::type        VertexExtstateMap;
//typedef property_map< Graph, vertex_myisstation_t >::type       VertexIsStationMap;
//// mark 1 if the station is removed
//typedef property_map< Graph, vertex_myremove_t >::type          VertexRemoveMap;
//typedef property_map< Graph, vertex_weight_t >::type            VertexWeightMap;
//typedef property_map< Graph, vertex_index_t >::type             VertexIndexMap;
//
//
////------------------------------------------------------------------------------
////  Customized edge maps
////------------------------------------------------------------------------------
//typedef property_map< Graph, edge_index_t >::type		EdgeIndexMap;
//typedef property_map< Graph, edge_myid_t >::type		EdgeIDMap;
//// array of metro line IDs
//typedef property_map< Graph, edge_mygeoangle_t >::type		EdgeGeoAngleMap;
//typedef property_map< Graph, edge_mysmoangle_t >::type		EdgeSmoAngleMap;
//typedef property_map< Graph, edge_mycurangle_t >::type		EdgeCurAngleMap;
//typedef property_map< Graph, edge_mytarget_t >::type		EdgeTargetMap;
//typedef property_map< Graph, edge_mylineid_t >::type		EdgeLineIDMap;
//typedef property_map< Graph, edge_myselectshift_t >::type	EdgeSelectShiftMap;
//typedef property_map< Graph, edge_myselectctrl_t >::type	EdgeSelectCtrlMap;
//typedef property_map< Graph, edge_myisline_t >::type	        EdgeIsLineMap;
//typedef property_map< Graph, edge_myisleader_t >::type	        EdgeIsLeaderMap;
//typedef property_map< Graph, edge_weight_t >::type		EdgeWeightMap;
//typedef property_map< Graph, edge_index_t >::type		EdgeIndexMap;
//
//// define the storage type for the planar embedding
//typedef vector< vector< EdgeDescriptor > >			EmbeddingStorage;
//typedef boost::iterator_property_map < EmbeddingStorage::iterator,
//                                       VertexIndexMap >		Embedding;
//
////------------------------------------------------------------------------------
////  Customized vertex function
////------------------------------------------------------------------------------
//void setVertexTexName( Graph & g, VertexDescriptor vd, const string & str );
//
////------------------------------------------------------------------------------
////	Customized Graph Functions
////------------------------------------------------------------------------------
//void printGraph( Graph & g );
//void clearGraph( Graph & g );
//void resetGraphCoord( Graph & g );
//void resetGraphWeight( Graph & g );
//int _longestLabelNum( VertexDescriptor vd, vector< unsigned int > & visit, Graph & g );
//void setGraphWeight( double x, double y, double meanNodeSize, double length, Graph & g );
//void shortestPath( VertexDescriptor sourceVD, VertexDescriptor targetVD, Graph & g );
//void geodesicDistance( VertexDescriptor focusVD, Graph & g );
//bool isEEOverlap( const Coord2 & As, const Coord2 & At,
//                  const Coord2 & Bs, const Coord2 & Bt );
#endif  // _Graph_H
