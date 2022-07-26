#ifndef _Label_H
#define _Label_H

#include <vector>
#include "Coord2.h"
#include "Grid2.h"

using namespace std;


class Label
{
private:
    unsigned int        _id;
    double              _width;
    double              _height;
    double              _leaderW;
    bool                _flag;
    Coord2              _curSite;
    Coord2              _geoSite;
    Coord2              _joint;
    Coord2              _leftTop;
    Coord2              _rightBottom;
    Grid2               _frame;
    Grid2               _base;
#ifdef  REWRITE
    VertexDescriptor    _ptrVertex;
#endif  // REWRITE

public:


    Label();
    Label( const Label & obj );

//------------------------------------------------------------------------------
//      Handling member variables
//------------------------------------------------------------------------------

    unsigned int & id( void )                   { return _id; }
    const unsigned int & id( void )  const      { return _id; }
    const unsigned int & getID( void ) const    { return _id; }
    void setID( unsigned int __id )             { _id = __id; }


    double & width( void )                     { return _width; }
    const double & width( void ) const         { return _width; }

    double & height( void )                     { return _height; }
    const double & height( void ) const         { return _height; }

    double & leaderWeight( void )               { return _leaderW; }
    const double & leaderWeight( void ) const   { return _leaderW; }

    const bool & flag( void ) const             { return _flag; }
    void setFlag( void )                        { _flag = true; }
    void clearFlag( void )                      { _flag = false; }

    Coord2 & geoSite( void )                    { return _geoSite; }
    const Coord2 & geoSite( void ) const        { return _geoSite; }
    void setGeoSite( const Coord2 & __geoSite ) { _geoSite = __geoSite; }

    Coord2 & curSite( void )                    { return _curSite; }
    const Coord2 & curSite( void ) const        { return _curSite; }
    void setCurSite( const Coord2 & __curSite ) { _curSite = __curSite; }

    Coord2 & joint( void )                      { return _joint; }
    const Coord2 & joint( void ) const          { return _joint; }
    void setJoint( const Coord2 & __joint )     { _joint = __joint; }

    Coord2 & lt( void )                         { return _leftTop; }
    const Coord2 & lt( void ) const             { return _leftTop; }
    void setLT( const Coord2 & __lt )           { _leftTop = __lt; }

    Coord2 & rb( void )                         { return _rightBottom; }
    const Coord2 & rb( void ) const             { return _rightBottom; }
    void setRB( const Coord2 & __rb )           { _rightBottom = __rb; }

    Grid2 & frame( void )                       { return _frame; }
    const Grid2 & frame( void ) const           { return _frame; }
    void setFrame( const Grid2 & __frame )      { _frame = __frame; }

    Grid2 & base( void )                        { return _base; }
    const Grid2 & base( void ) const            { return _base; }
    void setBase( const Grid2 & __base )        { _base = __base; }

#ifdef  REWRITE
    VertexDescriptor ptrVertex( void )                  { return _ptrVertex; }
    void setPtrVertex( VertexDescriptor __ptrVertex )   { _ptrVertex = __ptrVertex; }
#endif  // REWRITE

//------------------------------------------------------------------------------
//      Specific functions
//------------------------------------------------------------------------------

    double minimizeLeader( const Coord2 & curSite, const Coord2 & corner, const double & interval,
                           Coord2 & coord );

    const bool isOverlapped( const Coord2 & orig, const Coord2 & dest ) const;

    const bool isOverlapped( const Label & box ) const;

//------------------------------------------------------------------------------
//      Assignment opereators
//------------------------------------------------------------------------------

    Label & operator = ( const Label & obj );

};

#endif // _Label_H
