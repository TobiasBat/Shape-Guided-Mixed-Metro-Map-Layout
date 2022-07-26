//==============================================================================
// Label.cc
//	: program file for annotation labels
//
//------------------------------------------------------------------------------
//
//				Date: Sun Jul 22 03:01:14 2012
//
//==============================================================================

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "Common.h"
#include "Label.h"


//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Private Functions
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Protected Functions
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Public Functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//	Constructors
//------------------------------------------------------------------------------

//
//  Label::Label --	default constructor
//
//  Inputs
//	none
//
//  Outputs
//	none
//
Label::Label() 
{
    _id = NO_UNSIGNED_ID;
    _width = 0.0;
    _height = 0.0;
    _flag = false;
    _leaderW = 1.0;
    _geoSite.zero();
    _curSite.zero();
    _joint.zero();
    _leftTop.zero();
    _rightBottom.zero();
    _frame.zero();
    _base.zero();
}


//
//  Label::Label --	copy constructor
//
//  Inputs
//	obj	: object of this class
//
//  Outputs
//	none
//
Label::Label( const Label & obj )
{
    _id		= obj._id;
    _width	= obj._width;
    _height	= obj._height;
    _flag	= obj._flag;
    _geoSite	= obj._geoSite;
    _curSite	= obj._curSite;
    _leaderW	= obj._leaderW;
    _joint	= obj._joint;
    _leftTop	= obj._leftTop;
    _rightBottom= obj._rightBottom;
    _frame	= obj._frame;
    _base	= obj._base;
}


//------------------------------------------------------------------------------
//	Destructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Referring to members
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//	Functions for placing labels
//------------------------------------------------------------------------------

//
//  Label::minimizeLeader --	minimize the leader length by searching for the
//				best joint position 
//
//  Inputs
//	variable : explanation of a variable
//	variable : explanation of a variable
//
//  Outputs
//	leader length
//
double Label::minimizeLeader( const Coord2 & curSite, const Coord2 & corner, const double & interval,
			      Coord2 & coord )
{
    double minDist = INFINITY;
    Coord2 minCoord;

    int i, j;
    const int m = 5;

//------------------------------------------------------------------------------
//	top horizontal frame segment
//------------------------------------------------------------------------------
    j = 0;
    for ( i = m; i <= _frame.p()*10-m; ++i ) {
	double ii = ( double )i/10.0;
	Coord2 curCoord = corner + Coord2( ii*interval, (double)-j*interval );
	double curDist = distanceBetween( curSite, curCoord );
#ifdef DEBUG
	cerr << "[T] i = " << i << " j = " << j << " dist = " << curDist << " coord = " << curCoord;
#endif	// DEBUG
	if ( curDist < minDist ) {
	    minDist = curDist;
	    minCoord = curCoord;
	}
    }
//------------------------------------------------------------------------------
//	bottom horizontal frame segment
//------------------------------------------------------------------------------
    j = _frame.q();
    for ( i = m; i <= _frame.p()*10-m; ++i ) {
	double ii = ( double )i/10.0;
	Coord2 curCoord = corner + Coord2( ii*interval, (double)-j*interval );
	double curDist = distanceBetween( curSite, curCoord );
#ifdef DEBUG
	cerr << "[B] i = " << i << " j = " << j << " dist = " << curDist << " coord = " << curCoord;
#endif	// DEBUG
	if ( curDist < minDist ) {
	    minDist = curDist;
	    minCoord = curCoord;
	}
    }
//------------------------------------------------------------------------------
//	left vertical frame segment
//------------------------------------------------------------------------------
    i = 0;
    for ( j = m; j <= _frame.q()*10-m; ++j ) {
	double jj = ( double )j/10.0;
	Coord2 curCoord = corner + Coord2( (double)i*interval, -jj*interval );
	double curDist = distanceBetween( curSite, curCoord );
#ifdef DEBUG
	cerr << "[L] i = " << i << " j = " << j << " dist = " << curDist << " coord = " << curCoord;
#endif	// DEBUG
	if ( curDist < minDist ) {
	    minDist = curDist;
	    minCoord = curCoord;
	}
    }
//------------------------------------------------------------------------------
//	right vertical frame segment
//------------------------------------------------------------------------------
    i = _frame.p();
    for ( j = m; j <= _frame.q()*10-m; ++j ) {
	double jj = ( double )j/10.0;
	Coord2 curCoord = corner + Coord2( (double)i*interval, (double)-jj*interval );
	double curDist = distanceBetween( curSite, curCoord );
#ifdef DEBUG
	cerr << "[R] i = " << i << " j = " << j << " dist = " << curDist << " coord = " << curCoord;
#endif	// DEBUG
	if ( curDist < minDist ) {
	    minDist = curDist;
	    minCoord = curCoord;
	}
    }

    coord = minCoord;
    return minDist;
}


//
//  Label::isOverlapped --	check whether the input line segment intersects
//				with the label or not.
//
//  Inputs
//	orig, dest	: endpoints of the segment
//
//  Outputs
//	boolean value according to the existence of intersection
//
const bool Label::isOverlapped( const Coord2 & orig, const Coord2 & dest ) const
{
    if ( isIntersected( orig, dest, _curSite, _joint ) ) return true;
    const Coord2 & lt = _leftTop;
    const Coord2 & rb = _rightBottom;
    Coord2 lb = Coord2( _leftTop.x(), _rightBottom.y() );
    Coord2 rt = Coord2( _rightBottom.x(), _leftTop.y() );
    //if ( isIntersected( orig, dest, lt, rt ) ) return true;
    //if ( isIntersected( orig, dest, lt, lb ) ) return true;
    //if ( isIntersected( orig, dest, rb, rt ) ) return true;
    //if ( isIntersected( orig, dest, rb, lb ) ) return true;
    if ( doConflict( orig, dest, lt, rt ) ) return true;
    if ( doConflict( orig, dest, lt, lb ) ) return true;
    if ( doConflict( orig, dest, rb, rt ) ) return true;
    if ( doConflict( orig, dest, rb, lb ) ) return true;

    return false;
}


//
//  Label::isOverlapped --	check whether the input label intersects
//				with the label or not.
//
//  Inputs
//	box	: label
//
//  Outputs
//	boolean value according to the existence of intersection
//
const bool Label::isOverlapped( const Label & box ) const
{
    const double small = SHRINKAGE_RATIO;

    const Coord2 & lt = _leftTop;
    const Coord2 & rb = _rightBottom;
    const Coord2 & boxLT = box.lt();
    const Coord2 & boxRB = box.rb();

    if ( ( boxRB.x() <= lt.x() ) ||
	 ( rb.x() <= boxLT.x() ) ||
	 ( boxLT.y() <= rb.y() ) ||
	 ( lt.y() <= boxRB.y() ) ) {
	;			// do nothing
    }
    else return true;

    Coord2 boxO = ( 1.0-small )*box.curSite() + small*box.joint();
    Coord2 boxD = small*box.curSite()+( 1.0-small ) * box.joint();

    if ( isOverlapped( boxO, boxD ) ) return true;

    Coord2 orig = ( 1.0-small )*_curSite + small*_joint;
    Coord2 dest = small*_curSite + ( 1.0-small )*_joint;

    if ( box.isOverlapped( orig, dest ) ) return true;

    return false;
}


//------------------------------------------------------------------------------
//	Assignment opereators
//------------------------------------------------------------------------------

//
//  Label::operator = --	assignement
//				代入
//
//  Inputs
//	obj	: objects of this class
//
//  Outputs
//	this object
//
Label & Label::operator = ( const Label & obj )
{
    if ( this != &obj ) {
	_id	= obj._id;
	_width	= obj._width;
	_height	= obj._height;
	_flag	= obj._flag;
	_geoSite= obj._geoSite;
	_curSite= obj._curSite;
	_leaderW= obj._leaderW;
	_joint	= obj._joint;
	_leftTop	= obj._leftTop;
	_rightBottom	= obj._rightBottom;
	_frame	= obj._frame;
	_base	= obj._base;
    }
    return *this;
}


//------------------------------------------------------------------------------
//	I/O functions
//------------------------------------------------------------------------------



// end of header file
// Do not add any stuff under this line.
