#ifndef _Base_H
#define _Base_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "Coord2.h"

using namespace std;

class Base
{
	//------------------------------------------------------------------------------
	//  Defining Classes
	//------------------------------------------------------------------------------
	
private:

protected:

public:
	
	//------------------------------------------------------------------------------
	//      Constructors & Destructors
	//------------------------------------------------------------------------------
	// default constructor
	Base( void ) {}
	// copy constructor
	Base( const Base & v ) {}
	// destructor
	virtual ~Base( void ) {}
	
	//------------------------------------------------------------------------------
	//      Special functions
	//------------------------------------------------------------------------------
	virtual void init( int argc, char **argv ) = 0;
	virtual void run( void ) = 0;
	virtual void output( void ) = 0;

	// Added For MetroShapes
	virtual void computeAutoMatching( void ) = 0;
	virtual void computeMatching( void ) = 0; 
	virtual void computeSmooth( void ) = 0; 
	virtual void computeMixedlayout( void ) = 0;
	virtual void distortMetroMap( void ) = 0;
	virtual void scaleNonUniformlyBasedOnManuelPath() = 0;
	virtual void resetSmooth( void ) = 0;
	virtual bool removeNode( Coord2 ) = 0; 
	virtual bool addNode( Coord2 ) = 0;
	virtual void manuallyAddVertexToPath( Coord2 coord) = 0;
	virtual void selectPathFromMetroLine( int l ) = 0;
	
	//------------------------------------------------------------------------------
	//      I/O functions
	//------------------------------------------------------------------------------
	// output
	friend ostream & operator << ( ostream & stream, const Base & obj );
	// input
	friend istream & operator >> ( istream & stream, Base & obj );
	// class name
	virtual const char * className( void ) const { return "Base"; }
};

#endif // _Base_H
