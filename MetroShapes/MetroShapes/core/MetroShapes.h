/*
*	MetroShapes.h
*	@author Tobias Batik
*	@version 1.0  14/03/2021
*/

#ifndef _FocusContext_H
#define _FocusContext_H

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

#include "Base.h"
#include "Metro.h"
#include "Guide.h" 
#include "../optimization/Smooth.h"
#include "../optimization/Mixedlayout.h"
#include "../matching/Stationmatching.h"
#include "../matching/AutoMatching.h"
#include "../matching/Manuelpath.h"

class MetroShapes : public Base
{
	//------------------------------------------------------------------------------
	//  Defining Classes
	//------------------------------------------------------------------------------
	
	private:
		
		Stationmatching _matching;
		AutoMatching _autoMatching;
		Metro _metro;
		string _inputname;
		string _outputname;
		string _guidename;
		Guide _guide;
		Manuelpath _manuelPath;
        bool _automaticCase;

		bool _exportFrames = false;
		int _stepsSmooth = 100000;
		int _stepsMixed  = 100000;
	protected:

	public:
		
		Smooth 			_smooth;
		Mixedlayout 	_mixedlayout;

		//------------------------------------------------------------------------------
		//      Constructors & Destructors
		//------------------------------------------------------------------------------
		// default constructor
		MetroShapes( void ) {}
		// copy constructor
		MetroShapes( const MetroShapes & v ) {}
		// destructor
		virtual ~MetroShapes( void ) {}
		
		//------------------------------------------------------------------------------
		//      Special functions
		//------------------------------------------------------------------------------
		void init( int argc, char **argv ) override;
		void run( void ) override;
		void output( void ) override;
		void graphMatching( void ); 
		void computeSmooth( void );
		void computeMixedlayout( void );
		void computeMatching( void );
		void computeAutoMatching();
		void distortMetroMap( void );
		void scaleNonUniformlyBasedOnManuelPath();
		void resetSmooth( void ); 
		bool addNode( Coord2 ); 
		bool removeNode( Coord2 );
        void selectPathFromMetroLine(int l );
        void addVertexToSelectPath( int id);
        void exportSteps();

        void manuallyAddVertexToPath( Coord2 coord);


        //------------------------------------------------------------------------------
		//      Reference to elements
		//------------------------------------------------------------------------------
		Metro &getMetro( void )                 { return _metro; }
		Guide &getGuid( void )					{ return _guide; }
		AutoMatching &getAutoMatching(void)     { return _autoMatching; }
		void setMetro( const Metro m )          { _metro = m; }
		Stationmatching &getMatching( void ) 	{ return _matching; }
		
		
		//------------------------------------------------------------------------------
		//      I/O functions
		//------------------------------------------------------------------------------
		// output
		friend ostream & operator << ( ostream & stream, const MetroShapes & obj );
		// input
		friend istream & operator >> ( istream & stream, MetroShapes & obj );
		// class name
		virtual const char * className( void ) const { return "MetroShapes"; }

    void pathFromMetroLine(int l);
};

#endif // _MetroShapes_H
