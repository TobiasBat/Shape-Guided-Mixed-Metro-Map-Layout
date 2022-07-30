/*
*	MetroShapes.cpp
*	@author Tobias Batik
*	@version 1.0  14/03/2021
*/

#include "MetroShapes.h"

void MetroShapes::init( int argc, char **argv )
{

    _outputname = "../../../output-optimisation/";


    if( argc == 1 ) {
		_inputname = "../data/metro/taipei-new.txt";
        _guidename = "../data/guide/circle-guide.txt";
        _outputname += "taipei";
        _automaticCase = true;
    }
    else if ( argc == 2 ) {
		_inputname = "../data/metro/"; 
		_inputname += + argv[1];
        _guidename = "../data/guide/circle-guide.txt";
        _outputname += "taipei";
        _automaticCase = true;
    }
    else if ( argc == 3 ) {
		_inputname = "../data/metro/";
		_inputname += + argv[1];
		_guidename = "../data/guide/";
		_guidename += argv[2];
        _outputname += "taipei";
        _automaticCase = true;
    }
    else if (argc == 4 ) {
        cout << "Good right amount of arguments" << endl;
        _inputname = "../data/metro/";
        _inputname += + argv[1];
        _guidename = "../data/guide/";
        _guidename += argv[2];
        _outputname += argv[3];
        _automaticCase = true;
    }
    else if (argc == 5 ) {
        cout << "Good right amount of arguments" << endl;
        _inputname = "../data/metro/";
        _inputname += + argv[1];
        _guidename = "../data/guide/";
        _guidename += argv[2];
        _outputname += argv[3];
        string argv_4 = argv[4];
        if (argv_4.compare("1") == 0)
            _automaticCase = true;
        else
            _automaticCase = false;
    }

    _outputname += "/output";
    _metro.load( _inputname );
    _guide.load(_guidename);


	_metro.adjustsize( ZuKai::Base::Common::getMainwidgetWidth(),
	                   ZuKai::Base::Common::getMainwidgetHeight() );
	_matching = Stationmatching(); 
	_matching.prepare(&_metro, &_guide, 
		ZuKai::Base::Common::getMainwidgetWidth()/2,
		ZuKai::Base::Common::getMainwidgetHeight()/2 );

    _autoMatching = AutoMatching();
    _autoMatching.init(&_metro, &_guide, ZuKai::Base::Common::getMainwidgetWidth()/2,
                       ZuKai::Base::Common::getMainwidgetHeight()/2);



	cout << "---------- Guide: ----------" << endl; 
	cout << "  Num Edges:		" << _guide.getNumEdges() << endl; 
	cout << "  Num Vertex:		" << _guide.getNumVertex() << endl; 
	cout << "----------------------------" << "\n" << endl; 
	
	cout << "---------- Metro: ----------" << endl; 
	cout << "  Num Edges:		" << _metro.getNumEdges() << endl; 
	cout << "  Num Vertex:		" << _metro.getNumStations() << endl; 
	cout << "----------------------------" << "\n" << endl;

	cout << "number of unplanar edges: " << _metro.nonPlanarIntersections() << endl;



    _manuelPath = Manuelpath();
    _manuelPath.init(&_metro, &_guide);

    if (_automaticCase) {
        _smooth._exclusivlyPathMode = false;
        run();
    } else {
        _smooth._exclusivlyPathMode = true;
    }

}

void MetroShapes::exportSteps( void ) {
    _stepsMixed = std::max(0., _stepsSmooth - _metro.getNumStations());
    _stepsSmooth -= _stepsMixed;
    cout << "Creating Layout till step " << _stepsSmooth << ", " << _stepsMixed << endl;

    // finding Path
    computeAutoMatching();
    distortMetroMap();

    // Creating Smooth Layout
    if (_stepsSmooth > 0) {
        _metro.removeAdditionalEdges();
        _smooth.prepare( &_metro, &_guide,
                         ZuKai::Base::Common::getMainwidgetWidth()/2,
                         ZuKai::Base::Common::getMainwidgetHeight()/2 );
        int iter = _stepsSmooth;
        _smooth.ConjugateGradient( iter );
        _smooth.retrieve();
        _smooth.clear();
    }

    // Creating Mixed Layout
    if (_stepsMixed > 0) {
        int iter = _stepsMixed;
        _mixedlayout.prepare( &_metro, &_guide,
                              ZuKai::Base::Common::getMainwidgetWidth()/2,
                              ZuKai::Base::Common::getMainwidgetHeight()/2 );
        _mixedlayout.ConjugateGradient( iter );
        _mixedlayout.retrieve();
        // _mixedlayout.preparePostProcessing();
        // _mixedlayout.postProcessing();
        _mixedlayout.clear();
    }
}

/*
*	Calculates Matching, Smooth and Mixed Layout.
*	prints time to console
*
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::run( void )
{
	clock_t overall_time, match_time, smooth_time, mixed_time = 0;
	overall_time = clock(); 
	match_time = clock();

    if (!_smooth._exclusivlyPathMode) {
        computeAutoMatching();
        distortMetroMap();
    }

	match_time = clock() - match_time;

    smooth_time = clock();
	computeSmooth();
	smooth_time = clock() - smooth_time;

    resetSmooth();

	mixed_time = clock();
	computeMixedlayout();
	mixed_time = clock()  - mixed_time;

    overall_time = clock() - overall_time;

	cout << "------------------------------------------" << endl; 
	cout << "  match time:		" << double(match_time) / CLOCKS_PER_SEC << endl; 
	cout << "  smooth time:		" << double(smooth_time) / CLOCKS_PER_SEC << endl; 
	cout << "  mixed time:		" << double(mixed_time) / CLOCKS_PER_SEC << endl; 
	cout << "  overall time:		" << double(overall_time) / CLOCKS_PER_SEC << endl; 
	cout << "  sum:			" << double(match_time) / CLOCKS_PER_SEC + double(smooth_time) / CLOCKS_PER_SEC + double(mixed_time) / CLOCKS_PER_SEC << endl; 
	cout << "------------------------------------------" << endl;

    // save files
    output();
}

/*
*	Transforms guide based on matching 	
*
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::computeMatching( void )
{
	if (_matching.getNumberMatchings() > 0 ) {
		_matching.init();
		_matching.conjugateGradient(); 
		_matching.retrieve(); 
	} else {
		cout << "Matching: Not enough connections" << endl; 
	}
}

void MetroShapes::computeAutoMatching() {
    cout << "caled compute Auto Matching" << endl;
    _metro.enhanceMetro(20); // Berlin Bear 30 // paris 50? –> large// paris old 35
    _autoMatching.findPath();
    _autoMatching.translateGuide();
}

void MetroShapes::manuallyAddVertexToPath( Coord2 coord ) {
    if (_manuelPath.addVertex(coord))
        _manuelPath.alingGuide();
}
void MetroShapes::addVertexToSelectPath( int id ) {
    cout << "adding Vertex with id: " << id << "to path" << endl;
}

void MetroShapes::selectPathFromMetroLine(int l ) {
    _autoMatching.pathFromMetroLine(l);
}

/*
*	Computes Smooth Layout 	
*
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::computeSmooth( void )
{
	_metro.removeAdditionalEdges();
	// _metro.subdivideComplexEdges();
	// _metro.reorderID();
    _smooth.prepare( &_metro, &_guide,
		ZuKai::Base::Common::getMainwidgetWidth()/2,
		ZuKai::Base::Common::getMainwidgetHeight()/2 );

    int iter = _metro.nStations();
	_smooth.ConjugateGradient( iter );
	_smooth.retrieve();
	_smooth.clear();
	cout << "number of unplanar edges: " << _metro.nonPlanarIntersections() << endl;
}

/*
*	Computes mixed layout 	
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::computeMixedlayout( void )
{
    _metro.reorderID();
	int iter = _metro.nStations();
	// _metro.removeUnecessaryTempVertex();
	_mixedlayout.prepare( &_metro, &_guide,
					  ZuKai::Base::Common::getMainwidgetWidth()/2,
					  ZuKai::Base::Common::getMainwidgetHeight()/2 );
	_mixedlayout.ConjugateGradient( iter );
	//
	_mixedlayout.retrieve();
	// _metro.subdiveToOctilinear();
	_metro.reInsertRemovedStations();
	// _mixedlayout.complexPostProcessing();
	_mixedlayout.clear();
	cout << "number of unplanar edges: " << _metro.nonPlanarIntersections() << endl;

}

/*
*	Distorts metro layout not uniformly based on matching  	
*
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::distortMetroMap( void ) 
{
	// _matching.distortMap();
    _autoMatching.translateGuide();
}

void MetroShapes::scaleNonUniformlyBasedOnManuelPath() {
    _manuelPath.scaleMetroNonUniformly();
}

/*
*	outputs resulting layout  	
*
*	Inputs	none
*	Outputs	none
*/
void MetroShapes::output( void ) 
{
    cout << "out the files" << endl;
    string outputname_metro_graphml = _outputname + "_metro" + ".graphml";
    string outputname_guide_graphml = _outputname + "_guide" + ".graphml";
    if( !_outputname.empty() ){
        _metro.exportGraphML( outputname_metro_graphml );
        _guide.exportGraphML(outputname_guide_graphml );
    }
}

void MetroShapes::resetSmooth( void ) 
{
	_smooth._resetSmoothPtr(); 
}

bool MetroShapes::addNode( Coord2 cord ) 
{
	return _matching.addNewNode(cord); 
}

bool MetroShapes::removeNode( Coord2 cord )
{
	return _matching.removeNode(cord); 
}

// end of header file
// Do not add any stuff under this line.
