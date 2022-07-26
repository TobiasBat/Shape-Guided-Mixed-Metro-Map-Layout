#include "Window.h"

ofstream    ofs( "batch.txt" );

//----------------------------------------------------------
// Window
//----------------------------------------------------------
Window::Window( QGLWidget *parent )
    : QMainWindow( parent )
{
    format.setVersion( 3, 3 );
    format.setProfile( QGLFormat::CoreProfile ); // Requires >= Qt-4.8.0
    format.setSampleBuffers( true );

    _widget = new Widget( format, this );
    setCentralWidget( _widget );
    _widget->setMouseTracking( true );
    setMouseTracking( true );

    createActions();
    createMenus();
}

Window::~Window()
{
}

// void Window::init( Metro * __metro, Metro * __simmetro, RNGraph *__rng, Smooth * __smooth, Mixedlayout * __octilinear )
void Window::init( Base * __b_ptr )
{
//    _metro = __metro;
//    _simmetro = __simmetro;
//    _smooth = __smooth;
//    _octilinear = __octilinear;
//    _rng = __rng;
//
//    _widget->init( _metro, _simmetro );

    int lrmargin = LEFTRIGHT_MARGIN;
    int tdmargin = TOPBOTTOM_MARGIN;
    _content_width = width() - lrmargin;
    _content_height = height() - tdmargin;
}

void Window::createActions( void )
{
    // load
    selBerlinAct = new QAction( tr("B&erlin"), this );
    connect( selBerlinAct, SIGNAL( triggered() ), this, SLOT( selectBerlin() ) );
    selBostonAct = new QAction( tr("B&oston"), this );
    connect( selBostonAct, SIGNAL( triggered() ), this, SLOT( selectBoston() ) );
    selKashiwaAct = new QAction( tr("K&ashiwa"), this );
    connect( selKashiwaAct, SIGNAL( triggered() ), this, SLOT( selectKashiwa() ) );
    selLisbonAct = new QAction( tr("L&isbon"), this );
    connect( selLisbonAct, SIGNAL( triggered() ), this, SLOT( selectLisbon() ) );
    selLondonAct = new QAction( tr("L&ondon"), this );
    connect( selLondonAct, SIGNAL( triggered() ), this, SLOT( selectLondon() ) );
    selMontrealAct = new QAction( tr("M&ontreal"), this );
    connect( selMontrealAct, SIGNAL( triggered() ), this, SLOT( selectMontreal() ) );
    selMunichAct = new QAction( tr("M&unich"), this );
    connect( selMunichAct, SIGNAL( triggered() ), this, SLOT( selectMunich() ) );
    selNuernbergAct = new QAction( tr("N&uernberg"), this );
    connect( selNuernbergAct, SIGNAL( triggered() ), this, SLOT( selectNuernberg() ) );
    selParisAct = new QAction( tr("P&aris"), this );
    connect( selParisAct, SIGNAL( triggered() ), this, SLOT( selectParis() ) );
    selPragueAct = new QAction( tr("P&rague"), this );
    connect( selPragueAct, SIGNAL( triggered() ), this, SLOT( selectPrague() ) );
    selSingaporeAct = new QAction( tr("S&ingapore"), this );
    connect( selSingaporeAct, SIGNAL( triggered() ), this, SLOT( selectSingapore() ) );
    selTaipeiAct = new QAction( tr("T&aipei"), this );
    connect( selTaipeiAct, SIGNAL( triggered() ), this, SLOT( selectTaipei() ) );
    selTokyoAct = new QAction( tr("T&okyo"), this );
    connect( selTokyoAct, SIGNAL( triggered() ), this, SLOT( selectTokyo() ) );
    selViennaAct = new QAction( tr("V&ienna"), this );
    connect( selViennaAct, SIGNAL( triggered() ), this, SLOT( selectVienna() ) );

    // simplification
    selCloneGraphAct = new QAction( tr("C&loneGraph"), this );
    connect( selCloneGraphAct, SIGNAL( triggered() ), this, SLOT( selectCloneGraph() ) );
    selNeighborGraphAct = new QAction( tr("N&eighborGraph"), this );
    connect( selNeighborGraphAct, SIGNAL( triggered() ), this, SLOT( selectNeighborGraph() ) );
    selMinDistanceAct = new QAction( tr("M&inDistance"), this );
    connect( selMinDistanceAct, SIGNAL( triggered() ), this, SLOT( selectMinDistance() ) );
    selMovebackSmoothAct = new QAction( tr("M&ovebackSmooth"), this );
    connect( selMovebackSmoothAct, SIGNAL( triggered() ), this, SLOT( selectMovebackSmooth() ) );
    selMovebackOctilinearAct = new QAction( tr("M&ovebackOctilinear"), this );
    connect( selMovebackOctilinearAct, SIGNAL( triggered() ), this, SLOT( selectMovebackOctilinear() ) );

    // optimization
    // least square
    selSmoothLSAct = new QAction( tr("S&mooth Original (LS)"), this );
    connect( selSmoothLSAct, SIGNAL( triggered() ), this, SLOT( selectSmoothLS() ) );
    selOctilinearLSAct = new QAction( tr("O&ctilinear Original (LS)"), this );
    connect( selOctilinearLSAct, SIGNAL( triggered() ), this, SLOT( selectOctilinearLS() ) );
    // conjugate gradient
    selSmoothSmallCGAct = new QAction( tr("S&mooth Small (CG)"), this );
    connect( selSmoothSmallCGAct, SIGNAL( triggered() ), this, SLOT( selectSmoothSmallCG() ) );
    selSmoothCGAct = new QAction( tr("S&mooth Original (CG)"), this );
    connect( selSmoothCGAct, SIGNAL( triggered() ), this, SLOT( selectSmoothCG() ) );
    selOctilinearSmallCGAct = new QAction( tr("O&ctilinear Small (CG)"), this );
    connect( selOctilinearSmallCGAct, SIGNAL( triggered() ), this, SLOT( selectOctilinearSmallCG() ) );
    selOctilinearCGAct = new QAction( tr("O&ctilinear Original (CG)"), this );
    connect( selOctilinearCGAct, SIGNAL( triggered() ), this, SLOT( selectOctilinearCG() ) );
    selAllInOneAct = new QAction( tr("A&ll in One"), this );
    connect( selAllInOneAct, SIGNAL( triggered() ), this, SLOT( selectAllInOne() ) );
}

void Window::createMenus( void )
{
    // load
    loadMenu = menuBar()->addMenu( tr("&Load") );
    loadMenu->addAction( selBerlinAct );
    loadMenu->addAction( selBostonAct );
    loadMenu->addAction( selKashiwaAct );
    loadMenu->addAction( selLisbonAct );
    loadMenu->addAction( selLondonAct );
    loadMenu->addAction( selMontrealAct );
    loadMenu->addAction( selMunichAct );
    loadMenu->addAction( selNuernbergAct );
    loadMenu->addAction( selParisAct );
    loadMenu->addAction( selPragueAct );
    loadMenu->addAction( selSingaporeAct );
    loadMenu->addAction( selTaipeiAct );
    loadMenu->addAction( selTokyoAct );
    loadMenu->addAction( selViennaAct );

    // simplification
    simMenu = menuBar()->addMenu( tr("&Simplification") );
    simMenu->addAction( selCloneGraphAct );
    simMenu->addAction( selNeighborGraphAct );
    simMenu->addAction( selMinDistanceAct );
    simMenu->addAction( selMovebackSmoothAct );
    simMenu->addAction( selMovebackOctilinearAct );

    // optimization
    optMenu = menuBar()->addMenu( tr("&Optimization") );
    optMenu->addAction( selSmoothLSAct );
    optMenu->addAction( selOctilinearLSAct );
    optMenu->addAction( selSmoothSmallCGAct );
    optMenu->addAction( selSmoothCGAct );
    optMenu->addAction( selOctilinearSmallCGAct );
    optMenu->addAction( selOctilinearCGAct );
    optMenu->addAction( selAllInOneAct );
}

void Window::selectBerlin( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/berlin-ubahn.txt" );
    strcpy( labelfilename,    "image/Berlin/berlin-label.txt" );
    strcpy( imagefilename,    "image/Berlin/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectBoston( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/boston-subway.txt" );
    strcpy( labelfilename,    "image/Boston/boston-label.txt" );
    strcpy( imagefilename,    "image/Boston/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectKashiwa( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/kashiwa-railway.txt" );
    strcpy( labelfilename,    "image/Kashiwa/kashiwa-label.txt" );
    strcpy( imagefilename,    "image/Kashiwa/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectLisbon( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/lisbon-metro.txt" );
    strcpy( labelfilename,    "image/Lisbon/lisbon-label.txt" );
    strcpy( imagefilename,    "image/Lisbon/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectLondon( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/london-underground.txt" );
    strcpy( labelfilename,    "image/London/london-label.txt" );
    strcpy( imagefilename,    "image/London/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectMontreal( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/montreal-metro.txt" );
    strcpy( labelfilename,    "image/Montreal/montreal-label.txt" );
    strcpy( imagefilename,    "image/Montreal/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectMunich( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/munich-ubahn.txt" );
    strcpy( labelfilename,    "image/Munich/munich-label.txt" );
    strcpy( imagefilename,    "image/Munich/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectNuernberg( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/nuernberg-ubahn.txt" );
    strcpy( labelfilename,    "image/Neurnberg/nuernberg-label.txt" );
    strcpy( imagefilename,    "image/Neurnberg/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectParis( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/paris-metro.txt" );
    strcpy( labelfilename,    "image/Paris/paris-label.txt" );
    strcpy( imagefilename,    "image/Paris/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectPrague( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/prague-metro.txt" ); // <--
    //strcpy( datafilename, "metrodata/prague-label.txt" );
    //strcpy( datafilename, "metrodata/prague-simple.txt" );
    //strcpy( datafilename, "metrodata/test01.txt" );
    //strcpy( datafilename, "metrodata/test.txt" );   
    strcpy( labelfilename,    "image/Prague/prague-label.txt" );  
    strcpy( imagefilename,    "image/Prague/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectSingapore( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/singapore-MRT.txt" );
    strcpy( labelfilename,    "image/Singapore/singapore-label.txt" );
    strcpy( imagefilename,    "image/Singapore/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectTaipei( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/taipei-metro.txt" );
    strcpy( labelfilename,    "image/Taipei/taipei-label.txt" );
    //strcpy( datafilename, "metrodata/taipei-old.txt" );
    strcpy( imagefilename,    "image/Taipei/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectTokyo( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/tokyo-all.txt" );
    strcpy( labelfilename,    "image/Tokyo/tokyo-label.txt" );
    strcpy( imagefilename,    "image/Tokyo/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::selectVienna( void )
{
    // load file
    static char datafilename[ BUFFER_SIZE ];
    static char labelfilename[ BUFFER_SIZE ];
    static char imagefilename[ BUFFER_SIZE ];
    strcpy( datafilename, "metrodata/vienna-ubahn.txt" );
    strcpy( labelfilename,    "image/Vienna/vienna-label.txt" );
    strcpy( imagefilename,    "image/Vienna/" );
    _metro->load( datafilename );
    _metro->loadLabel( labelfilename );
    _widget->initTextures( imagefilename );

    postLoad();
}

void Window::postLoad( void )
{
    //cerr << " width()/2 = " << width()/2 << " height()/2 = " << height()/2 << endl;
    //_metro->adjustsize( width()/2, height()/2 );
    _metro->adjustsize( _content_width/2, _content_height/2 );
    _metro->clearConflicts();
    _widget->clearPath();
    _widget->repaint();
    update();
}

void Window::selectCloneGraph( void )
{
    _simmetro->cloneLayout( *_metro );
    _simmetro->cloneLabel();
    _simmetro->clearConflicts();
}

void Window::selectNeighborGraph( void )
{
//    _rng->init( _simmetro );
//    //_rng->constrainedTriangulation();
//    _rng->relativeNeighborhoodGraph( _last_pointer_x, _last_pointer_y );
//    _widget->setSimMetroWeight();
}

void Window::selectMinDistance( void )
{
    _simmetro->simplifyLayout();
    _widget->repaint();
    update();
}

void Window::selectMovebackSmooth( void )
{
    bool isFinished = true;
    isFinished = _simmetro->movebackSmooth( *_metro );
    // if( isFinished == true ) cerr << "All stations are moved back!!" << endl;
    _widget->repaint();
    update();
}

void Window::selectMovebackOctilinear( void )
{
    bool isFinished = true;
    isFinished = _simmetro->movebackOctilinear( *_metro );
    // if( isFinished == true ) cerr << "All stations are moved back!!" << endl;
    _widget->repaint();
    update();
}

void Window::selectSmoothSmallCG( void )
{    
//    // run coarse smooth optimization
//    clock_t start_time = clock();
//    double err = 0.0;
//    unsigned int nLabels = _simmetro->nLabels();
//    _smooth->prepare( _simmetro, _content_width, _content_height );
//    err = _smooth->ConjugateGradient( 5 * _simmetro->nStations() );
//    _smooth->retrieve();
//
//    // cerr << "simNStation = " << _simmetro->nStations() << endl;
//    // ofs << "    Coarse CG: " << clock() - start_time << " err = " << err << " iter = " << 2 * _simmetro->nStations() << endl;
//
//    // add sub labels
//    //_simmetro->cloneSubLabel();
//
//    // run smooth optimization
//    while( true ) {
//
//        clock_t start = clock();
//        int iter = 0;
//
//        // check if all nodes are moved back
//#ifdef  DEBUG
//        cerr << " num_vertices( _simmetro->g() ) = " << num_vertices( _simmetro->g() ) << endl
//             << " num_vertices( _metro->g() ) = " << num_vertices( _metro->g() ) << endl
//             << " nLabels = " << nLabels << endl;
//#endif  // DEBUG
//        if( num_vertices( _simmetro->g() ) == ( num_vertices( _metro->g() ) + nLabels ) ) {
//            break;
//        }
//        else {
//            for( int i = 0; i < REMOVEBACKNUM; i++ ){
//                selectMovebackSmooth();
//            }
//        }
//
//        _smooth->prepare( _simmetro, _content_width/2, _content_height/2 );
//        if( num_vertices( _simmetro->g() ) == ( num_vertices( _metro->g() ) + nLabels ) ) {
//            iter = MAX( 2 * ( _metro->nStations() + nLabels - _simmetro->nStations() ), 30 );
//        }
//        else{
//            iter = MAX( 2 * ( _metro->nStations() + nLabels - _simmetro->nStations() ), 30 );
//        }
//        err = _smooth->ConjugateGradient( iter );
//        _smooth->retrieve();
//
//        cerr << "    time each loop = " << clock() - start << " err = " << err << " iter = " << iter << endl;
//        // ofs << "    time each loop = " << clock() - start << " err = " << err << " iter = " << iter << endl;
//    }
//
//    _smooth->clear();
//    _metro->cloneSmooth( *_simmetro );
//
//    // render the result
//    // printGraph( g );
//    //_widget->repaint();
//    //update();
//
//    // cerr << "Total Time CG = " << clock() - start_time << endl;
//    // ofs << "Total Time CG = " << clock() - start_time << endl;
}

void Window::selectSmooth( OPTTYPE opttype )
{
//    // run smooth optimization
//    int iter = 0;
//    _smooth->prepare( _metro, _content_width/2, _content_height/2 );
//    switch( opttype ){
//        case LEAST_SQUARE:
//            iter = 2 * _metro->nStations();
//            //iter = 20;
//            _smooth->LeastSquare( iter );
//            break;
//        case CONJUGATE_GRADIENT:
//            iter = 2 * _metro->nStations();
//            _smooth->ConjugateGradient( iter );
//            break;
//        default:
//            break;
//    }
//    _smooth->retrieve();
//    _smooth->clear();
//
//    // printGraph( g );
//    _widget->repaint();
//    update();
}

void Window::selectSmoothLS( void )
{
    clock_t start_time = clock();
    selectSmooth( LEAST_SQUARE );
    cerr << "Time of selectSmoothLS = " << clock() - start_time << endl;
}

void Window::selectSmoothCG( void )
{
    clock_t start_time = clock();
    selectSmooth( CONJUGATE_GRADIENT );
    cerr << "Time of selectSmoothCG = " << clock() - start_time << endl;
}

void Window::selectOctilinearSmallCG( void )
{
//    // run coarse octilinear optimization
//    clock_t start_time = clock();
//    double err = 0.0;
//    unsigned int nLabels = _simmetro->nLabels();
//
//    _octilinear->prepare( _simmetro, _content_width/2, _content_height/2 );
//    err = _octilinear->ConjugateGradient( 5 * _simmetro->nStations() );
//    _octilinear->retrieve();
//    //ofs << "    Coarse CG: " << clock() - start_time << " err = " << err << " iter = " << 5 * _simmetro->nStations() << endl;
//
//    // run octilinear optimization
//    while( true ) {
//
//        clock_t start = clock();
//        int iter = 0;
//
//        // check if all nodes are moved back
//        if( num_vertices( _simmetro->g() ) == ( num_vertices( _metro->g() ) + nLabels ) ) {
//            break;
//        }
//        else {
//            for( int i = 0; i < REMOVEBACKNUM; i++ ){
//                selectMovebackOctilinear();
//            }
//        }
//
//        _octilinear->prepare( _simmetro, _content_width/2, _content_height/2 );
//        if( num_vertices( _simmetro->g() ) == ( num_vertices( _metro->g() ) + nLabels ) ) {
//            iter = MAX( 2 * ( _metro->nStations() + nLabels - _simmetro->nStations() ), 30 );
//        }
//        else{
//            iter = MAX( 2 * ( _metro->nStations() + nLabels - _simmetro->nStations() ), 30 );
//        }
//        err = _octilinear->ConjugateGradient( iter );
//        _octilinear->retrieve();
//
//        // ofs << "    time each loop = " << clock() - start << " err = " << err << " iter = " << iter << endl;
//    }
//
//    _octilinear->clear();
//    _metro->cloneOctilinear( *_simmetro );
//
//    // render the result
//    // printGraph( g );
//    //_widget->repaint();
//    //update();
//
//    // cerr << "Total Time CG = " << clock() - start_time << endl;
//    // ofs << "Total Time CG = " << clock() - start_time << endl;
}

void Window::selectOctilinear( OPTTYPE opttype )
{
//    // run octilinear optimization
//    double iter = 0;
//    _octilinear->prepare( _metro, _content_width/2, _content_height/2 );
//    switch( opttype ) {
//        case LEAST_SQUARE:
//            iter = 2 * _metro->nStations();
//            _octilinear->LeastSquare( iter );
//            break;
//        case CONJUGATE_GRADIENT:
//            iter = 2 * _metro->nStations();
//            _octilinear->ConjugateGradient( iter );
//            break;
//        default:
//            break;
//    }
//    _octilinear->retrieve();
//    _octilinear->clear();

//    // printGraph( g );
//    _widget->repaint();
//    update();
}

void Window::selectOctilinearLS( void )
{
    clock_t start_time = clock();
    selectOctilinear( LEAST_SQUARE );
    cerr << "Time of selectOctilinearLS = " << clock() - start_time << endl;
}

void Window::selectOctilinearCG( void )
{
    clock_t start_time = clock();
    selectOctilinear( CONJUGATE_GRADIENT );
    cerr << "Time of selectOctilinearCG = " << clock() - start_time << endl;
}

void Window::selectAllInOne( void )
{
    // Key_C
    selectCloneGraph();
    //cerr << "test1" << endl;

    // Key_N
    selectNeighborGraph();
    //cerr << "test2" << endl;

    // Key_M
    //selectMinDistance();
    //cerr << "test3" << endl;

    // Key_S
    selectSmoothSmallCG();
    //_widget->animationSmooth();
    //cerr << "test4" << endl;

    // Key_O
    selectCloneGraph();
    _widget->start( SMOOTH );
    //selectNeighborGraph();
    //selectMinDistance();
    //QCoreApplication::processEvents();
    //usleep( 1000000 );
    //while( _widget->timer()->isActive() );
    selectOctilinearSmallCG();
    _widget->start( OCTILINEAR );
    //_widget->animationOctilinear();
}

void Window::keyPressEvent( QKeyEvent *e )
{
    char buf[ BUFFER_SIZE ], filename[ BUFFER_SIZE ], *ptr;
    // l -> c -> n -> m -> s -> c -> m -> o
    switch( e->key() ){
//        case Qt::Key_Plus:
//            magnified_radius += 10;
//            cerr << "magnified_radius = " << magnified_radius << endl;
//            break;
//        case Qt::Key_Minus:
//            magnified_radius -= 10;
//            cerr << "magnified_radius = " << magnified_radius << endl;
//            break;
//        case Qt::Key_Less:
//            rng_samples--;
//            cerr << "rng_samples = " << rng_samples << endl;
//            break;
//        case Qt::Key_Greater:
//            rng_samples++;
//            cerr << "rng_samples = " << rng_samples << endl;
//            break;
//        case Qt::Key_BraceLeft:
//            if( doi_max > doi_min ) doi_max--;
//            cerr << "doi_max = " << doi_max << endl;
//            break;
//        case Qt::Key_BraceRight:
//            doi_max++;
//            cerr << "doi_max = " << doi_max << endl;
//            break;
//        case Qt::Key_Colon:
//            doi_min--;
//            cerr << "doi_min = " << doi_min << endl;
//            break;
//        case Qt::Key_QuoteDbl:
//            if( doi_min < doi_max ) doi_min++;
//            cerr << "doi_min = " << doi_min << endl;
//            break;
        case Qt::Key_I:
            {
            Graph & g = _metro->g();
            resetGraphCoord( g );
            resetGraphWeight( g );
            _widget->isDeformed() = false;
            }
            break;
        case Qt::Key_A:
            selectAllInOne();
            break;
        case Qt::Key_P:
            fprintf( stderr, " Input file name : " );
            ptr = fgets( buf, sizeof( buf ), stdin );
            sscanf( buf, "%s", filename );
            _widget->capture( filename );
            break;
        case Qt::Key_L:
            //selectPrague();
            //selectTaipei();
            selectVienna();
            break;
        case Qt::Key_C:
            selectCloneGraph();
            break;
        case Qt::Key_N:
            selectNeighborGraph();
            break;
        case Qt::Key_M:
            selectMinDistance();
            break;
        case Qt::Key_S:
            _metro->updateTempCoord();
            selectSmoothSmallCG();
            _widget->start( SMOOTH );
            //_widget->simFlag() = true;
            break;
        case Qt::Key_O:
            selectOctilinearSmallCG();
            _widget->start( OCTILINEAR );
            break;
        case Qt::Key_B:
            selectMovebackOctilinear();
            break;
        case Qt::Key_F:
            _widget->simFlag() = !_widget->simFlag();
            if ( _widget->simFlag() ) cerr << "Draw simplified metro layout" << endl;
            else cerr << "Draw original metro layout" << endl;
            break;
        case Qt::Key_Escape:
            close();
            break;
        default:
            // cerr << "key = " << hex << e->key() << endl;
            QWidget::keyPressEvent( e );
    }
    _widget->repaint();
    update();
}

void Window::mouseMoveEvent( QMouseEvent *event )
{
    _now_pointer_x = event->x() - width()/2;
    _now_pointer_y = -event->y() + height()/2;
    _widget->setMouseMovingPointer( _now_pointer_x, _now_pointer_y );
    switch( event->buttons() ){
    case Qt::LeftButton:
        if ( event->type() == QEvent::MouseMove ){
#ifdef  SKIP
            if( event->modifiers() == Qt::ShiftModifier ) {
                _widget->pickVertexSet( event->x(), event->y(), Qt::LeftButton, Qt::ShiftModifier );
            }
#endif  // SKIP
        }
        break;
    case Qt::RightButton:
        if ( event->type() == QEvent::MouseMove ){
            _now_pointer_x = event->x() - width()/2;
            _now_pointer_y = -event->y() + height()/2;
        }
        break;
    case Qt::MiddleButton:
        if ( event->type() == QEvent::MouseMove ){
            _now_pointer_x = event->x() - width()/2;
            _now_pointer_y = -event->y() + height()/2;
        }
        break;
    default:
        QWidget::mouseMoveEvent( event );
    }
    update();
}

void Window::mousePressEvent( QMouseEvent *event )
{
    bool isShortest = false;

    switch( event->buttons() ){
    case Qt::LeftButton:
        if ( event->type() == QEvent::MouseButtonPress ){
            _now_pointer_x = event->x() - width()/2;
            _now_pointer_y = -event->y() + height()/2;
            _last_pointer_x = _now_pointer_x;
            _last_pointer_y = _now_pointer_y;
            // cerr << "eventX = " << event->x() << " point_x = " << _now_pointer_x << endl;
            if( event->modifiers() == Qt::ShiftModifier ) {
                _widget->setMousePointer( _now_pointer_x, _now_pointer_y );
                _widget->selectMagnification( event->x(), event->y(), Qt::LeftButton, Qt::ShiftModifier );
                _widget->selectTopology();
                _widget->setMetroWeight();
                //selectAllInOne();
                //_widget->calcGeodesicDistance();
                //_widget->pickVertex( event->x(), event->y(), Qt::LeftButton, Qt::ShiftModifier );
                //_widget->pickEdge( event->x(), event->y(), Qt::LeftButton, Qt::ShiftModifier );
                //isShortest = _widget->calcGeodesicDistance();
            }
            if( event->modifiers() == Qt::ControlModifier ) {   // "command" button of Mac
                _widget->pickVertex( event->x(), event->y(), Qt::LeftButton, Qt::ControlModifier );
                isShortest = _widget->calcShortestPath();
            }
        }
        break;
    case Qt::RightButton:
        if ( event->type() == QEvent::MouseButtonPress ){
            _now_pointer_x = event->x() - width()/2;
            _now_pointer_y = -event->y() + height()/2;
        }
        break;
    case Qt::MiddleButton:
        if ( event->type() == QEvent::MouseButtonPress ){
            //_now_pointer_x = event->x() - width()/2;
            //_now_pointer_y = -event->y() + height()/2;
            if( event->modifiers() == Qt::ControlModifier ) {   // "command" button of Mac
                _widget->pickVertex( event->x(), event->y(), Qt::MiddleButton, Qt::ControlModifier );
                isShortest = _widget->calcShortestPath();
            }
        }
        break;
    default:
        QWidget::mousePressEvent( event );
    }

    //if( isShortest == true ) selectAllInOne();

    update();
}

void Window::mouseReleaseEvent( QMouseEvent *event )
{
    switch( event->buttons() ){
    case Qt::LeftButton:
        if ( event->type() == QEvent::MouseButtonRelease ){
        }
        break;
    case Qt::RightButton:
        if ( event->type() == QEvent::MouseButtonRelease ){
        }
        break;
    case Qt::MiddleButton:
        if ( event->type() == QEvent::MouseButtonPress ){
        }
        break;
    default:
        QWidget::mouseReleaseEvent( event );
    }
    update();
}


