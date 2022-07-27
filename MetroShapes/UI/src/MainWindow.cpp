//******************************************************************************
// MainWindow.cpp
//	: program file for main window
//
//------------------------------------------------------------------------------
//
//	Ver 2.00		Date: Sun Mar 14 20:00:00 2021
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <QtWidgets>
#include "MainWindow.h"

namespace Ui {


	
	//------------------------------------------------------------------------------
    //	Private functions
    //------------------------------------------------------------------------------
    //
    //  MainWindow::_init -- initialization
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_init( void )
    {
		 createMainView();

        _createActions();
        _createStatusBar();
        _createDockWindows();
		
		setWindowTitle( tr("MetroShapes") );
		setMouseTracking( false );
		setUnifiedTitleAndToolBarOnMac(true );
    }

    //
    //  MainWindow::init -- initialization
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::init( Metro *__m_ptr, Base *__b_ptr )
    {
    	_metroPtr = __m_ptr;
    	_basePtr = __b_ptr;
	    _init();
    }

    void MainWindow::init( Metro * __m_ptr, Base * __b_ptr, Smooth *__smooth_ptr, Mixedlayout *__octi_ptr, Guide * __g_ptr, Stationmatching * __matchingPtr, AutoMatching * __autoMatchingPtr)
    {
        _metroPtr = __m_ptr; 
        _basePtr = __b_ptr; 
        _guidePtr = __g_ptr; 
        _smoothPtr = __smooth_ptr; 
        _mixedPtr = __octi_ptr; 
        _matchingPtr = __matchingPtr;
        _autoMatchPtr = __autoMatchingPtr;
        _init();
    }

    //
    //  MainWindow::_createActions -- create actions
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_createActions()
    {
        QMenu *fileMenu = menuBar()->addMenu(tr("&File"));

        // tool box setting
        QToolBar *fileToolBar = addToolBar(tr("File"));
        fileToolBar->setIconSize(QSize(25, 25));
        //fileToolBar->setFixedHeight(25);
	
	    //------------------------------------------------------------------------------
	    //	tool box
	    //------------------------------------------------------------------------------
	    const QIcon newIcon = QIcon::fromTheme("", QIcon(":/icons/app.png"));
	    const QIcon newCircleIcon = QIcon::fromTheme("", QIcon(":/icons/circle.png"));
	    const QIcon newPolyIcon = QIcon::fromTheme("", QIcon(":/icons/poly.png"));
        QAction *newFileAct = new QAction(newIcon, tr("&Open"), this );
        newFileAct->setShortcuts( QKeySequence::New) ;
        newFileAct->setStatusTip(tr("Open a new file"));
        connect( newFileAct, &QAction::triggered, this, &MainWindow::newFile );
        fileMenu->addAction( newFileAct );
        fileToolBar->addAction( newFileAct );

        fileToolBar->addSeparator();
        QAction *autoMatch = new QAction(newIcon, tr("Auto Match"), this);
        fileToolBar->addAction(autoMatch);
        connect(autoMatch, &QAction::triggered, this, &MainWindow::autoMatch);

        QAction *smoothAct = new QAction(newCircleIcon, tr("QSmooth"), this);
        fileToolBar->addAction(smoothAct);
        connect( smoothAct, &QAction::triggered, this, &MainWindow::smooth );

        QAction *octiAct = new QAction(newPolyIcon, tr("Mixed Layout"), this);
        fileToolBar->addAction(octiAct);
        connect( octiAct, &QAction::triggered, this, &MainWindow::octi);






	    //------------------------------------------------------------------------------
	    //	file menu
	    //------------------------------------------------------------------------------

        fileMenu->addSeparator();

        QAction *quitAct = fileMenu->addAction( tr("&Quit"), this, &QWidget::close );
        quitAct->setShortcuts(QKeySequence::Quit);
        quitAct->setStatusTip(tr("Quit the application"));

        QAction *exportPngAct = fileMenu->addAction( tr("&Save PNG"), this , SLOT(savePNG()), Qt::CTRL + Qt::Key_P );
        exportPngAct->setStatusTip(tr("Saving PNG")); 

        QAction *exportSvgAct = fileMenu->addAction( tr("&Save SVG"), this, SLOT(saveSVG()), Qt::CTRL + Qt::Key_S );
        exportSvgAct->setStatusTip(tr("Saving SVG"));  

        QAction *outData = fileMenu->addAction(tr("&output Data"), this, SLOT(outputData()), Qt::CTRL + Qt::Key_O );
        outData->setStatusTip(tr("&output Data"));

        QMenu *editMenu = menuBar()->addMenu(tr("&Edit"));
        QToolBar *editToolBar = addToolBar(tr("Edit"));
        editToolBar->setIconSize(QSize(25, 25));
        editToolBar->setFixedHeight(50);
        
        menuBar()->addSeparator();

        _viewMenu = menuBar()->addMenu(tr("&View"));
        menuBar()->addSeparator();

        QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));

        QAction *aboutAct = helpMenu->addAction(tr("&About"), this, &MainWindow::about );
        aboutAct->setStatusTip(tr("Show the application's About box"));

        QAction *aboutQtAct = helpMenu->addAction(tr("About &Qt"), qApp, &QApplication::aboutQt );
        aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
    }

    //
    //  MainWindow::_createStatusBar -- create status bar
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_createStatusBar()
    {
        statusBar()->showMessage(tr("Ready"));
    }

    //
    //  MainWindow::_createDockWindows -- create dock windows
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_createDockWindows( void )
    {
        // ---------------------------------------------------------------------------
        // Metro Lines
        // ---------------------------------------------------------------------------
        _linesDock = new QDockWidget(tr("MetroLines"), this);
        _linesDock->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

        _lines = new QWidget( _linesDock );
        _lines->setGeometry( QRect(0,0,ZuKai::Base::Common::getDockWidgetWidth(),
									 ZuKai::Base::Common::getMainwidgetHeight()/2.0 ) );
        _lines->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
										 ZuKai::Base::Common::getMainwidgetHeight()/3.0 ) );
        _linesDock->setWidget( _lines );
        addDockWidget(Qt::LeftDockWidgetArea, _linesDock );
        _viewMenu->addAction(_linesDock->toggleViewAction() );
        _linesDock->hide();
        vector<bool> shapeLines = _metroPtr->getShapeLines();

        _linesLayout = new QGridLayout;

        QSignalMapper *signalMapper = new QSignalMapper(this);
        for(int i = 0; i != shapeLines.size(); i++) {
            QCheckBox *checkbox = new QCheckBox("",this);
            QPushButton *selBut = new QPushButton("&Select", this);
            connect(selBut, SIGNAL(clicked()), signalMapper, SLOT(map()));
            signalMapper->setMapping(selBut, i);

            if (shapeLines[i]) checkbox->setChecked(true);

            QLabel* lab = new QLabel(QLatin1String(" "));
            double r = _metroPtr->lineColor(i)[0] * 255;
            QString bgc{"background-color: rgb( "};
            bgc += QString::number(r) + ", ";
            bgc += QString::number(_metroPtr->lineColor(i)[1] * 255 ) + ", ";
            bgc += QString::number(_metroPtr->lineColor(i)[2]* 255) + ")";
            lab->setStyleSheet(bgc);

            _linesLayout->addWidget(checkbox, i,0, 1, 1);
            _linesLayout->addWidget(lab, i,6, 1 , 10);
            _linesLayout->addWidget(selBut, i, 1, 1, 5 );
            _boxesLines.push_back(checkbox);
        }
        connect(signalMapper, SIGNAL(mapped(int)),
            this, SLOT(selectLine(int)));
        _lines->setLayout(_linesLayout);

        // ---------------------------------------------------------------------------
        // Smooth varibles / Slider
        // ---------------------------------------------------------------------------
        _settingsDock = new QDockWidget(tr("Smooth Settings"), this );
        _settingsDock->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

	    _setting = new QWidget( _settingsDock );
        _setting->setGeometry( QRect(0,0,ZuKai::Base::Common::getDockWidgetWidth(),
									 ZuKai::Base::Common::getMainwidgetHeight()/2.0 ) );
        _setting->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
										 ZuKai::Base::Common::getMainwidgetHeight()/3.0 ) );
        _settingsDock->setWidget( _setting );

        addDockWidget(Qt::RightDockWidgetArea, _settingsDock );
        _viewMenu->addAction(_settingsDock->toggleViewAction() );

        _sliderEdge = new QDoubleSpinBox;
        _sliderEdge->setValue(_smoothPtr->_w_alongEdge);
        _sliderEdge->setSingleStep(0.2);
        QLabel* labelEdge = new QLabel(QLatin1String("w_c : Shape Approximation"));

        // _spinEdgeDist = new QDoubleSpinBox;
        // _spinEdgeDist->setValue(_smoothPtr->_constSmooth);
        // _spinEdgeDist->setSingleStep(0.01);
        // QLabel* labelEdgeDist = new QLabel(QLatin1String("Edge Dist"));

        _spinWangle = new QDoubleSpinBox;
        _spinWangle->setValue(_smoothPtr->_w_angle);
        _spinWangle->setSingleStep(0.01);
        QLabel* labelangle = new QLabel(QLatin1String("w_a : Max. Angle Weight"));

        _spinWposition = new QDoubleSpinBox;
        _spinWposition->setValue(_smoothPtr->_w_position);
        _spinWposition->setSingleStep(0.01);
        QLabel* labelposition = new QLabel(QLatin1String("w_p : Position Weight"));


        _spinLength = new QDoubleSpinBox; 
        _spinLength->setValue(_smoothPtr->_w_contextlength); 
        QLabel* labellength = new QLabel(QLatin1String("w_l : Uniform Edge Lenght Weight"));


        _spinWcrossing = new QDoubleSpinBox;
        _spinWcrossing->setValue(_smoothPtr->_w_crossing);
        _spinWcrossing->setSingleStep(0.01);
        QLabel* labelCrossing = new QLabel(QLatin1Literal("Crossing Weight"));

        // _spinGamma = new QDoubleSpinBox;
        // _spinGamma->setValue(_smoothPtr->_gamma);
        // _spinGamma->setSingleStep(0.01);
        // QLabel* labelGamma = new QLabel(QLatin1Literal("Gamma"));

        // _spinWoverlap = new QDoubleSpinBox;
        // _spinWoverlap->setValue(_smoothPtr->_w_overlapping);
        // _spinWoverlap->setSingleStep(0.01);
        // QLabel* labelWoverlap = new QLabel(QLatin1String("Overlapping weight"));

       _wLayout = new QGridLayout;
       _wLayout->addWidget(labelEdge, 1,0);
       _wLayout->addWidget(_sliderEdge,1,1);
       // _wLayout->addWidget(labelEdgeDist, 2, 0);
       // _wLayout->addWidget(_spinEdgeDist, 2, 1);
       _wLayout->addWidget(labellength, 2, 0);
       _wLayout->addWidget(_spinLength, 2, 1); 
       _wLayout->addWidget(labelangle, 3, 0);
       _wLayout->addWidget(_spinWangle, 3,1);
       _wLayout->addWidget(labelposition, 4, 0);
       _wLayout->addWidget(_spinWposition, 4, 1);
       _wLayout->addWidget(labelCrossing, 5, 0);
       _wLayout->addWidget(_spinWcrossing, 5, 1);
       // _wLayout->addWidget(labelGamma, 6, 0);
       // _wLayout->addWidget(_spinGamma, 6, 1);
       // _wLayout->addWidget(labelWoverlap, 6, 0);
       // _wLayout->addWidget(_spinWoverlap, 6, 1);

       _setting->setLayout(_wLayout);

        // ---------------------------------------------------------------------------
        // Mixedlayout variables / Slider
        // ---------------------------------------------------------------------------
        _octiDock = new QDockWidget(tr("Mixed Settings"), this);
        _octiDock->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        _octi = new QWidget(_octiDock);
        _octi->setGeometry( QRect(0,0,ZuKai::Base::Common::getDockWidgetWidth(),
										 ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );
        _octi->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
											 ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );

        _spinOcti = new QDoubleSpinBox;
        _spinOcti->setValue(_mixedPtr->_w_mixed);
        _spinOcti->setSingleStep(0.2);
        QLabel* labelOcti = new QLabel(QLatin1String("w_o : Octilinear Weight"));

        _spinOctiClose = new QDoubleSpinBox; 
        _spinOctiClose->setValue(_mixedPtr->_w_position_path); 
        _spinOctiClose->setSingleStep(0.2);
        QLabel* labelOctiClose = new QLabel(QLatin1String("w_c : Shape Approximation"));

        _spinOctiPos = new QDoubleSpinBox;
        _spinOctiPos->setDecimals(5);
        _spinOctiPos->setValue(_mixedPtr->_w_position);
        QLabel* labelOctiPos = new QLabel(QLatin1String("w_p : Position Weight"));

        _spinOctiCros = new QDoubleSpinBox;
        _spinOctiCros->setValue(_mixedPtr->_w_crossing);
        QLabel* labelOctiCros = new QLabel(QLatin1String("Crossing Weight"));

        // _spinOctiGama = new QDoubleSpinBox;
        // _spinOctiGama->setValue(_mixedPtr->_gama);
        // QLabel* labelOctiGama = new QLabel(QLatin1String("min Dist nodes"));

        // _spinOctiOver = new QDoubleSpinBox;
        // _spinOctiOver->setValue(_mixedPtr->_w_overlap);
        // QLabel* labelOctiOver = new QLabel(QLatin1String("Overlay weight"));

        // _spinMixedConst = new QDoubleSpinBox;
        // _spinMixedConst->setValue(_mixedPtr->_constMixed);
        // QLabel* labelMixedConst = new QLabel(QLatin1String("C_m"));

        // _spinMixedMinAngle = new QDoubleSpinBox;
        // _spinMixedMinAngle->setValue(_mixedPtr->_minAngle);
        // QLabel* labelMixedAngle = new QLabel(QLatin1String("beta"));

        _mixedLayout = new QGridLayout;
        _mixedLayout->addWidget(labelOcti, 0,0);
        _mixedLayout->addWidget(_spinOcti, 0, 1);
        _mixedLayout->addWidget(labelOctiPos, 1, 0);
        _mixedLayout->addWidget(_spinOctiPos, 1, 1);
        _mixedLayout->addWidget(labelOctiClose, 2, 0); 
        _mixedLayout->addWidget(_spinOctiClose, 2, 1); 
        _mixedLayout->addWidget(labelOctiCros, 3,0);
        _mixedLayout->addWidget(_spinOctiCros, 3, 1);

        // _mixedLayout->addWidget(labelOctiGama, 3, 0);
        // _mixedLayout->addWidget(_spinOctiGama, 3, 1);
        // _mixedLayout->addWidget(labelOctiOver, 4,0);
        // _mixedLayout->addWidget(_spinOctiOver, 4, 1);
        // _mixedLayout->addWidget(labelMixedConst, 5, 0);
        // _mixedLayout->addWidget(_spinMixedConst, 5, 1);
        // _mixedLayout->addWidget(labelMixedAngle, 6, 0);
        // _mixedLayout->addWidget(_spinMixedMinAngle, 6, 1);
        _octi->setLayout(_mixedLayout);

        _octiDock->setWidget( _octi );
        addDockWidget(Qt::RightDockWidgetArea, _octiDock );
        _viewMenu->addAction( _octiDock->toggleViewAction() );

        // ---------------------------------------------------------------------------
        // Matching Doc
        // ---------------------------------------------------------------------------
        _matchDock = new QDockWidget(tr("Route Matching"), this);
        _matchDock->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        _match = new QWidget(_matchDock);
        _match->setGeometry( QRect(0,0,ZuKai::Base::Common::getDockWidgetWidth(),
                                         ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );
        _match->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
                                             ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );

        _wMetro = new QDoubleSpinBox;
        _wMetro->setValue(_autoMatchPtr->getCostMetro());
        QLabel* labelWMetro = new QLabel(QLatin1String("cost Metro"));

        _wAdd = new QDoubleSpinBox;
        _wAdd->setValue(_autoMatchPtr->getCostAdd());
        QLabel* labelWAdd = new QLabel(QLatin1String("cost Add Edges"));

        _wColor = new QDoubleSpinBox;
        _wColor->setValue(_autoMatchPtr->getWeightColor());
        QLabel* labelWColor = new QLabel(QLatin1String("weigth Color"));

        _machtLayout = new QGridLayout;
        _machtLayout->addWidget(labelWMetro, 0,0);
        _machtLayout->addWidget(_wMetro, 0, 1);
        _machtLayout->addWidget(labelWAdd, 1, 0);
        _machtLayout->addWidget(_wAdd, 1, 1);
        _machtLayout->addWidget(labelWColor, 3, 0);
        _machtLayout->addWidget(_wColor, 3, 1);

        _match->setLayout(_machtLayout);
        _matchDock->setWidget(_match);

        addDockWidget(Qt::RightDockWidgetArea, _matchDock );
        _viewMenu->addAction( _matchDock->toggleViewAction() );

        // ---------------------------------------------------------------------------
        // Interaction Doc
        // ---------------------------------------------------------------------------
	    _interactionDock = new QDockWidget(tr("Interaction"), this );
	    _interactionDock->setAllowedAreas( Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
        _interaction = new QWidget( _interactionDock );
        _interaction->setGeometry( QRect(0,0,ZuKai::Base::Common::getDockWidgetWidth(),
										 ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );
        _interaction->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
											 ZuKai::Base::Common::getDockWidgetWidth()*ZuKai::Base::Common::getMainwidgetHeight()/ZuKai::Base::Common::getMainwidgetWidth() ) );


        _sliderOpacity = new QSlider(Qt::Horizontal);
        _sliderOpacity->setValue(_mainGV->_guideOpacity);
        _sliderOpacity->setRange(0, 255);
        QLabel* labelOpacity = new QLabel(QLatin1String("Guide Opacity"));
        connect(_sliderOpacity, SIGNAL(valueChanged(int)), this, SLOT(reDraw()));

        _sliderScale = new QSlider(Qt::Horizontal);
        _lastScale = 100;
        _sliderScale->setValue(_lastScale);
        _sliderScale->setRange(0.3, 200.0);
        QLabel* labelScale = new QLabel(QLatin1String("Translate Mode: Guide Scale"));
        connect(_sliderScale, SIGNAL(valueChanged(int)), this, SLOT(scaleGuide()));

        connect(_sliderOpacity, SIGNAL(valueChanged(int)), this, SLOT(reDraw()));

        _buttonKeyPoints = new QRadioButton("Key Points", this);
        _buttonKeyPoints->setChecked((_mainGV->_selectionMode == 0));
        connect(_buttonKeyPoints, SIGNAL(clicked()), this, SLOT(resetInteractionMode()));

        _buttonManuelPath = new QRadioButton("Manuel Path", this);
        _buttonManuelPath->setChecked((_mainGV->_selectionMode == 1));
        connect(_buttonManuelPath, SIGNAL(clicked()), this, SLOT(resetInteractionMode()));

        _buttonTranslatePath =  new QRadioButton("Translate Path", this);
        _buttonTranslatePath->setChecked((_mainGV->_selectionMode == 2));
        connect(_buttonTranslatePath, SIGNAL(clicked()), this, SLOT(resetInteractionMode()));

        _buttonReCalc = new QPushButton("&Calculate Based On User Path (Manual Case)", this);
        connect(_buttonReCalc, SIGNAL(clicked()), this, SLOT(reCalcManually()));
        // _buttonSmooth = new QPushButton("Smooth", this);
        // connect(_buttonSmooth, SIGNAL(clicked()), this, SLOT(smooth()));

        _buttonCalcManually  = new QPushButton("&Calculate Automatic Case", this);
        connect(_buttonCalcManually, SIGNAL(clicked()), this, SLOT(reCalc()));


        // _buttonNonUniformlyManuel = new QPushButton("Non Unifoly \n man. Path", this);
        // connect(_buttonNonUniformlyManuel, SIGNAL(clicked()), this, SLOT(scaleNonUniformlyBasedOnManuelPath()));

        _buttonMatch = new QPushButton("Keypoints Mode: Align", this);
        // _buttonDistort = new QPushButton("Distord based on Path", this);
        // _buttonOcti = new QPushButton("Mixedlayout", this);
        connect(_buttonMatch, SIGNAL(clicked()), this, SLOT(match()));
        // connect(_buttonDistort, SIGNAL(clicked()), this, SLOT(distort()));
        // connect(_buttonOcti, SIGNAL(clicked()), this, SLOT(octi()));

        QGridLayout *_interLayout = new QGridLayout;
        _interLayout->addWidget(labelOpacity, 0,0);
        _interLayout->addWidget(_sliderOpacity, 0, 1);

        // _interLayout->addWidget(_buttonNonUniformlyManuel, 1, 0);
        _interLayout->addWidget(_buttonCalcManually, 2, 0);
        _interLayout->addWidget(_buttonReCalc, 3,0);
        _interLayout->addWidget(_buttonMatch, 4,0);
        // _interLayout->addWidget(_buttonDistort, 4,1);
        // _interLayout->addWidget(_buttonSmooth, 5,0);
        // _interLayout->addWidget(_buttonOcti, 5,1);
        // _interLayout->addWidget(_boxHighlightMatch, 6, 0);

        _interLayout->addWidget(labelScale, 5, 0);
        _interLayout->addWidget(_sliderScale, 5, 1);
        _interLayout->addWidget(_buttonKeyPoints, 6, 0);
        _interLayout->addWidget(_buttonManuelPath, 7, 0);
        _interLayout->addWidget(_buttonTranslatePath, 8, 0);


        _interaction->setLayout(_interLayout);

	    _interactionDock->setWidget( _interaction );
        addDockWidget(Qt::LeftDockWidgetArea, _interactionDock );
        _viewMenu->addAction( _interactionDock->toggleViewAction() );

    }

    //
    //  MainWindow::_updateVaribles -- updates variables from spins
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_updateVaribles( void ) 
    {
        // Metro Lines Update 
        vector<bool> shapeLines = _metroPtr->getShapeLines(); 
        for(vector<int>::size_type i = 0; i != shapeLines.size(); i++) {
            QCheckBox* box = _boxesLines.value(i);
            _metroPtr->setShapeLine(i, box->isChecked());
        }
        _metroPtr->updateEffectedMetroLines();

        // _mainGV->highlightMatchPath = _boxHighlightMatch->isChecked();
        // cerr << _boxHighlightMatch->isChecked() << endl;

        // Setting Varibles
        _smoothPtr->_w_alongEdge = _sliderEdge->value(); 
        // _smoothPtr->_constSmooth = _spinEdgeDist->value();
        _smoothPtr->_w_angle = _spinWangle->value(); 
        _smoothPtr->_w_contextlength = _spinLength->value();
        _smoothPtr->_w_position = _spinWposition->value(); 
        _smoothPtr->_w_crossing = _spinWcrossing->value();  

        // Octi Varibles 
        _mixedPtr->_w_mixed = _spinOcti->value();  
        _mixedPtr->_w_position = _spinOctiPos->value(); 
        _mixedPtr->_w_crossing = _spinOctiCros->value(); 
        // _mixedPtr->_gama = _spinOctiGama->value(); 
        _mixedPtr->_w_position_path = _spinOctiClose->value();
        // _mixedPtr->_w_overlap = _spinOctiOver->value(); 
        // _mixedPtr->_constMixed = _spinMixedConst->value(); 
        // _mixedPtr->_minAngle = _spinMixedMinAngle->value();

        // Matching
        _autoMatchPtr->setCostAdd(_wAdd->value());
        _autoMatchPtr->setCostMetro(_wMetro->value());
        _autoMatchPtr->setWeightColor(_wColor->value());
        // _autoMatchPtr->setTolerancePartInprovment(_tolPartInp->value());
        reDraw();
    }

    //
    //  MainWindow::_reCalc -- recalculate Metro Map
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::reCalc() {
        _updateVaribles();

        cout << "  Calculating a path (Automatic Case)" << endl;
        _smoothPtr->_exclusivlyPathMode = false;

        _basePtr->run();
        _mainGV->_highlightPath = false;
        _mainGV->_highlightSmoothPath = true;
        reDraw();
        cout << "... Recalcution Done" << endl; 
    }

    void MainWindow::reCalcManually() {
        scaleNonUniformlyBasedOnManuelPath();
        _updateVaribles();

        cout << "  Taking user Path (Manual Case)" << endl;
        _smoothPtr->_exclusivlyPathMode = true;

        _basePtr->run();
        _mainGV->_highlightPath = false;
        _mainGV->_highlightSmoothPath = true;
        reDraw();
        cout << "... Recalcution Done" << endl;
    }



    void MainWindow::resetInteractionMode() {
        cout << "resetInteractionMode" << endl;
        int i = 0;
        if (_buttonKeyPoints->isChecked())
            i = 0;
        else if (_buttonManuelPath->isChecked())
            i = 1;
        else if (_buttonTranslatePath->isChecked())
            i = 2;
        _mainGV->_selectionMode = i;
        _buttonKeyPoints->setChecked((_mainGV->_selectionMode == 0));
        _buttonManuelPath->setChecked((_mainGV->_selectionMode == 1));
        _buttonTranslatePath->setChecked((_mainGV->_selectionMode == 2));
    }

    void MainWindow::scaleGuide() {
        if (_mainGV->_selectionMode == 2) {
            double delta = _sliderScale->value() / 100.0 - _lastScale / 100.0 + 1.0;
            if (delta  != 1.0 )
                _guidePtr->scale(delta);
            _lastScale = _sliderScale->value();
            reDraw();
        }
    }

    //
    //  MainWindow::match -- calculate Matching
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::match() {
        if (_mainGV->_selectionMode == 0) {
            _basePtr->computeMatching();
            // _basePtr->alignKeyPointsWithoutRotation();
        } else if (_mainGV->_selectionMode == 1 ) {
            _basePtr->manuallyAddVertexToPath(Coord2(10000.0, 100000.0));
        }
        reDraw();
    }

    void MainWindow::autoMatch() {
        _updateVaribles();
        _basePtr->computeAutoMatching();
        reDraw();
    }

    void MainWindow::selectLine(int l) {
        cerr << l << endl;
        _basePtr->selectPathFromMetroLine(l);
        reDraw();
    }

    //
    //  MainWindow::distort -- distort Map not uniformly, based on matching
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::distort() 
    {
        _basePtr->distortMetroMap();
        reDraw();
    }

    void MainWindow::scaleNonUniformlyBasedOnManuelPath() {
        _basePtr->scaleNonUniformlyBasedOnManuelPath();
        reDraw();
    }

    //
    //  MainWindow::smooth -- calculate Smooth Layout
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::smooth() {
        _updateVaribles(); 
        _basePtr->computeSmooth(); 
        // _smoothPtr->_resetSmoothPtr(); 
        _basePtr->resetSmooth();
        _mainGV->_highlightPath = false;
        _mainGV->_highlightSmoothPath = true;
        reDraw(); 
    }

    void MainWindow::octi() {
        _updateVaribles();
        _basePtr->computeMixedlayout();
        _mainGV->_highlightPath = false;
        _mainGV->_highlightSmoothPath = true;
        reDraw();
    }

    //
    //  MainWindow::reDraw -- draw Graphics view
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::reDraw() {
        _mainGV->_guideOpacity = _sliderOpacity->value();
        // _mainGV->highlightMatchPath = _boxHighlightMatch->isChecked();
        _mainGV[0]._item_metro();
    }
	
    void MainWindow::savePNG() {
        _mainGV->exportPNG(3840, 2160 ); 
    }

    void MainWindow::savePNG(double w, double h, string path) {
        _mainGV->exportPNG(w, h, path);
    }

    //
    //  MainWindow::saveSVG -- export as S
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::saveSVG() {
        _mainGV->exportSVG(1200, 1200 ); 
    }

    void MainWindow::saveSVG(double w, double h, string path) {
        _mainGV->exportSVG(w, h, path);
    }

    void MainWindow::outputData() {
        _basePtr->output();
    }




    //------------------------------------------------------------------------------
    //	Protected functions
    //------------------------------------------------------------------------------
    //
    //  MainWindow::updateSetting -- update setting
    //
    //  Inputs
    //  setting: text string
    //
    //  Outputs
    //  none
    //
    void MainWindow::_updateSetting( const QString &setting )
    {
        // _setting->updateSceneItems();
        if( setting.isEmpty() ) return;
    }

    //
    //  MainWindow::updateInteraction -- update interaction
    //
    //  Inputs
    //  interaction: text string
    //
    //  Outputs
    //  none
    //
    void MainWindow::_updateInteraction( const QString &interaction )
    {
//        _interaction->updateSceneItems();
        if( interaction.isEmpty() ) return;
    }

    //
    //  MainWindow::_updateAllDocks -- update all dock widgets
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::_updateAllDocks( void )
    {
        _updateSetting( "" );
        _updateInteraction( "" );
    }

    //------------------------------------------------------------------------------
    //	Public functions
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //	Constructors & Destructors
    //------------------------------------------------------------------------------
    //
    //  MainWindow::MainWindow -- constructor
    //
    //  Inputs
    //  parent: parent widget
    //
    //  Outputs
    //  none
    //
    MainWindow::MainWindow( QWidget *parent )
        : QMainWindow( parent )
    {
	    //------------------------------------------------------------------------------
	    // initialization
	    //------------------------------------------------------------------------------
	    _metroPtr = nullptr;
	    _basePtr = nullptr;
	
	    //------------------------------------------------------------------------------
        // clear stored images
        //------------------------------------------------------------------------------
        QString path = "../../svg/";
        QDir dir( path );
        dir.setNameFilters( QStringList() << "*.*" );
        dir.setFilter( QDir::Files );
        for( const QString &dirFile: dir.entryList() ) {
            dir.remove( dirFile );
        }

        //------------------------------------------------------------------------------
        // configuration file
        //------------------------------------------------------------------------------
       
        string configFilePath = qApp->applicationDirPath().toStdString() + "/../config/MainWindow.conf";
        ZuKai::Base::Config conf( configFilePath );

        int icon_width = 0,
            icon_height = 0;

        if ( conf.has( "mainwidget_width" ) ){
            string parammainwidgetWidth = conf.gets( "mainwidget_width" );
            ZuKai::Base::Common::setMainwidgetWidth( stoi( parammainwidgetWidth ) );
        }
        if ( conf.has( "mainwidget_height" ) ){
            string parammainwidgetHeight = conf.gets( "mainwidget_height" );
	        ZuKai::Base::Common::setMainwidgetHeight( stoi( parammainwidgetHeight ) );
        }
        if ( conf.has( "dockwidget_width" ) ){
            string paramDockWidgetWidth = conf.gets( "dockwidget_width" );
	        ZuKai::Base::Common::setDockWidgetWidth( stoi( paramDockWidgetWidth ) );
       }
        if ( conf.has( "menubar_height" ) ){
            string paramMenuBarHeight = conf.gets( "menubar_height" );
	        ZuKai::Base::Common::setMenubarHeight( stoi( paramMenuBarHeight ) );
        }
        if ( conf.has( "icon_width" ) ){
            string paramIconWidth = conf.gets( "icon_width" );
            icon_width = stoi( paramIconWidth );
        }
        if ( conf.has( "icon_height" ) ){
            string paramIconHeight = conf.gets( "icon_height" );
            icon_height = stoi( paramIconHeight );
        }

        setGeometry( QRect( 50, 50,
		 					ZuKai::Base::Common::getMainwidgetWidth()+ZuKai::Base::Common::getDockWidgetWidth(),
		 					ZuKai::Base::Common::getMainwidgetHeight() + ZuKai::Base::Common::getMenubarHeight() ) );
	    setMinimumSize(QSize( ZuKai::Base::Common::getMainwidgetWidth()+ZuKai::Base::Common::getDockWidgetWidth(),
						   ZuKai::Base::Common::getMainwidgetHeight() + ZuKai::Base::Common::getMenubarHeight() ) );
    }


    void MainWindow::createMainView(void)
    {
	    _mainGV = new Vector::GraphicsView( this );
        _mainGV->setStyleSheet("background: white; border: transparent;");
	    _mainGV->setGeometry( QRect( 0, 0, ZuKai::Base::Common::getDockWidgetWidth(),
	                                 ZuKai::Base::Common::getMainwidgetHeight() ) );
	    _mainGV->setMinimumSize( QSize( ZuKai::Base::Common::getDockWidgetWidth(),
	                                    ZuKai::Base::Common::getMainwidgetHeight() ) );
	    _mainGV->setMouseTracking( true );
	    _mainGV->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
	    _mainGV->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        if (_guidePtr != nullptr) {   //no pointer has been asigned
            _mainGV->init( _metroPtr, _basePtr, _guidePtr, _matchingPtr); 
        } else {
	        _mainGV->init( _metroPtr, _basePtr );
        }
	    _mainGV->initSceneItems();

	    setCentralWidget( _mainGV );
    }
    

    //------------------------------------------------------------------------------
    //	Event handlers
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Special functions
    //------------------------------------------------------------------------------
    //
    //  MainWindow::newFile -- open new file
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::newFile( void )
    {
        _mainGV->simulateKey( Qt::Key_L );
    }

    //
    //  MainWindow::print -- print
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::print( void )
    {
    }

    //
    //  MainWindow::undo -- undo
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //Base
    void MainWindow::undo( void )
    {
    }

    //
    //  MainWindow::save -- save
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::save( void )
    {
    #ifdef REVISE
        QMimeDatabase mimeDatabase;
        QString fileName = QFileDialog::getSaveFileName(this,
                            tr("Choose a file name"), ".",
                            mimeDatabase.mimeTypeForName("text/html").filterString());
        if (fileName.isEmpty())
            return;
        QFile file(fileName);
        if (!file.open(QFile::WriteOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("Dock Widgets"),
                                 tr("Cannot write file %1:\n%2.")
                                 .arg(QDir::toNativeSeparators(fileName), file.errorString()));
            return;
        }

        QTextStream out(&file);
        QApplication::setOverrideCursor(Qt::WaitCursor);
        out << textEdit->toHtml();
        QApplication::restoreOverrideCursor();

        statusBar()->showMessage(tr("Saved '%1'").arg(fileName), 2000);
    #endif // REVISE
    }

    //
    //  MainWindow::about -- about
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void MainWindow::about( void )
    {
       QMessageBox::about(this, tr("About KeiRo"),
                tr( "Map-based biological pathway diagram"));
    }
	
} // namespace Ui