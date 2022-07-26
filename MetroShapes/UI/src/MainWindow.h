//******************************************************************************
// MainWindow.h
//	: header file for main window
//
//------------------------------------------------------------------------------
//
//	Ver 2.00		Date: Sun Mar 14 20:00:00 2021
//
//******************************************************************************

#ifndef Ui_MainWindow_H
#define Ui_MainWindow_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------
#include "istream"

#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMainWindow>
#include <QSlider>
#include <QDoubleSpinBox>
#include <QtWidgets>
#include <QGridLayout>

QT_BEGIN_NAMESPACE
//class QAction;
//class QMenu;
QT_END_NAMESPACE

#ifndef Q_MOC_RUN
#include "GraphicsView.h"
#include "Config.h"
#include "../../MetroShapes/optimization/Smooth.h"
#include "../../MetroShapes/optimization/Mixedlayout.h"
#include "../../MetroShapes/matching/Stationmatching.h"
#include "../../MetroShapes/matching/AutoMatching.h"
#endif // Q_MOC_RUN

//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace Ui {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class MainWindow : public QMainWindow
    {
        Q_OBJECT

    private:
    	
    	Metro                                       *_metroPtr;
        Base                                        *_basePtr;
        Guide                                        *_guidePtr;
        Stationmatching                             *_matchingPtr;
        AutoMatching                                *_autoMatchPtr;
        Smooth                                      *_smoothPtr;
        Mixedlayout                                  *_mixedPtr;

        Vector::GraphicsView                        *_mainGV;
	    QDockWidget                                 *_settingsDock;
	    QDockWidget                                 *_interactionDock;
        QDockWidget                                 *_octiDock; 
        QDockWidget                                 *_linesDock;
        QDockWidget                                 *_matchDock;
        QWidget                                     *_match;
        QWidget                                     *_lines;  
        QWidget                                     *_setting;
	    QWidget                                     *_interaction;
        QWidget                                     *_octi; 
        QMenu                                       *_viewMenu;
        

        QPushButton                                 *_buttonReCalc; 
        QPushButton                                 *_buttonMatch; 
        QPushButton                                 *_buttonSmooth; 
        QPushButton                                 *_buttonOcti; 
        QPushButton                                 *_buttonDistort;
        QPushButton                                 *_buttonNonUniformlyManuel;
        QGridLayout                                  *_wLayout;    
        QGridLayout                                 *_mixedLayout;
        QGridLayout                                 *_machtLayout;

        // -----------------------------------------------------------------------------
        // Shape Lines 
        // -----------------------------------------------------------------------------
        //vector<QCheckBox>                           *_boxesLines; 
        QList<QCheckBox*>                            _boxesLines; 
        QGridLayout                                 *_linesLayout; 

        // -----------------------------------------------------------------------------
        //  Smooth Weights 
        // -----------------------------------------------------------------------------
        QDoubleSpinBox                              *_spinWangle; 
        QDoubleSpinBox                              *_spinWposition; 
        QDoubleSpinBox                              *_spinWboundary; 
        QDoubleSpinBox                              *_spinWcrossing; 
        QDoubleSpinBox                              *_sliderParallel; 
        QDoubleSpinBox                              *_sliderEdge;
        QDoubleSpinBox                              *_spinEdgeDist;
        QDoubleSpinBox                              *_spinLength; 
        QDoubleSpinBox                              *_spinGamma; 
        QDoubleSpinBox                              *_spinWoverlap; 

        // -----------------------------------------------------------------------------
        //  Mixedlayout Weigts 
        // -----------------------------------------------------------------------------
        QDoubleSpinBox                              *_spinOcti;  
        QDoubleSpinBox                              *_spinOctiThress;
        QDoubleSpinBox                              *_spinOctiPos; 
        QDoubleSpinBox                              *_spinOctiCros; 
        QDoubleSpinBox                              *_spinOctiBound;  
        QDoubleSpinBox                              *_spinOctiGama; 
        QDoubleSpinBox                              *_spinOctiOver; 
        QDoubleSpinBox                              *_spinMixedConst; 
        QDoubleSpinBox                              *_spinMixedMinAngle;

        // -----------------------------------------------------------------------------
        //  Matching Weigts
        // -----------------------------------------------------------------------------
        QDoubleSpinBox                              *_wMetro;
        QDoubleSpinBox                              *_wAdd;
        QDoubleSpinBox                              *_tolPartInp;
        QDoubleSpinBox                              *_wColor;

        // ---------------------------------------------------------------------------
        // Interaction Doc
        // ---------------------------------------------------------------------------
        QSlider                                     *_sliderOpacity;   
        QSlider                                     *_sliderDot; 
        QSlider                                     *_sliderLineWidth;
        QCheckBox                                   *_boxHighlightMatch;
        QRadioButton                                *_buttonKeyPoints;
        QRadioButton                                *_buttonManuelPath;
        QRadioButton                                *_buttonTranslatePath;
        QSlider                                     *_sliderScale;
        double                                      _lastScale;

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void _init  ( void );
        void _createActions     ( void );
        void _createStatusBar   ( void );
        void _createDockWindows ( void );
        void _updateVaribles    ( void );



    public slots: 
        void reCalc( void ); 
        void reDraw( void ); 
        void match( void );
        void autoMatch( void );
        void distort( void );
        void scaleNonUniformlyBasedOnManuelPath();
        void smooth( void );
        void octi( void ); 
        void savePNG( void );
        void savePNG(double w, double h, string path);
        void saveSVG( void );
        void saveSVG(double w, double h, string path);
        void outputData( void );
        void selectLine( int i );
        void resetInteractionMode();
        void scaleGuide();
    signals:
        void clicked(int i);

    protected:

        void mouseMoveEvent     ( QMouseEvent *event )  Q_DECL_OVERRIDE {
        };


    public:

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        explicit MainWindow( QWidget *parent = Q_NULLPTR );
        // destructor
        ~MainWindow( void ) {}

        //------------------------------------------------------------------------------
        //	Specific functions
        //------------------------------------------------------------------------------
        void init  ( Metro *__m_ptr, Base *__b_ptr );

        void init ( Metro * __m_ptr, Base * __b_ptr, Smooth *__smooth_ptr, Mixedlayout *__octi_ptr, Guide * __g_ptr, Stationmatching * __matchingPtr, AutoMatching * __autoMatchingPtr);
	    //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const MainWindow & m );
        // input
        friend istream &	operator >> ( istream & s, MainWindow & m );
        // class name
        virtual const char * className( void ) const { return "MainWindow"; }

    private slots:

        void newFile( void );
//        void loadFileToTab( void );
        void save   ( void );
        void print  ( void );
        void undo   ( void );
        void about  ( void );
        void createMainView(void);
//        void createTabBasicLayoutView(void);
//        void createTabSbgnViewer(void);
        void _updateSetting      ( const QString &setting );
        void _updateInteraction  ( const QString &interaction );
        void _updateAllDocks    ( void );
//        void onTabChanged(int index);
//        void onTabCloseRequest(int index);
       
       
    };

} // namespace Ui

#endif // Ui_MainWindow_H
