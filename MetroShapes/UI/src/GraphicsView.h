//******************************************************************************
// GraphicsView.h
//	: header file for graphics scene
//
//------------------------------------------------------------------------------
//
//	Ver 2.00		Date: Sun Mar 14 20:00:00 2021
//
//******************************************************************************

#ifndef Ui_GraphicsView_H
#define Ui_GraphicsView_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <thread>

using namespace std;

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QApplication>
#include <QtWidgets/QGraphicsSceneDragDropEvent>
#include <QtWidgets/QFileDialog>
#include <QtGui/QMouseEvent>
#include <QtCore/QMimeData>
#include <QtCore/QDir>
#include <QtCore/QTimer>
#include <QtSvg/QSvgGenerator>

#ifndef Q_MOC_RUN
#include "Base.h"
#include "Metro.h"
#include "Guide.h"
#include "../../MetroShapes/matching/Stationmatching.h"
#include "GraphicsVertexItem.h"
#include "GraphicsBallItem.h"
#include "GraphicsEdgeItem.h"
#endif // Q_MOC_RUN

#define DOTTEDLINES
// #define BLACKWHITE

// #define HIGHLIGHTSMOOTHPATH //final style remove for debug
#define HIGHLIGHTMATCHPATH

//------------------------------------------------------------------------------
//	Defining Macros
//------------------------------------------------------------------------------

namespace Ui{
namespace Vector {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class GraphicsView : public QGraphicsView
    {
        Q_OBJECT

    private:

        // picking testing
//        QPainterPath                    _selectionArea;

        //------------------------------------------------------------------------------
        //	Data
        //------------------------------------------------------------------------------
        // scene
        QGraphicsScene                  *_scenePtr;
	    Metro                           *_metroPtr;
	    Base                            *_basePtr;
        Guide                            *_guidePtr;
        Stationmatching                 *_matchingPtr;

        bool posterStyle = true;

        //------------------------------------------------------------------------------
        //	UI
        //------------------------------------------------------------------------------
        double      _min_point_distance;

        QPoint      _oldCursor, _cursor;
        bool        _left_button_pressed,
                    _middle_button_pressed,
                    _right_button_pressed;
        QPointF lastTrackedPointerPosition;

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void _init              ( void );
        void _clear             ( void );
        void _calcLineWidth     ( double nStations ); 

    protected:

        //------------------------------------------------------------------------------
        //	Event handlers
        //------------------------------------------------------------------------------
        void keyPressEvent      ( QKeyEvent *event )    Q_DECL_OVERRIDE;
        void mousePressEvent    ( QMouseEvent *event )  Q_DECL_OVERRIDE;
        void mouseMoveEvent     ( QMouseEvent *event )  Q_DECL_OVERRIDE;
        void mouseReleaseEvent  ( QMouseEvent *event )  Q_DECL_OVERRIDE;
       

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------

    public:
        //------------------------------------------------------------------------------
        //Drawing Varibles 
        //------------------------------------------------------------------------------
        double                          _guideOpacity = 0;
        double                          _stationDot = 2; 
        double                          _lineWidth = 0.5;
        double                          _lineWidthHighlight = 1.0;
        bool                            _highlightPath = true;
        bool                            _highlightSmoothPath = false;
        bool                            highlightMatchPath = true;
        int  _selectionMode = 1; // 0 = keyPoints, 1 = manually
        //------------------------------------------------------------------------------
        //	Drawing functions
        //------------------------------------------------------------------------------
        void _item_metro                ( void );
        void _update_item_metro         ( void );

        //settings
        // QWidget* getSettings() override;

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        explicit GraphicsView( QWidget *parent = Q_NULLPTR );
        // copy constructor
        explicit GraphicsView( GraphicsView *parent = Q_NULLPTR ) {}
        // destructor
        ~GraphicsView( void ) {}

        //------------------------------------------------------------------------------
        //      Reimplementation
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //      Reference to elements
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //	    Special functions
        //------------------------------------------------------------------------------
        void init               ( Metro *__m_ptr, Base *__b_ptr ) {
	        _metroPtr = __m_ptr;
	        _basePtr = __b_ptr;
            _guidePtr = nullptr;
            _init();
        }
        
        void init               ( Metro *__m_ptr, Base *__b_ptr ,Guide * __g_ptr, Stationmatching * __matchingPtr ) {     //ADD smooth
	        _metroPtr = __m_ptr;
	        _basePtr = __b_ptr;
            _guidePtr = __g_ptr;
            _matchingPtr = __matchingPtr;
             _calcLineWidth(_metroPtr->getNumStations()); 
            _init();
        }
        
        void simulateKey        ( Qt::Key key );
        void initSceneItems     ( void );
        void updateSceneItems     ( void );

        void exportPNG ( double w, double h );
        void exportPNG ( double w, double h, string path);
        void exportSVG ( double w, double h );

        //------------------------------------------------------------------------------
        //	I/O functions
        //------------------------------------------------------------------------------
        // output
        friend ostream &	operator << ( ostream & s, const GraphicsView & m );
        // input
        friend istream &	operator >> ( istream & s, GraphicsView & m );
        // class name
        virtual const char * className( void ) const { return "GraphicsView"; }

    Q_SIGNALS:

    public Q_SLOTS:

    };

} // namespace Vector
} // namespace Ui

#endif // Ui_Vector_GraphicsView_H
