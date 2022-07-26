#ifndef WIDGET_H
#define WIDGET_H

#include <QtOpenGL/QGLWidget>
#include <QtCore/QBasicTimer>
#include <QtCore/QThread>
#include <GLUT/glut.h>
#include <unistd.h>

// OpenCV
#include <opencv/cv.h>
#include <opencv/highgui.h>

// For preventing conflict with Qt and Boost library
#ifndef Q_MOC_RUN
#include "Common.h"
#include "Metro.h"
//#include "RNGraph.h"
//#include "Smooth.h"
//#include "Mixedlayout.h"
#endif


class Widget : public QGLWidget
{
    Q_OBJECT
private:
    Metro       *_metro;
    Metro       *_simmetro;

    QBasicTimer *       _timer;
    int         _time_step;
    METROTYPE   _time_type;

    // mouse
    int         _last_pointer_x,    // Mouse positions with last pressed
                _last_pointer_y;
    int         _now_pointer_x,     // Mouse positions 
                _now_pointer_y;

    bool        _isDeformed;
    
    // focus point 
    VertexDescriptor focusVD;
    // start and target point of the shortest path
    VertexDescriptor sourceVD;
    VertexDescriptor targetVD;

    // ui
    static bool _sim_flag;  // o = original, 1 = simplified

protected:
    void initializeGL( void );
    void paintGL( void );
    void _paintDisk( const double & radius );
    void _paintCircle( const double & radius );
    void _paintStation( VertexDescriptor vd, METROTYPE type );
    void _paintMetro( METROTYPE type );
    void _paintName( void );
    void _paintLabels( void );
    void _paintMagnification( void );

    void _plotStations( void );
    void _plotLines( void );
    bool _handleVertex( int nHits, unsigned int * buf, const int button, const int modifier );
    bool _handleEdge( int nHits, unsigned int * buf );

    void timerEvent      ( QTimerEvent * e );

public:
    Widget( const QGLFormat& format, QWidget *parent = 0 );
    ~Widget();

    void init( Metro * __metro, Metro * __simmetro );
    void initTextures( const char * filename );
    void capture( char * name );

    // I/O
    void setMouseMovingPointer( int __last_pointer_x, int __last_pointer_y ) {
        _now_pointer_x = __last_pointer_x;
        _now_pointer_y = __last_pointer_y;
        //cerr << "point.x = " << _now_pointer_x << endl;
        //cerr << "point.y = " << _now_pointer_y << endl;
    }
    void setMousePointer( int __last_pointer_x, int __last_pointer_y ) {
        _last_pointer_x = __last_pointer_x;
        _last_pointer_y = __last_pointer_y;
    }

    // UI
    void selectMagnification( double x, double y, int button, int modifier );
    void selectTopology( void );
    void setMetroWeight( void );
    void setSimMetroWeight( void );
    void pickVertex( double x, double y, int button, int modifier );
    void pickEdge( double x, double y, int button, int modifier );
    bool calcShortestPath( void );
    bool calcGeodesicDistance( void );
    void clearPath( void ) {
        sourceVD = NULL;
        targetVD = NULL;
    }
    void animationSmooth( void );
    void animationOctilinear( void );
    bool & isDeformed( void )               { return _isDeformed; }

    // static variable access
    static bool & simFlag( void )           { return _sim_flag; }

    // timer
    const QBasicTimer * timer( void ) const { return _timer; }
    void start( METROTYPE type )                      {
        _time_step = 0;
        _time_type = type;
        _isDeformed = false;
        _timer->start( 30, this );
    }
    void stop( void )                       {
        _timer->stop();
    }

signals:

public slots:

};

#endif // WIDGET_H
