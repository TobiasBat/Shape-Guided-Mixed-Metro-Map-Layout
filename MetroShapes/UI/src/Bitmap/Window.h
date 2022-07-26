#ifndef WINDOW_H
#define WINDOW_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <list>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

#include <QtWidgets/QAction>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtGui/QMouseEvent>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QApplication>
#include <QtOpenGL/QGLFormat>

#include "Widget.h"
#include "Base.h"

#define REMOVEBACKNUM   (15)

class Window : public QMainWindow
{
    Q_OBJECT
private:
    Widget      *_widget;
    
    Metro       *_metro;
    Metro       *_simmetro;
//    RNGraph     *_rng;
//    Smooth      *_smooth;
//    Mixedlayout  *_octilinear;

    // display
    int         _content_width;
    int         _content_height;

    // mouse 
    int         _now_pointer_x,    // Mouse positions and its button status
                _now_pointer_y;
    int         _last_pointer_x,    // Last mouse positions after pressing
                _last_pointer_y;

    // OpenGL format
    QGLFormat format;

    // menu
    // load
    QMenu *loadMenu;
    QAction *selBerlinAct;
    QAction *selBostonAct;
    QAction *selKashiwaAct;
    QAction *selLisbonAct;
    QAction *selLondonAct;
    QAction *selMontrealAct;
    QAction *selMunichAct;
    QAction *selNuernbergAct;
    QAction *selParisAct;
    QAction *selPragueAct;
    QAction *selSingaporeAct;
    QAction *selTaipeiAct;
    QAction *selTokyoAct;
    QAction *selViennaAct;
    // simplification
    QMenu *simMenu;
    QAction *selCloneGraphAct;
    QAction *selNeighborGraphAct;
    QAction *selMinDistanceAct;
    QAction *selMovebackSmoothAct;
    QAction *selMovebackOctilinearAct;
    // optimization
    QMenu *optMenu;
    QAction *selSmoothLSAct;
    QAction *selOctilinearLSAct;
    QAction *selSmoothSmallCGAct;
    QAction *selSmoothCGAct;
    QAction *selOctilinearSmallCGAct;
    QAction *selOctilinearCGAct;
    QAction *selAllInOneAct;

    void createActions( void );
    void createMenus( void );

    void postLoad( void );

public slots:
//private slots:
    // load
    void selectBerlin( void );
    void selectBoston( void );
    void selectKashiwa( void );
    void selectLisbon( void );
    void selectLondon( void );
    void selectMontreal( void );
    void selectMunich( void );
    void selectNuernberg( void );
    void selectParis( void );
    void selectPrague( void );
    void selectSingapore( void );
    void selectTaipei( void );
    void selectTokyo( void );
    void selectVienna( void );

    // simplification
    void selectCloneGraph( void );
    void selectNeighborGraph( void );
    void selectMinDistance( void );
    void selectMovebackSmooth( void );
    void selectMovebackOctilinear( void );

    // optimization
    void selectSmoothLS( void );
    void selectOctilinearLS( void );
    void selectSmoothSmallCG( void );
    void selectSmoothCG( void );
    void selectOctilinearSmallCG( void );
    void selectOctilinearCG( void );
    void selectSmooth( OPTTYPE opttype = CONJUGATE_GRADIENT );
    void selectOctilinear( OPTTYPE opttype = CONJUGATE_GRADIENT );
    void selectAllInOne( void );

public:
    explicit Window( QGLWidget *parent = 0 );
    ~Window();

    void initializeGL( void );
    void paintGL( void );
	void init( Base * __b_ptr );
//    void init( Metro * __metro, Metro * __simmetro, RNGraph * __rng,
//               Smooth * __smooth, Mixedlayout * __octilinear );

protected:
    void keyPressEvent( QKeyEvent *event );
    void mouseMoveEvent( QMouseEvent *event );
    void mousePressEvent( QMouseEvent *event );
    void mouseReleaseEvent( QMouseEvent *event );

};


#endif // WINDOW_H
