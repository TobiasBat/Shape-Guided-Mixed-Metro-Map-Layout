//******************************************************************************
// GraphicsPolygonItem.h
//	: header file for polygon items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

#ifndef Ui_Vector_GraphicsPolygonItem_H
#define Ui_Vector_GraphicsPolygonItem_H

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <string>

using namespace std;

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsItem>
#include <QtGui/QPainter>
#include <QtCore/QString>
#include <QtWidgets/QGraphicsSceneMouseEvent>

//#ifndef Q_MOC_RUN
//#include "base/Common.h"
//#include "base/Config.h"
//#endif // Q_MOC_RUN

namespace Ui {
namespace Vector {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class GraphicsPolygonItem : public QGraphicsPolygonItem
    {
    private:

        unsigned int    _id;
        QString         _text;
        bool            _textOn;

        int             _font_size;
        QFont           _font;
        QRect           _bbox;

    protected:

    public:

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        explicit GraphicsPolygonItem( QGraphicsItem *parent = Q_NULLPTR );
        // copy constructor
        explicit GraphicsPolygonItem( QPolygonF &polygon, QGraphicsItem *parent = Q_NULLPTR ) {}
        // destructor
        ~GraphicsPolygonItem( void ) {}

        // source from the qt library
        //QPainterPath path() const;
        //void setPath(const QPainterPath &path);

        //QPainterPath shape() const Q_DECL_OVERRIDE;
        //bool contains(const QPointF &point) const Q_DECL_OVERRIDE;

        //------------------------------------------------------------------------------
        //      Reimplementation
        //------------------------------------------------------------------------------
        QRectF boundingRect() const Q_DECL_OVERRIDE;
        int type( void ) const Q_DECL_OVERRIDE;

        void paint( QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget = Q_NULLPTR ) Q_DECL_OVERRIDE;

        //------------------------------------------------------------------------------
        //      Reference to elements
        //------------------------------------------------------------------------------

        unsigned int &	        id( void ) 	        { return _id; }
        const unsigned int &	id( void ) const	{ return _id; }

        QString &	            text( void )        { return _text; }
        const QString &	        text( void ) const	{ return _text; }

        bool &	                textOn( void ) 	    { return _textOn; }
        const bool &	        textOn( void ) const{ return _textOn; }

        QRect &	                bbox( void )        { return _bbox; }
        const QRect &	        bbox( void ) const	{ return _bbox; }

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------

    };

} // namespace Vector
} // namespace Ui

#endif // Ui_Vector_GraphicsPolygonItem_H
