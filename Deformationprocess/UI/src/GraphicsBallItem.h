//******************************************************************************
// GraphicsBallItem.h
//	: header file for ball items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

#ifndef Ui_Vector_GraphicsBallItem_H
#define Ui_Vector_GraphicsBallItem_H

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

using namespace std;

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsItem>
#include <QtGui/QPainter>
#include <QtCore/QString>
#include <QtWidgets/QGraphicsSceneMouseEvent>

//#ifndef Q_MOC_RUN
//#include "base/Coord2.h"
//#include "base/Common.h"
//#endif // Q_MOC_RUN

namespace Ui {
namespace Vector {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class GraphicsBallItem : public  QGraphicsRectItem
    {
    private:

        unsigned int    _id;
        QString         _name;
        QString         _text;
        bool            _textOn;

        int             _font_size;
        QFont           _font;
        QPen            _textpen;
        double          _radius;

    protected:

    public:

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        explicit GraphicsBallItem( QGraphicsItem *parent = Q_NULLPTR );
        // parameterized constructor
        explicit GraphicsBallItem( qreal x, qreal y, qreal w, qreal h, QGraphicsItem *parent = Q_NULLPTR );
        // copy constructor
        explicit GraphicsBallItem( const QRectF &rect, QGraphicsItem *parent = Q_NULLPTR );
        // destructor
        ~GraphicsBallItem( void ) {}

        // source from the qt library
        //QRectF rect() const;
        //void setRect(const QRectF &rect);
        //inline void setRect(qreal x, qreal y, qreal w, qreal h);

        //QPainterPath shape() const Q_DECL_OVERRIDE;
        //bool contains(const QPointF &point) const Q_DECL_OVERRIDE;

        //bool isObscuredBy(const QGraphicsItem *item) const Q_DECL_OVERRIDE;
        //QPainterPath opaqueArea() const Q_DECL_OVERRIDE;

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

        int &	                fontSize( void ) 	    { return _font_size; }
        const int &	            fontSize( void ) const	{ return _font_size; }

        bool &	                textOn( void ) 	    { return _textOn; }
        const bool &	        textOn( void ) const{ return _textOn; }

        double &	            radius( void ) 	    { return _radius; }
        const double &	        radius( void ) const{ return _radius; }

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void    init      ( void );
        void setRadius( double d ); 

    };

} // namespace Vector
} // namespace Ui

#endif // Ui_Vector_GraphicsBallItem_H
