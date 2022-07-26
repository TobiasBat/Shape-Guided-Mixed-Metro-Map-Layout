//******************************************************************************
// GraphicsVertexItem.h
//	: header file for vertex items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

#ifndef Ui_Vector_GraphicsVertexItem_H
#define Ui_Vector_GraphicsVertexItem_H

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
#include <cassert>

using namespace std;

//#ifndef Q_MOC_RUN
//#include "Base/Coord2.h"
//#include "Base/Common.h"
//#endif // Q_MOC_RUN

#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsItem>
#include <QtGui/QPainter>
#include <QtCore/QString>
#include <QtWidgets/QGraphicsSceneMouseEvent>

namespace Ui {
namespace Vector {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class GraphicsVertexItem : public  QGraphicsRectItem
    {
    private:

        unsigned int    _id;
        QString         _name;
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
        explicit GraphicsVertexItem( QGraphicsItem *parent = Q_NULLPTR );
        // parameterized constructor
        explicit GraphicsVertexItem( qreal x, qreal y, qreal w, qreal h, QGraphicsItem *parent = Q_NULLPTR );
        // copy constructor
        explicit GraphicsVertexItem( const QRectF &rect, QGraphicsItem *parent = Q_NULLPTR );
        // destructor
        ~GraphicsVertexItem( void ) {}

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

        QString &	            name( void )        { return _name; }
        const QString &	        name( void ) const	{ return _name; }

        int &	                fontSize( void ) 	    { return _font_size; }
        const int &	            fontSize( void ) const	{ return _font_size; }

        bool &	                textOn( void ) 	    { return _textOn; }
        const bool &	        textOn( void ) const{ return _textOn; }

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void    init      ( void );

    private:

    };

} // namespace Vector
} // namespace Ui

#endif // Ui_Vector_GraphicsVertexItem_H
