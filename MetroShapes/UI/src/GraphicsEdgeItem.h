//******************************************************************************
// GraphicsEdgeItem.h
//	: header file for edge items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

#ifndef Ui_Vector_GraphicsEdgeItem_H
#define Ui_Vector_GraphicsEdgeItem_H

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

#ifndef Q_MOC_RUN
//#include "Base/Common.h"
#endif // Q_MOC_RUN

namespace Ui {
namespace Vector {

    //------------------------------------------------------------------------------
    //	Class definition
    //------------------------------------------------------------------------------
    class GraphicsEdgeItem : public  QGraphicsPathItem
    {
    private:

        unsigned int    _id;
        double          _weight;
        QString         _text;
        bool            _textOn;

    protected:

    public:

        //------------------------------------------------------------------------------
        //	Constructors & Destructors
        //------------------------------------------------------------------------------
        // default constructor
        explicit GraphicsEdgeItem( QGraphicsItem *parent = Q_NULLPTR );
        // copy constructor
        explicit GraphicsEdgeItem( const QPainterPath &path, QGraphicsItem *parent = Q_NULLPTR );
        // destructor
        ~GraphicsEdgeItem( void ) {}

        //------------------------------------------------------------------------------
        //      Reimplementation
        //------------------------------------------------------------------------------
        // source from the qt library
        //QPainterPath path() const;
        //void setPath(const QPainterPath &path);

        //QPainterPath shape() const Q_DECL_OVERRIDE;
        //bool contains(const QPointF &point) const Q_DECL_OVERRIDE;

        QRectF boundingRect() const Q_DECL_OVERRIDE;
        int type( void ) const Q_DECL_OVERRIDE;

        void paint( QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget = Q_NULLPTR ) Q_DECL_OVERRIDE;

        //------------------------------------------------------------------------------
        //      Reference to elements
        //------------------------------------------------------------------------------

        unsigned int &	        id( void ) 	        { return _id; }
        const unsigned int &	id( void ) const	{ return _id; }

        double &	            weight( void ) 	    { return _weight; }
        const double &	        weight( void ) const{ return _weight; }

        QString &	            text( void )        { return _text; }
        const QString &	        text( void ) const	{ return _text; }

        bool &	                textOn( void ) 	    { return _textOn; }
        const bool &	        textOn( void ) const{ return _textOn; }

        //------------------------------------------------------------------------------
        //	Special functions
        //------------------------------------------------------------------------------
        void    init      ( void );

    };

} // namespace Vector
} // namespace Ui

#endif // Ui_Vector_GraphicsEdgeItem_H
