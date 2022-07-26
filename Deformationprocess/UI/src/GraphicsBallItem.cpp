//******************************************************************************
// GraphicsBallItem.h
//	: header file for ball items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "GraphicsBallItem.h"

namespace Ui {
namespace Vector {

    //------------------------------------------------------------------------------
    //	Private functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Protected functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    //	Public functions
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    //	Constructors & Destructors
    //------------------------------------------------------------------------------
    //
    //  GraphicsBallItem::GraphicsBallItem -- constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsBallItem::GraphicsBallItem( QGraphicsItem *parent )
    {
        //setFlag( QGraphicsItem::ItemIsSelectable );
        //setFlag( QGraphicsItem::ItemIsMovable );
        //setFlag( QGraphicsItem::ItemSendsGeometryChanges );
        //setAcceptDrops( true );

        _textOn = false;
        _radius = 5;
    }

    void GraphicsBallItem::setRadius(double d ) 
    {
        _radius = d; 
    }

    //
    //  GraphicsBallItem::GraphicsBallItem -- parameterized constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsBallItem::GraphicsBallItem( qreal x, qreal y, qreal w, qreal h, QGraphicsItem *parent )
    {
        setRect( QRectF( x, y, w, h ) );
    }

    //
    //  GraphicsBallItem::GraphicsBallItem -- copy constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsBallItem::GraphicsBallItem( const QRectF &rect, QGraphicsItem *parent )
    {
        setRect( rect );
    }

    //
    //  GraphicsBallItem::init -- initialization
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void GraphicsBallItem::init ( void )
    {
    }

    //
    //  GraphicsBallItem::boundingRect -- find bounding box
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the bounding box
    //
    QRectF GraphicsBallItem::boundingRect( void ) const
    {
        return rect();
    }

    //
    //  GraphicsBallItem::paint -- paint the object
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void GraphicsBallItem::paint( QPainter *painter, const QStyleOptionGraphicsItem *option,
                                  QWidget *widget )
    {
        _font_size = 12;
        _font = QFont( "Arial", _font_size, QFont::Bold, false );

        // draw boundary
        //rect().setX( rect().x() - sx );
        //rect().setY( rect().y() - sy );
        QRectF fineRect( rect() );
        fineRect.setX( fineRect.x()-_radius );
        fineRect.setY( fineRect.y()-_radius );
        fineRect.setWidth( 2.0*_radius );
        fineRect.setHeight( 2.0*_radius );
        painter->setRenderHints( QPainter::Antialiasing );
        painter->setPen( pen() );
        painter->setBrush( brush() );
        painter->setFont( _font );

        painter->drawEllipse( fineRect );

        //cerr << "id = " << _id << endl;
        if( _textOn == true ){
            painter->drawText( rect().x()+10, rect().y()-10, _text );
            //painter->drawText( rect().x()+5, rect().y()-5, QString::fromStdString( to_string( _id) ) );
            //painter->drawText( rect().x()+10, rect().y()-10, _name );
        }

        //cerr << "paint x = " << pos().x() << " y = " << pos().y() << endl;

        // Qt function
        //if ( option->state & QStyle::State_Selected )
        //	qt_graphicsItem_highlightSelected( this, painter, option );
        // cerr << "painting ball..." << endl;
    }

    //
    //  GraphicsBallItem::paint -- find the object type
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the object type
    //
    int GraphicsBallItem::type( void ) const
    {
        return 0; //GRAPHICS_BALL+QGraphicsItem::UserType;
    }

} // namespace Vector
} // namespace Ui
