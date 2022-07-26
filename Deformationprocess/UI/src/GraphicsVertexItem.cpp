//******************************************************************************
// GraphicsVertexItem.cpp
//	: program file for vertex items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "GraphicsVertexItem.h"

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
    //  GraphicsVertexItem::GraphicsVertexItem -- constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsVertexItem::GraphicsVertexItem( QGraphicsItem *parent )
    {
        //setFlag( QGraphicsItem::ItemIsSelectable );
        //setFlag( QGraphicsItem::ItemIsMovable );
        //setFlag( QGraphicsItem::ItemSendsGeometryChanges );
        //setAcceptDrops( true );

        _radius = 10;
    }

    //
    //  GraphicsVertexItem::GraphicsVertexItem -- parameterized constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsVertexItem::GraphicsVertexItem( qreal x, qreal y, qreal w, qreal h, QGraphicsItem *parent )
    {
        setRect( QRectF( x, y, w, h ) );
    }

    //
    //  GraphicsVertexItem::GraphicsVertexItem -- copy constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsVertexItem::GraphicsVertexItem( const QRectF &rect, QGraphicsItem *parent )
    {
        setRect( rect );
    }

    //
    //  GraphicsVertexItem::init -- initialization
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void GraphicsVertexItem::init ( void )
    {
    }

    //
    //  GraphicsVertexItem::boundingRect -- find bounding box
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the bounding box
    //
    QRectF GraphicsVertexItem::boundingRect( void ) const
    {
        return rect();
    }

    //
    //  GraphicsVertexItem::paint -- paint the object
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  none
    //
    void GraphicsVertexItem::paint( QPainter *painter, const QStyleOptionGraphicsItem *option,
                                    QWidget *widget )
    {
        _font_size = 12;
        _font = QFont( "Arial", _font_size, QFont::Bold, false );

        QFontMetrics metrics( _font );
        double sx = metrics.width( _name );
        double sy = 0.5*metrics.height();

        QRectF fineRect( rect() );
        fineRect.setX( rect().x() - 0.5*sx );
        fineRect.setY( rect().y() - sy );
        fineRect.setWidth( sx );
        fineRect.setHeight( 2.0*sy );

        painter->setRenderHints( QPainter::Antialiasing );
        painter->setPen( pen() );
        painter->setBrush( brush() );

        painter->drawRoundedRect( fineRect, 5, 5, Qt::AbsoluteSize );

        painter->setFont( _font );
        //cerr << "id = " << _id << endl;
        //painter->drawText( rect().x()+10, rect().y()-10, QString::fromStdString( to_string( _id ) ) );
        painter->drawText( fineRect.x()+0.5*( fineRect.width()-sx ),
                           fineRect.y()+0.5*( fineRect.height()+0.5*sy ), _name ); // (x ,y) must be left-upper corner
        //painter->drawText( fineRect.x()+0.5*( fineRect.width()-sx ) - 0.5*MIN_NEIGHBOR_DISTANCE,
        //                   fineRect.y()+0.5*( fineRect.height()+0.5*sy ) + MIN_NEIGHBOR_DISTANCE, _name );

        //cerr << "paint x = " << pos().x() << " y = " << pos().y() << endl;

        // Qt function
        //if ( option->state & QStyle::State_Selected )
        //	qt_graphicsItem_highlightSelected( this, painter, option );
        // cerr << "painting ball..." << endl;
    }

    //
    //  GraphicsVertexItem::type -- find type
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the type of the object
    //
    int GraphicsVertexItem::type( void ) const
    {
        return 0; //GRAPHICS_VERTEX+QGraphicsItem::UserType;
    }

} // namespace Vector
} // namespace Ui