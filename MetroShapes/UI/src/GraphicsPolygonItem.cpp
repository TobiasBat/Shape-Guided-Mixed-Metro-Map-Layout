//******************************************************************************
// GraphicsPolygonItem.cpp
//	: program file for polygon items
//
//------------------------------------------------------------------------------
//
//	Ver 1.00		Date: Tue Jun 19 02:36:37 2019
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "GraphicsPolygonItem.h"

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
    //  GraphicsPolygonItem::GraphicsPolygonItem -- constructor
    //
    //  Inputs
    //  parent: parent object
    //
    //  Outputs
    //  none
    //
    GraphicsPolygonItem::GraphicsPolygonItem( QGraphicsItem *parent )
    {
        _font_size = 12;
        _font = QFont( "Arial", _font_size, QFont::Bold, false );

        setFlag( QGraphicsItem::ItemIsSelectable );
        setFlag( QGraphicsItem::ItemIsMovable );
        setFlag( QGraphicsItem::ItemSendsGeometryChanges );
        //setAcceptDrops( true );

        //pen().setJoinStyle( Qt::MiterJoin );
        pen().setJoinStyle( Qt::RoundJoin );
        _textOn = false;
    }

    //
    //  GraphicsPolygonItem::boundingRect -- find bounding box
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the bounding box
    //
    QRectF GraphicsPolygonItem::boundingRect( void ) const
    {
#ifdef GRAPHICSPOLYGONITEM_DEBUG
        cerr << "id = " << _id
             << setw(10 )<< " bx = " << polygon().boundingRect().x()
             << setw(10 )<< " by = " << polygon().boundingRect().y()
             << setw(10 )<< " bw = " << polygon().boundingRect().width()
             << setw(10 )<< " bh = " << polygon().boundingRect().height()
             << endl;
#endif // GRAPHICSPOLYGONITEM_DEBUG
        return polygon().boundingRect();
    }

    //
    //  GraphicsPolygonItem::paint -- paint scene
    //
    //  Inputs
    //  painter: Qpainter
    //  option: QStyleOptionGraphicsItem
    //  widget: QWidget
    //
    //  Outputs
    //  none
    //
    void GraphicsPolygonItem::paint( QPainter *painter, const QStyleOptionGraphicsItem *option,
                                     QWidget *widget )
    {
        QFontMetrics metrics( _font );
        double sx = metrics.width( _text );
        double sy = 0.5*metrics.height();

        // draw boundary
        painter->setRenderHints( QPainter::Antialiasing );
        painter->setPen( pen() );
        painter->setBrush( brush() );
        painter->drawPolygon( polygon() );

        const QPolygonF &p = polygon();
        if( _textOn == true ){

            painter->setPen( QPen( QColor( 0 ,0, 0, 255 ), 4 ) );
#ifdef SKIP
            for( unsigned int i = 0; i < p.size(); i++ ){
                painter->drawText( p.at(i).x()+5, p.at(i).y()-5, QString::fromStdString( to_string( _id ) ) );
            }
#endif // SKIP

            // painter->drawText( p.at(0).x(), p.at(0).y(), _text );
            painter->drawText(
                    _bbox.x() + 0.5 * _bbox.width() - 0.5 * sx,
                    -_bbox.y() - 0.5 * _bbox.height() + 0.5 * sy,
                    _text
                    );
        }

        //cerr << "paint x = " << pos().x() << " y = " << pos().y() << endl;
        // Qt function
        //if ( option->state & QStyle::State_Selected )
        //	qt_graphicsItem_highlightSelected( this, painter, option );
    }

    //
    //  GraphicsPolygonItem::type -- find type
    //
    //  Inputs
    //  none
    //
    //  Outputs
    //  reference to the type of the object
    //
    int GraphicsPolygonItem::type( void ) const
    {
        return 0; //GRAPHICS_POLYGON+QGraphicsItem::UserType;
    }

} // namespace Vector
} // namespace Ui