//******************************************************************************
// GraphicsView.cpp
//	: program file for graphics view
//
//------------------------------------------------------------------------------
//
//	Ver 2.00		Date: Sun Mar 14 20:00:00 2021
//
//******************************************************************************

//------------------------------------------------------------------------------
//	Including Header Files
//------------------------------------------------------------------------------

#include "GraphicsBallItem.h"
#include "GraphicsEdgeItem.h"
#include "GraphicsView.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Ui {
	namespace Vector {
		
		//------------------------------------------------------------------------------
		//	Private functions
		//------------------------------------------------------------------------------
		//
		//  GraphicsView::_init -- initialize data
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::_init( void )
		{
		}
		
		//
		//  GraphicsView::_clear -- clear data
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::_clear( void )
		{
		}

		
		//------------------------------------------------------------------------------
		//	Protected functions
		//------------------------------------------------------------------------------
		//
		//  GraphicsView::_item_metro -- draw metro
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::_item_metro( void )
		{
			if( _metroPtr == nullptr ) return;
			_scenePtr->clear();
			
			UndirectedGraph &g = _metroPtr->g();

			//Draw Guide
			if (_guidePtr != nullptr) {
				UndirectedGraph &guideG = _guidePtr->g();

				//Draw Matching Lines
				vector<VertexDescriptor> matchGuide = _matchingPtr->getGuides(); 
				vector<VertexDescriptor> matchStat = _matchingPtr->getStations(); 

				if (matchStat.size() > 0 ) {
					QPainterPath matchPath; 
					
					for (int i = 0; i < matchStat.size(); i++) {
						Coord2 gCoord = *guideG[matchGuide[i]].coordPtr; 
						Coord2 sCoord = *g[matchStat[i]].coordPtr; 
						matchPath.moveTo(gCoord.x(), -gCoord.y()); 
						matchPath.lineTo(sCoord.x(), -sCoord.y()); 
					}

					//Draw line from point to cloestst point on Guide;
                    BGL_FORALL_VERTICES(vertex, g, UndirectedGraph) {
					    Coord2 s = *g[vertex].coordPtr;
					    // Coord2 t = g[vertex].closestPointOnEdge;
					    Coord2 t = g[vertex].closOnMetro;
					    if ((t - Coord2(0,0)).norm() > 0.5 ) {
                            matchPath.moveTo(s.x(), - s.y());
                            matchPath.lineTo(t.x(), - t.y());
					    }
					}
				
					GraphicsEdgeItem *itemptr = new GraphicsEdgeItem;
					QPen guidePen = QPen( QColor( 150, 150, 150, _guideOpacity), 2 ); 
					guidePen.setStyle(Qt::DashLine); 
					itemptr->setPen( guidePen );
					itemptr->setPath( matchPath );
					itemptr->weight() = 1.0;
					
					_scenePtr->addItem( itemptr );
				}

				//draw Guides
				QPainterPath path; 
				BGL_FORALL_EDGES ( ed, guideG, UndirectedGraph) {
					VertexDescriptor s = source(ed, guideG); 
					VertexDescriptor t = target(ed, guideG); 
					Coord2 sCord = *guideG[s].coordPtr; 
					Coord2 tCord = *guideG[t].coordPtr; 
					path.moveTo(sCord.x(), -sCord.y());
					path.lineTo(tCord.x(), -tCord.y());
				}
				GraphicsEdgeItem *itemptr = new GraphicsEdgeItem;
				itemptr->setPen( QPen( QColor(100, 100, 100, _guideOpacity), 2 ) );
				itemptr->setBrush( QBrush( QColor( 255, 255, 255, 0 ), Qt::SolidPattern ) );
				itemptr->setPath( path );
				itemptr->weight() = 1.0;
				
				_scenePtr->addItem( itemptr );

				BGL_FORALL_VERTICES(vd, guideG, UndirectedGraph) {
					GraphicsBallItem *verItemptr = new GraphicsBallItem;
					DegreeSizeType degrees = out_degree( vd, guideG );
					verItemptr->setRadius(1.5);
					if (degrees > 2) {
					    verItemptr->setBrush( QBrush( QColor( 255, 100, 255, _guideOpacity ), Qt::SolidPattern ) );
					} else {
					    auto outerPath = _guidePtr->getOutherPath();
					    if(std::find(outerPath.begin(), outerPath.end(), vd) != outerPath.end()) {
					        verItemptr->setBrush( QBrush( QColor( 100, 100, 255, _guideOpacity ), Qt::SolidPattern ) );
					    } else {
					        verItemptr->setBrush( QBrush( QColor( 100, 255, 100, _guideOpacity ), Qt::SolidPattern ) );
					    }
					}
					verItemptr->setPen( QPen( QColor( 100, 100, 100, 0 ), 0 ) );
					verItemptr->setRect( QRectF( guideG[ vd ].coordPtr->x(), - guideG[ vd ].coordPtr->y(), 12, 12 ) );
					_scenePtr->addItem( verItemptr );
				}
			} 
			
			// draw edges
			BGL_FORALL_EDGES( ed, g, UndirectedGraph ) {
				
				VertexDescriptor vdS = source( ed, g );
				VertexDescriptor vdT = target( ed, g );
				Coord2 toward = *g[vdT].coordPtr - *g[vdS].coordPtr;
				toward.normalize();
				Coord2 cross( toward.y(), -toward.x() );
				
				const unsigned int & N = g[ ed ].lineID.size();
				// const double lineWidth = 10.0/(double)N;
				for ( unsigned int k = 0; k < N; ++k ) {

					GraphicsEdgeItem *itemptr = new GraphicsEdgeItem;
					double shift = ( ( ( double ) k + .5 ) / ( double ) N - .5 ) * _lineWidth;
					if ( (double)k > ((double) N) * 0.5) {
					    shift *= -1;
					}

					QPen pen = QPen( QColor( _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 0 ]*255,
                                             _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 1 ]*255,
                                             _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*255, 120 ), _lineWidth * 0.65 );
                                                // _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*255, 220. ), _lineWidth * 1. );
                    if (posterStyle) {
                        pen.setColor(QColor( _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 0 ]*255,
                                         _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 1 ]*255,
                                         _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*255, 180 ));
                    }
                    double n = N;
                    // if (n > 1 ) n *= .75;
#ifdef DOTTEDLINES
					if (!g[ed].matchPath) {
                        if (!g[ed].isMetro) {
                            pen.setStyle(Qt::DashLine);
                            pen.setWidth(_lineWidth * 0.005);
                        } else {
                            if (highlightMatchPath) pen.setWidth(_lineWidth * 0.2);
                        }
					}
#endif
                    pen.setWidth(_lineWidth * 1.0);
#ifdef HIGHLIGHTMATCHPATH
					if (_highlightPath) {
                        if (g[ed].matchPath) {
                            pen.setColor(QColor( _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 0 ]*255,
                                                 _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 1 ]*255,
                                                 _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*255, 255. ));
                            pen.setWidth(_lineWidthHighlight);
                            shift = ( ( ( double ) k + .5 ) / ( double ) N - .5 ) *_lineWidthHighlight;
                            if ( (double)k > ((double) N) * 0.5) {
                                shift *= -1;
                            }
                        } else if (!g[ed].isMetro) {
                            pen.setWidth(_lineWidth * .2 );
                        }
                        else {
                            pen.setWidth(_lineWidth * 1.0);
                        }
					}
#endif

#ifdef HIGHLIGHTSMOOTHPATH
                    if (_highlightSmoothPath) {
					    if (g[ed].isSmoothPath) {
					        pen.setColor(QColor( _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 0 ]*255,
                                                 _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 1 ]*255,
                                                 _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*255, 255. ));
					        pen.setWidth(_lineWidthHighlight);
					        shift = ( ( ( double ) k + .5 ) / ( double ) N - .5 ) *_lineWidthHighlight;
					        if ( (double)k > ((double) N ) * 0.5) {
					            shift *= -1;
					        }
					    } else if (!g[ed].isMetro) {
					        pen.setWidth(_lineWidth * .2);
					    }
					    else {
					        pen.setWidth(_lineWidth * 1.0);
					    }
					}

#endif
                    if (!g[ed].isVisible) {
                        pen.setWidth(_lineWidth * 0.0);
                        pen.setColor(QColor( _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 0 ]*0,
                                             _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 1 ]*0,
                                             _metroPtr->lineColor( g[ ed ].lineID[ k ] )[ 2 ]*0, 0 ));
                    }
#ifdef BLACKWHITE
                    pen.setColor(QColor(0, 0, 0, 130));
					if (!g[ed].matchPath) pen.setColor(QColor(0,0,0,100));
                    if (!g[ed].isMetro && !g[ed].matchPath) pen.setColor(QColor(0,0,0, 80));
#endif
					pen.setCapStyle(Qt::RoundCap); 
					itemptr->setPen( pen );
					
					Coord2 coordS = *g[ vdS ].coordPtr;
					Coord2 coordT = *g[ vdT ].coordPtr;
					Coord2 src = coordS + shift*cross;
					Coord2 tar = coordT + shift*cross;

					QPainterPath path;
					path.moveTo( src.x(), -src.y() ); 
					QPointF i( src.x(), -src.y() ); 
					QPointF j( tar.x(), -tar.y() );

					Coord2 mCo = src; 
					Coord2 kCo = tar; 
					Coord2 t1 = mCo; 
					Coord2 t1_a = mCo; 
					Coord2 t1_b = mCo; 
					Coord2 t2 = kCo;
					Coord2 t2_a = kCo; 
					Coord2 t2_b = kCo;  

					Coord2 ij = tar - src; 
					Coord2 ji = src - tar; 
					double magIJ = sqrt(ij.x() * ij.x() + ij.y() * ij.y()); 
					
					//Fing c1 
					DegreeSizeType degrees = out_degree( vdS, g );
        			if( degrees <= 2 && g[ed].isCurved && false ){
						double magIM = 1.0; 
						OutEdgeIterator e, e_end;
						for ( tie( e, e_end ) = out_edges( vdS, g ); e != e_end; ++e ) {
							EdgeDescriptor edI = *e;
							VertexDescriptor vS = source( edI, g );
							VertexDescriptor vT = target( edI, g );
							if ( edI != ed ) {
								if ( vdS == vS ) {
									mCo = *g[ vT ].coordPtr; 
								} else {
									mCo = *g[ vS ].coordPtr; 
								}
								Coord2 im = src - mCo; 
								magIM = sqrt(im.x() * im.x() + im.y() * im.y()); 
								break; 
							}		
						}
						double dIM = (magIM + magIJ * 0.4) / magIM; 
						t1_a = src - mCo; 
						t1_a = Coord2(t1_a.x() * dIM, t1_a.y() * dIM); 
						t1_a = t1_a + mCo; 
					}
					t1_b = Coord2(ij.x() * 0.4, ij.y() * 0.4); 
					t1_b = t1_b + src;
					t1 = t1_b - t1_a; 
					t1 = Coord2(t1.x() * 0.5, t1.y()* 0.5); 
					t1 = t1_a + t1;  

					//Find c2 
					degrees = out_degree( vdT, g ); 
					if (degrees <= 2 && g[ed].isCurved && false ) {
						double magJK = 1.0; 
						OutEdgeIterator e, e_end;
						for ( tie( e, e_end ) = out_edges( vdT, g ); e != e_end; ++e ) {
							EdgeDescriptor edJ = *e;
							VertexDescriptor vS = source( edJ, g );
							VertexDescriptor vT = target( edJ, g );
							if ( edJ != ed ) {
								if ( vdT == vT ) {
									kCo = *g[ vS ].coordPtr; 
								} else {
									kCo = *g[ vT ].coordPtr; 
								}
								Coord2 jk = tar - kCo; 
								magJK = sqrt(jk.x() *jk.x() + jk.y() * jk.y()); 
								break; 
							}		
						}
						double dJK = (magJK + magIJ * 0.4 ) / magJK; 
						t2_a = tar - kCo; 
						t2_a = Coord2(t2_a.x() * dJK, t2_a.y() * dJK ); 
						t2_a = t2_a + kCo;
						t2_b = Coord2(ji.x() * 0.4, ji.y() * 0.4 );
						t2_b = t2_b + tar; 
						t2 = t2_b - t2_a; 
						t2 = Coord2(t2.x() * 0.5, t2.y() * 0.5 ); 
						t2 = t2_a + t2; 
					}
					
					QPointF t1Q(t1.x(), -t1.y()); 
					QPointF t2Q(t2.x(), -t2.y()); 

					path.cubicTo(t1Q, t2Q,j); 
					itemptr->id() = g[ ed ].id;
					itemptr->setPath( path );
					itemptr->weight() = g[ ed ].weight;
					
					_scenePtr->addItem( itemptr );
				}
			}
			
			// draw stations
			BGL_FORALL_VERTICES( vd, g, UndirectedGraph ) {
			    DegreeSizeType degrees = out_degree( vd, g );

			    auto dummyStrign = "Dummystation";
			    bool isDummyNode = false;

			    // isDummyNode = false;
			    if (g[vd].isStation and degrees > 0) { // and _metroPtr->isOctolinearTurningPoint(vd)
			        bool isDummyNode = false;
			        if ((*g[vd].namePtr).find("Dummystation") == 0)
			            isDummyNode = true;
			        if (!isDummyNode) {
			            GraphicsBallItem *itemptr = new GraphicsBallItem;
			            itemptr->id() = g[ vd ].id;
			            double calculatedRadius = (_lineWidth);


			            if (posterStyle) {
			                if (g[vd].interstate) {
			                    GraphicsBallItem *outlinePtr = new GraphicsBallItem;
			                    outlinePtr->id() = g[vd].id;
			                    outlinePtr->setRadius(calculatedRadius * 1);
			                    outlinePtr->setPen(QPen(QColor(0,0,0,255), .7));
			                    outlinePtr->setBrush(QBrush(QColor(255,255,255,255), Qt::SolidPattern));
			                    outlinePtr->setRect( QRectF( g[ vd ].coordPtr->x(), -g[ vd ].coordPtr->y(), 50, 50 ) );
			                    // _scenePtr->addItem(outlinePtr);

                                itemptr->setRadius(calculatedRadius);
                                itemptr->setPen(QPen(QColor(0,0,0,255), _lineWidth * 0.4));
                                itemptr->setBrush( QBrush( QColor( 255, 255, 255, 255 ), Qt::SolidPattern )  );
                                itemptr->setRect( QRectF( g[ vd ].coordPtr->x(), -g[ vd ].coordPtr->y(), 50, 50 ) );
                                _scenePtr->addItem( itemptr );



			                } else {
			                    OutEdgeIterator e, e_end;
			                    tie(e, e_end) = out_edges(vd, g);
			                    EdgeDescriptor ed = *e;
			                    QColor strokeColor;
			                    if (_metroPtr->lineColor( g[ed].lineID.size() == 1)) {
			                        strokeColor = QColor(_metroPtr->lineColor( g[ed].lineID[0])[0] * 100,
                                                         _metroPtr->lineColor( g[ed].lineID[0])[1] * 100,
                                                         _metroPtr->lineColor( g[ed].lineID[0])[2] * 100,
                                                         255);
			                    } else {
			                        strokeColor = QColor(0,0,0,255);
			                    }

			                    itemptr->setRadius(calculatedRadius);
                                itemptr->setPen(QPen(QColor(0,0,0,255), _lineWidth * 0.4));
			                    itemptr->setBrush( QBrush( QColor( 255, 255, 255, 255 ), Qt::SolidPattern )  );
			                    itemptr->setRect( QRectF( g[ vd ].coordPtr->x(), -g[ vd ].coordPtr->y(), 50, 50 ) );
			                    _scenePtr->addItem( itemptr );
			                }
			            } else { // Not poster style for paper
			                itemptr->setRadius(calculatedRadius);
			                itemptr->setPen( QPen( QColor( 0, 0, 0, 0 ), 0 ) );

			                // if (g[vd].smoothPath) {
			                // if (g[vd].intersection) {
			                // if (g[vd].collision ) {
			                if (false) {
			                    itemptr->setBrush( QBrush( QColor( 255, 0, 0, 255 ), Qt::SolidPattern ) );
			                } else {
			                    // itemptr->setRadius(calculatedRadius * 0.65);
			                    itemptr->setBrush( QBrush( QColor( 0, 0, 0, 180 ), Qt::SolidPattern ) );
			                }

			                itemptr->setRect( QRectF( g[ vd ].coordPtr->x(), -g[ vd ].coordPtr->y(), 50, 50 ) );
			                _scenePtr->addItem( itemptr ); // ADD black stroke to stations
			                GraphicsBallItem *whiteptr = new GraphicsBallItem;
			                whiteptr->id() = g[vd].id;
			                whiteptr->setRadius(calculatedRadius * 1);
			                whiteptr->setPen(QPen(QColor(255, 255, 255, 0), 0));

			                whiteptr->setBrush( QBrush( QColor( 255, 255, 255, 255 ), Qt::SolidPattern ) );
			                whiteptr->setRect( QRectF( g[ vd ].coordPtr->x(), -g[ vd ].coordPtr->y(), 50, 50 ) );
			                _scenePtr->addItem( whiteptr );
			            }
			        }

                } else {
			        // cout << "not a station" << endl;
			    }
			}

		}
		void GraphicsView::_calcLineWidth( double nStations ) 
		{
            cout << "Calculating Line Width" << endl;
			_lineWidth = 13.0 / pow((nStations / 50.0), (1.0 / 1.8 ));
			_lineWidth = 0.08 * 1000 / pow(nStations, 1.0 / 2.0 );
			_lineWidth *= 0.7;
			// _lineWidth *= 0.7; // For tokyo –> thinner lines needed
			if (!posterStyle) _lineWidth *= 0.8;
			else _lineWidth *= 0.8;


            _lineWidth = 3.0;
            _lineWidthHighlight = 6.0;
		}
		
		//
		//  GraphicsView::_item_metro -- update metro
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::_update_item_metro( void )
		{
		}
		
		//------------------------------------------------------------------------------
		//	Public functions
		//------------------------------------------------------------------------------
		
		//------------------------------------------------------------------------------
		//	Constructors & Destructors
		//------------------------------------------------------------------------------
		//
		//  GraphicsView::GraphicsView -- constructor
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		GraphicsView::GraphicsView( QWidget *parent )
				: QGraphicsView( parent ) {
			//------------------------------------------------------------------------------
			// configuration file
			//------------------------------------------------------------------------------
			_metroPtr = nullptr;
			_basePtr = nullptr;
			_scenePtr = new QGraphicsScene;
            _scenePtr->setSceneRect( -ZuKai::Base::Common::getMainwidgetWidth() / 2.0, -ZuKai::Base::Common::getMainwidgetHeight() / 2.0,
                                     ZuKai::Base::Common::getMainwidgetWidth(), ZuKai::Base::Common::getMainwidgetHeight() );  // x, y, w, h

			this->setScene( _scenePtr );
		}
		

		//------------------------------------------------------------------------------
		//	Event handlers
		//------------------------------------------------------------------------------
		void GraphicsView::simulateKey( Qt::Key key ) {

			// press the key
			QKeyEvent eventP( QEvent::KeyPress, key, Qt::NoModifier );
			QApplication::sendEvent( this, &eventP );
			// release the key
			QKeyEvent eventR( QEvent::KeyRelease, key, Qt::NoModifier );
			QApplication::sendEvent( this, &eventR );
		}
		
		//
		//  GraphicsView::keyPressEvent -- key press event
		//
		//  Inputs
		//  event: key event
		//
		//  Outputs
		//  none
		//
		void GraphicsView::keyPressEvent( QKeyEvent *event )
		{
			QGraphicsView::keyPressEvent( event );
			updateSceneItems();
		}
		
		//
		//  GraphicsView::mousePressEvent -- mouse press event
		//
		//  Inputs
		//  event: mouse event
		//
		//  Outputs
		//  none
		//
		void GraphicsView::mousePressEvent( QMouseEvent *event ) {
            // TODO
			QPointF coord = mapToScene( event->pos() );
			lastTrackedPointerPosition = coord;
			if (_selectionMode == 0) {
			    bool removed = _basePtr->removeNode(Coord2(coord.x(), -coord.y()));
			    if (!removed) _basePtr->addNode( Coord2(coord.x(), -coord.y()) );
			}
            else if ( _selectionMode == 1 ) {
                _basePtr->manuallyAddVertexToPath(Coord2(coord.x(), -coord.y()));
            }

			_item_metro(); 
			double square = 10.0;
			switch( event->buttons() ) {
			
			case Qt::RightButton:
				break;
			case Qt::LeftButton:
				break;
			case Qt::MiddleButton:
				break;
			default: {
				break;
			}
			}
			
			QGraphicsView::mousePressEvent( event );
		}
		
		//
		//  GraphicsView::mouseMoveEvent -- mouse move event
		//
		//  Inputs
		//  event: mouse event
		//
		//  Outputs
		//  none
		//
		void GraphicsView::mouseMoveEvent( QMouseEvent *event ) {
			double norm = QLineF( _cursor, _oldCursor ).length();
			switch( event->buttons() ) {
			    case Qt::LeftButton:
			        if (_selectionMode == 2 ) {
                        QPointF screenCoord = mapToScene( event->pos() );
                        QPointF delta = screenCoord - lastTrackedPointerPosition;

                        Coord2 translate = Coord2(delta.x(), -delta.y());
                        _guidePtr->translate(translate);
                        lastTrackedPointerPosition = screenCoord;
                        _item_metro();
			        }
			        break;
			}
			// if( norm > 1.0 * _min_point_distance ) {
			// 	simulateKey( Qt::Key_Z );
			//	_oldCursor = _cursor;
			// }
			QGraphicsView::mouseMoveEvent( event );
		}
		
		//
		//  GraphicsView::mouseReleaseEvent -- mouse release event
		//
		//  Inputs
		//  event: mouse event
		//
		//  Outputs
		//  none
		//
		void GraphicsView::mouseReleaseEvent( QMouseEvent *event ) {
			_right_button_pressed = false;
			_left_button_pressed = false;
			_middle_button_pressed = false;
			
			QGraphicsView::mouseReleaseEvent( event );
		}
		
		
		//------------------------------------------------------------------------------
		//	Special functions
		//------------------------------------------------------------------------------
		//
		//  GraphicsView::initSceneItems -- initialize SceneItems
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::initSceneItems( void ) {
			_scenePtr->clear();
			
            _item_metro();
		}
		
		//
		//  GraphicsView::updateSceneItems -- initialize SceneItems
		//
		//  Inputs
		//  none
		//
		//  Outputs
		//  none
		//
		void GraphicsView::updateSceneItems( void ) {
			scene()->update();
		}
		
		//
		//  GraphicsView::exportPNG -- export png
		//
		//  Inputs
		//  x, y, w, h: center, width, height
		//
		//  Outputs
		//  none
		//
		void GraphicsView::exportPNG( double w, double h ) {
			
			QString path = QString("Metro.png"); 

			QString newPath = QFileDialog::getSaveFileName(this, tr("Save PNG"),
        			path, tr("PNG files (*.png)"));

    		if (newPath.isEmpty())
        		return;

    		path = newPath;
			QImage screenshot( w, h, QImage::Format_RGB32 ); // maximum 32767x32767
			
			QPainter painter( &screenshot );
			painter.setRenderHint( QPainter::Antialiasing );
			painter.fillRect( 0, 0, w, h, Qt::white );
			_scenePtr->render( &painter );
			screenshot.save( path );
		}

		void GraphicsView::exportPNG ( double w, double h, string path) {
            cout << "Export png to path: " << path << endl;
		    // exportPNG(w, h);
		    QString qPath = QString(QLatin1String(path.c_str()));
		    QImage screenshot( w, h, QImage::Format_RGB32 ); // maximum 32767x32767

		    QPainter painter( &screenshot );
		    painter.setRenderHint( QPainter::Antialiasing );
		    painter.fillRect( 0, 0, w, h, Qt::white );
		    _scenePtr->render( &painter );
		    screenshot.save( qPath );
		}
		
		//
		//  GraphicsView::exportSVG -- export svg
		//
		//  Inputs
		//  x, y, w, h: center, width, height
		//
		//  Outputs
		//  none
		//
		void GraphicsView::exportSVG(double w, double h ) {
			QString path = QString("Metro.svg");

			QString newPath = QFileDialog::getSaveFileName(this, tr("Save svg"),
        			path, tr("SVG files (*.svg)"));

    		if (newPath.isEmpty())
        		return;

    		path = newPath;

			QSvgGenerator generator; 
			generator.setFileName(path);
    		generator.setSize(QSize(w, h));
   	 		generator.setViewBox(QRect(0, 0, w, h));

			QPainter painter;
			painter.begin( &generator );
			_scenePtr->render( &painter );
			painter.end();
		}
		
	} // namespace Vector
} // namespace Ui