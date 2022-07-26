#include "Widget.h"

// static variables
bool Widget::_sim_flag = false;

Widget::Widget( const QGLFormat& format, QWidget *parent )
    : QGLWidget( format, parent )
{
    _timer = new QBasicTimer();
}

Widget::~Widget()
{
}

void Widget::init( Metro * __metro, Metro * __simmetro )
{
    _metro = __metro;
    _simmetro = __simmetro;
    focusVD = NULL;
    sourceVD = NULL;
    targetVD = NULL;
    _last_pointer_x = 0;
    _last_pointer_y = 0;
    _time_step = 0;
    _isDeformed = false;
}

void Widget::timerEvent( QTimerEvent *e )
{
    _time_step++;
    update();
    //cerr << "time_type = " << _time_type << ", time_step = " << _time_step << endl;
    if( _time_step >= TIMER_STEP ) {
        _timer->stop();
        _isDeformed = true;
    }
}

void Widget::initializeGL( void )
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
    glClearColor( 1.0, 1.0, 1.0, 1.0 );

    cerr << " Support OpenGL: " << glGetString( GL_VERSION ) << endl;
    cerr << " Support GLSL: " << glGetString( GL_SHADING_LANGUAGE_VERSION ) << endl;
}

void Widget::paintGL( void )
{
    //for( int i = 0; i < 100; i++ ){
    // initialization
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glClearColor( 1.0, 1.0, 1.0, 1.0 );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glPushMatrix();
    glOrtho( -width()/2, width()/2, -height()/2, height()/2, -2.0, 2.0 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    //_paintMetro( GEOGRAPHY );
    //_paintMetro( SMOOTH );
    //animationSmooth();
    //animationOctilinear();
    _paintMetro( OCTILINEAR );
    if( _isDeformed == true ) _paintLabels();
    //_paintName();
    _paintMagnification();

    glFlush();
    //glutSwapBuffers();
    //usleep( 10000 );
    //}
}

void Widget::capture( char * name )
{
    static unsigned char        *data = NULL;
    int                         h = 800;
    int                         w = 600;
    static IplImage*            ptrImage = NULL;        // Mesh image

    paintGL();
    //display();

    if ( data == NULL ) data = new unsigned char [ w * h * 3 ];
    glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, data );

    if ( ptrImage != NULL ) {
        cvReleaseImage( &ptrImage );
        ptrImage = NULL;
    }
    ptrImage = cvCreateImage( cvSize( w, h ), IPL_DEPTH_8U, 3 );
    memcpy( ptrImage->imageData, data, ptrImage->imageSize );

    cvCvtColor( ptrImage, ptrImage, CV_BGR2RGB );
    cvFlip( ptrImage, NULL, 0 );
    cvSaveImage( name, ptrImage );


    cerr << "Capturing the main window ... done." << endl;
}

void Widget::_paintDisk( const double & radius )
{
    const int nDiv = 36;
    glBegin( GL_POLYGON );
    for ( int k = 0; k <= nDiv; k++ ) {
        double theta = 2.0*M_PI*(double)k/(double)nDiv;
        glVertex2d( radius*cos( theta ), radius*sin( theta ) );
    }
    glEnd();
}

void Widget::_paintCircle( const double & radius )
{
    const int nDiv = 36;
    glBegin( GL_LINE_LOOP );
    for ( int k = 0; k <= nDiv; k++ ) {
        double theta = 2.0*M_PI*(double)k/(double)nDiv;
        glVertex2d( radius*cos( theta ), radius*sin( theta ) );
    }
    glEnd();
}

void Widget::_paintStation( VertexDescriptor vd, METROTYPE type )
{
    Metro * ptrM = _metro;
    Graph * g = & _metro->g();
    if( _sim_flag == true ) {
        ptrM = _simmetro;
        g = & _simmetro->g();
    }
    else {
        ptrM = _metro;
        g = & _metro->g();
    }

    VertexSmoothMap         vertexSmooth        = get( vertex_mysmooth, *g );
    VertexTempMap           vertexTemp          = get( vertex_mytemp, *g );
    VertexCoordMap          vertexCoord         = get( vertex_mycoord, *g );
    VertexScaleMap          vertexScale         = get( vertex_myscale, *g );
    VertexSelectMagMap      vertexSelectMag     = get( vertex_myselectmag, *g );
    DegreeSizeType          degrees             = out_degree( vd, *g );

    Coord2 coord;
    if( _time_type == SMOOTH ) {
        coord = vertexTemp[ vd ] + (double)_time_step*( vertexSmooth[ vd ] - vertexTemp[ vd ] )/(double)( TIMER_STEP-1 );
    }
    else{
        coord = vertexSmooth[ vd ] + (double)_time_step*( vertexCoord[ vd ] - vertexSmooth[ vd ] )/(double)( TIMER_STEP-1 );
    }

    double x = coord.x();
    double y = coord.y();
    double radius = 4.0; // * ( double ) min_grid_size / DEFAULT_GRIDSIZE;
    double scale = vertexScale[ vd ];

    // cerr << "VID = " <<  vertexID[ vd ] << " Coord = " << vertexCoord[ vd ];

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glTranslated( x, y, 0.0 );
    glScaled( radius, radius, 1.0 );

    if ( degrees > 2 ) {
        // draw the outermost disk
        glColor4d( 0.2, 0.2, 0.2, 1.0 );
        _paintDisk( OUTDISK_LARGE * scale );
        // draw inner disk
        glColor4d( 0.8, 0.8, 0.8, 1.0 );
        _paintDisk( 1.9 * scale );
    }

    // draw the outermost disk
    glColor4d( 0.2, 0.2, 0.2, 1.0 );
    _paintDisk( OUTDISK_SMALL * scale );
    // draw inner disk
    if ( vertexSelectMag[ vd ] == true ) 
        glColor4d( 0.8, 0.4, 0.4, 1.0 );
    else 
        glColor4d( 0.8, 0.8, 0.8, 1.0 );
    _paintDisk( 1.0 * scale );

    glPopMatrix();
}

void Widget::animationSmooth( void )
{
    cerr << "animationSmooth!" << endl;
    repaint();
}

void Widget::animationOctilinear( void )
{
    cerr << "animationOctilinear!" << endl;
    repaint();
}

void Widget::_paintMetro( METROTYPE type )
{
    Metro * ptrM = _metro;
    Graph * g = & _metro->g();
    if( _sim_flag == true ) {
        ptrM = _simmetro;
        g = & _simmetro->g();
    }
    else { 
        ptrM = _metro;
        g = & _metro->g();
    }

    VertexIDMap         vertexID        = get( vertex_myid, *g );
    VertexNameMap       vertexName      = get( vertex_myname, *g );
    VertexCoordMap      vertexCoord     = get( vertex_mycoord, *g );
    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, *g );
    VertexGeoMap        vertexGeo       = get( vertex_mygeo, *g );
    VertexTempMap       vertexTemp      = get( vertex_mytemp, *g );
    VertexLineIDMap     vertexLineID    = get( vertex_mylineid, *g );
    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, *g );
    VertexSelectTopMap  vertexSelectTop = get( vertex_myselecttop, *g );
    VertexExternalMap   vertexExternal  = get( vertex_myexternal, *g );
    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, *g );
    VertexIsStationMap  vertexIsStation = get( vertex_myisstation, *g );
    EdgeIDMap           edgeIndex       = get( edge_myid, *g );
    EdgeLineIDMap       edgeLineID      = get( edge_mylineid, *g );
    EdgeSelectShiftMap  edgeSelectShift = get( edge_myselectshift, *g );
    EdgeSelectCtrlMap   edgeSelectCtrl  = get( edge_myselectctrl, *g );

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

#define RENDER
#ifdef  RENDER
    // draw edges
    const double lineWidth = 10.0;
    BGL_FORALL_EDGES( edge, *g, Graph )
    {
        VertexDescriptor vS = source( edge, *g );
        VertexDescriptor vT = target( edge, *g );

        Coord2 toward = vertexCoord[vT] - vertexCoord[vS];
        toward.normalize();
        Coord2 cross( toward.y(), -toward.x() );

        const unsigned int & N = edgeLineID[ edge ].size();

        // Drawing edges
        glLineWidth( lineWidth/(double)N );
        for ( unsigned int k = 0; k < N; ++k ) {
            double shift = ( (( double )k+0.5)/( double )N - 0.5 ) * lineWidth;

            if( edgeSelectShift[ edge ] == true )
                glColor4d( 0.8, 0.4, 0.4, 0.7 );
            else if( edgeSelectCtrl[ edge ] == true )
                glColor4d( ptrM->lineColor( edgeLineID[edge][ k ] )[ 0 ],
                           ptrM->lineColor( edgeLineID[edge][ k ] )[ 1 ],
                           ptrM->lineColor( edgeLineID[edge][ k ] )[ 2 ],
                        1.0 );
            else
                glColor4d( ptrM->lineColor( edgeLineID[edge][ k ] )[ 0 ],
                           ptrM->lineColor( edgeLineID[edge][ k ] )[ 1 ],
                           ptrM->lineColor( edgeLineID[edge][ k ] )[ 2 ],
                        0.5 );

            Coord2 coordS, coordT;
            if( _time_type == SMOOTH ) {
                coordS = vertexTemp[ vS ] + (double)_time_step*( vertexSmooth[ vS ] - vertexTemp[ vS ] )/(double)(TIMER_STEP-1);
                coordT = vertexTemp[ vT ] + (double)_time_step*( vertexSmooth[ vT ] - vertexTemp[ vT ] )/(double)(TIMER_STEP-1);
            }
            else{
                coordS = vertexSmooth[ vS ] + (double)_time_step*( vertexCoord[ vS ] - vertexSmooth[ vS ] )/(double)(TIMER_STEP-1);
                coordT = vertexSmooth[ vT ] + (double)_time_step*( vertexCoord[ vT ] - vertexSmooth[ vT ] )/(double)(TIMER_STEP-1);
            }
            Coord2 src = coordS + shift*cross;
            Coord2 tar = coordT + shift*cross;

            glBegin( GL_LINES );
            glVertex2dv( src.element() );
            glVertex2dv( tar.element() );
            glEnd();
        }
    }

    // draw vertices
    BGL_FORALL_VERTICES( vertex, *g, Graph )
    {
        if( vertexIsStation[ vertex ] == true )
            _paintStation( vertex, type );
    }
#endif  // RENDER

#ifndef  RENDER
    // input 
    glLineWidth( 5.0 );
    glPointSize( 10.0 );
    switch( type ){
      case GEOGRAPHY:
        glColor4d( 0.6, 0.8, 0.4, 1.0 );
        break;
      case SMOOTH: 
        glColor4d( 0.0, 0.0, 0.6, 1.0 );
        break;
      case OCTILINEAR:
        glColor4d( 0.6, 0.0, 0.0, 1.0 );
        break;
      default:
        cerr << "something wrong here =_=, in " << __LINE__ << __FILE__ << endl;
        break;
    }

    // draw edges
    BGL_FORALL_EDGES( edge, *g, Graph ){

        VertexDescriptor vdS = source( edge, *g );
        VertexDescriptor vdT = target( edge, *g );

        Coord2 coordS, coordT;
        switch( type ){
          case GEOGRAPHY:
            coordS = vertexGeo[ vdS ];
            coordT = vertexGeo[ vdT ];
            break;
          case SMOOTH:
            coordS = vertexSmooth[ vdS ];
            coordT = vertexSmooth[ vdT ];
            break;
          case OCTILINEAR:
            if( vertexSelectMag[ vdS ] == true && vertexSelectMag[ vdT ] == true )
                glColor4d( 0.6, 0.0, 0.0, 1.0 );
            //else if( vertexSelectTop[ vdS ] == true && vertexSelectTop[ vdT ] == true )
            //    glColor4d( 0.6, 0.0, 0.0, 0.5 );
            else
                glColor4d( 0.6, 0.0, 0.0, 0.2 );
            if( _time_type == SMOOTH ) {
                coordS = vertexTemp[ vdS ] + (double)_time_step*( vertexSmooth[ vdS ] - vertexTemp[ vdS ] )/(double)(TIMER_STEP-1);
                coordT = vertexTemp[ vdT ] + (double)_time_step*( vertexSmooth[ vdT ] - vertexTemp[ vdT ] )/(double)(TIMER_STEP-1);
            }
            else{
                coordS = vertexSmooth[ vdS ] + (double)_time_step*( vertexCoord[ vdS ] - vertexSmooth[ vdS ] )/(double)(TIMER_STEP-1);
                coordT = vertexSmooth[ vdT ] + (double)_time_step*( vertexCoord[ vdT ] - vertexSmooth[ vdT ] )/(double)(TIMER_STEP-1);
            }
            break;
          default:
            cerr << "something wrong here =_=, in " << __LINE__ << __FILE__ << endl;
            break;
        }

        glBegin( GL_LINES );
        glVertex2dv( coordS.element() );
        glVertex2dv( coordT.element() );
        glEnd();
    }

    // draw vertices
    BGL_FORALL_VERTICES( vertex, *g, Graph ){
        Coord2 coord;
        switch( type ){
          case GEOGRAPHY:
            coord = vertexGeo[ vertex ];
            break;
          case SMOOTH:
            coord = vertexSmooth[ vertex ];
            break;
          case OCTILINEAR:
            glColor4d( 0.6, 0.0, 0.0, 1.0 );
            //coord = vertexCoord[ vertex ];
            if( _time_type == SMOOTH ) {
                coord = vertexTemp[ vertex ] + (double)_time_step*( vertexSmooth[ vertex ] - vertexTemp[ vertex ] )/(double)(TIMER_STEP-1);
            }
            else{
                coord = vertexSmooth[ vertex ] + (double)_time_step*( vertexCoord[ vertex ] - vertexSmooth[ vertex ] )/(double)(TIMER_STEP-1);
            }        
            break;
          default:
            cerr << "something wrong here =_=, in " << __LINE__ << __FILE__ << endl;
            break;
        }

        glBegin( GL_POINTS );
        glVertex2dv( coord.element() );
        glEnd();
    }

    // draw labels
    BGL_FORALL_VERTICES( vertex, *g, Graph ){
        Coord2 coord;
        if( vertexSelectMag[ vertex ] == true && vertexExtstate[ vertex ] == true ) {
            switch( type ){
              case GEOGRAPHY:
                coord = vertexGeo[ vertex ];
                break;
              case SMOOTH:
                coord = vertexSmooth[ vertex ];
                break;
              case OCTILINEAR:
                glColor4d( 0.0, 0.0, 1.0, 1.0 );
                coord = vertexExternal[ vertex ].curSite();
                break;
              default:
                cerr << "something wrong here =_=, in " << __LINE__ << __FILE__ << endl;
                break;
            }

            if( _sim_flag == false ) {

                glBegin( GL_POINTS );
                    glVertex2dv( coord.element() );
                glEnd();
                glBegin( GL_LINES );
                    glVertex2dv( vertexCoord[ vertex ].element() );
                    glVertex2dv( coord.element() );
                glEnd();

#ifdef  SKIP
            Coord2 w( vertexExternal[ vertex ].width(), 0.0 );
            Coord2 h( 0.0, vertexExternal[ vertex ].height() );
            glBegin( GL_LINE_LOOP );
                glVertex2dv( ( coord + w + h ).element() );
                glVertex2dv( ( coord - w + h ).element() );
                glVertex2dv( ( coord - w - h ).element() );
                glVertex2dv( ( coord + w - h ).element() );
            glEnd();
#endif  // SKIP

                glMatrixMode( GL_MODELVIEW );
                glPushMatrix();
                glTranslated( coord.x(), coord.y(), 0.0 );
                glColor4d( 0.0, 0.0, 1.0, 0.2 );
                _paintDisk( vertexExternal[ vertex ].width() );
                glPopMatrix();
            }
        }
    }

#endif  // RENDER

    glDisable( GL_TEXTURE_2D );
    glDisable( GL_ALPHA_TEST );
    glDisable( GL_BLEND );

    //float end = clock()/CLOCKS_PER_SEC + 1.0;
    //while((clock()/CLOCKS_PER_SEC) < end);
}


void Widget::_paintLabels( void )
{
    Graph & g = _metro->g();

    VertexIDMap         vertexID        = get( vertex_myid, g );
    VertexCoordMap      vertexCoord     = get( vertex_mycoord, g );
    VertexSmoothMap     vertexSmooth    = get( vertex_mysmooth, g );
    VertexTempMap       vertexTemp      = get( vertex_mytemp, g );
    VertexSelectMagMap  vertexSelectMag = get( vertex_myselectmag, g );
    VertexScaleMap      vertexScale     = get( vertex_myscale, g );
    VertexExternalMap   vertexExternal  = get( vertex_myexternal, g );
    VertexExtstateMap   vertexExtstate  = get( vertex_myextstate, g );
    VertexTextureIDMap  vertexTextureID = get( vertex_mytextureid, g );


    glLineWidth( 3.0 );
    BGL_FORALL_VERTICES( vd, g, Graph )
    {
        if( vertexExtstate[ vd ] == true ){
        //if( vertexExtstate[ vd ] == true && vertexSelectMag[ vd ] == true ){
            Label  &        label   = vertexExternal[ vd ];
            Coord2          orig;
            if( _time_type == SMOOTH ) {
                orig = vertexTemp[ vd ] + (double)_time_step*( vertexSmooth[ vd ] - vertexTemp[ vd ] )/(double)(TIMER_STEP-1);
            }
            else{
                orig = vertexSmooth[ vd ] + (double)_time_step*( vertexCoord[ vd ] - vertexSmooth[ vd ] )/(double)(TIMER_STEP-1);
            }
            Coord2          dest    = label.curSite();
            Coord2          center  = ( dest-orig )/2.0 + orig;
            double &        width   = label.width();
            double &        height  = label.height();
            double          scale   = vertexScale[ vd ];
            if( vertexSelectMag[ vd ] == false ) {
                dest = orig;
                width = height = 5.0;
            }

            // ---------------------------------------
            // Draw leaders
            if( vertexSelectMag[ vd ] == true ){
                glColor4d( 0.2, 0.2, 0.2, 1.0 );
                glBegin( GL_LINES );
                glVertex2dv( orig.element() );
                glVertex2dv( center.element() );
                glEnd();
            }

            // ---------------------------------------
            // setting
            double ratio = sqrt( 2.0 ) / 2.0;
            //Coord2 lt = dest + Coord2( -ratio*width*scale,  ratio*height*scale );
            //Coord2 lb = dest + Coord2( -ratio*width*scale, -ratio*height*scale );
            //Coord2 rb = dest + Coord2(  ratio*width*scale, -ratio*height*scale );
            //Coord2 rt = dest + Coord2(  ratio*width*scale,  ratio*height*scale );
            Coord2 lt = dest + Coord2( -ratio*width,  ratio*height );
            Coord2 lb = dest + Coord2( -ratio*width, -ratio*height );
            Coord2 rb = dest + Coord2(  ratio*width, -ratio*height );
            Coord2 rt = dest + Coord2(  ratio*width,  ratio*height );
            double l = lt.x();
            double r = rb.x();
            double b = lt.y();
            double t = rb.y();

            // ---------------------------------------
            // Draw bounding boxes
#ifdef  SKIP
            glColor4d( 0.2, 0.2, 0.2, 1.0 );
            glBegin( GL_QUADS );
            glVertex2dv( lt.element() );
            glVertex2dv( rt.element() );
            glVertex2dv( rb.element() );
            glVertex2dv( lb.element() );
            glEnd();
#endif  // SKIP

            // ---------------------------------------
            // Draw texture photos
            glColor4d( 1.0, 0.5, 0.0, 0.8 );
            glBindTexture       ( GL_TEXTURE_2D, vertexTextureID[vd] );
            glTexParameteri     ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
            glTexParameteri     ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
            glTexParameteri     ( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
            glTexParameteri     ( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
            glTexEnvi           ( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

            glEnable( GL_DEPTH_TEST );
            glDepthFunc( GL_ALWAYS );
            glEnable( GL_BLEND );
            glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

            glDisable( GL_LIGHTING );
            glEnable( GL_TEXTURE_2D );


            glBegin( GL_QUADS );
            glTexCoord2d( 0.0, 0.0 );
            glVertex2i( l, b );
            glTexCoord2d( 0.0, 1.0 );
            glVertex2i( l ,t );
            glTexCoord2d( 1.0, 1.0 );
            glVertex2i( r, t );
            glTexCoord2d( 1.0, 0.0 );
            glVertex2i( r, b );
            glEnd();

            glDisable( GL_TEXTURE_2D );
            glDisable( GL_ALPHA_TEST );
            glDisable( GL_BLEND );

#ifdef DEBUG
            glEnable( GL_DEPTH_TEST );
            glDepthFunc( GL_ALWAYS );
            glEnable( GL_BLEND );
            glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

            Coord2 coord = dest;
            glMatrixMode( GL_MODELVIEW );
            glPushMatrix();
            glTranslated( coord.x(), coord.y(), 0.0 );
            glColor4d( 0.0, 0.0, 1.0, 0.2 );
            _paintDisk( vertexExternal[ vd ].width() );
            glPopMatrix();

            glDisable( GL_ALPHA_TEST );
            glDisable( GL_BLEND );
#endif  // DEBUG

#ifdef  SKIP        
            // ---------------------------------------
            // Draw text
            unsigned int imageC = getVertexCommentWidth( g, vd );

            int frameW = abs( r - l );
            //glColor4d( 1.0, 1.0, 1.0, 1.0 );
            //glColor4d( 0.2, 0.2, 0.2, 1.0 );
            //glColor4d( 0.8, 1.0, 1.0, 1.0 );    // <-- default
            //glColor4d( 0.8, 1.0, 0.8, 1.0 );
            double scale = min( 1.5, ( double )( frameW/(5.0*min_grid_size) ) );

            if ( imageC > 8 ) scale *= 8.0/( double )imageC;

            // cerr << " scale = " << scale << endl;
            // int nameH = ptrV->height();
            double xpos = 0.50 * l + 0.50 * r;
            double ypos = 0.70 * t + 0.30 * b;

            drawString( xpos, ypos, vertexComment[vd].c_str(), scale, 0.0, CENTERING );
            //drawName( vd, 0.0 );

            // ---------------------------------------
            // Draw leaders
            glEnable( GL_DEPTH_TEST );
            glDepthFunc( GL_LESS );
            glEnable( GL_LINE_SMOOTH );
            glEnable( GL_BLEND );
            glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
            glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

            glLineWidth( 4.0 );
            //glLineStipple( 4, 0xAAAA );
            //glEnable( GL_LINE_STIPPLE );
            glColor4d( 0.6, 0.6, 0.6, 1.0 );
            glBegin( GL_LINES );
            glVertex2dv( orig.element() );
            glVertex2dv( joint.element() );
            glEnd();
            glBegin( GL_LINES );
            glVertex2dv( joint.element() );
            glVertex2dv( dest.element() );
            glEnd();

            //glLineStipple( 1, 0xFFFF );
            glDisable( GL_LINE_SMOOTH );
            glDisable( GL_BLEND );
            glDisable( GL_DEPTH_TEST );
#endif  // SKIP
        }
    }
}


void Widget::_paintName( void )
{
    Graph * g = & _metro->g();
    if( _sim_flag == true ) g = & _simmetro->g();
    else g = & _metro->g();

    double shift = 5;

    VertexIndexMap          vertexIndex         = get( vertex_index, *g );
    VertexIDMap             vertexID            = get( vertex_myid, *g );
    VertexNameMap           vertexName          = get( vertex_myname, *g );
    VertexCoordMap          vertexCoord         = get( vertex_mycoord, *g );
    VertexSmoothMap         vertexSmooth        = get( vertex_mysmooth, *g );
    VertexGeoMap            vertexGeo           = get( vertex_mygeo, *g );
    VertexSelectMagMap      vertexSelectMag     = get( vertex_myselectmag, *g );
    VertexSelectTopMap      vertexSelectTop     = get( vertex_myselecttop, *g );
    VertexScaleMap          vertexScale         = get( vertex_myscale, *g );
    VertexGeodesicMap       vertexGeodesic      = get( vertex_mygeodesic, *g );
    VertexZoneMap           vertexZone          = get( vertex_myzone, *g );
    VertexExternalMap       vertexExternal      = get( vertex_myexternal, *g );
    EdgeWeightMap           edgeWeight          = get( edge_weight, *g );
    EdgeIDMap               edgeID              = get( edge_myid, *g );
    EdgeGeoAngleMap         edgeGeoAngle        = get( edge_mygeoangle, *g );

//#ifdef  SKIP
    // draw edge weight
    glColor4d( 0.0, 0.0, 1.0, 1.0 );
    BGL_FORALL_EDGES( edge, *g, Graph ){
        VertexDescriptor vdS = source( edge, *g );
        VertexDescriptor vdT = target( edge, *g );
        Coord2 coord = ( vertexCoord[ vdS ] + vertexCoord[ vdT ] ) / 2.0;
        string idstr = to_string( edgeID[ edge ] );
        //string idstr = to_string( edgeWeight[ edge ] );
        //string idstr = to_string( edgeGeoAngle[ edge ] );
        renderText( coord.x() + shift, coord.y() + shift, 0.0, idstr.c_str(),
                    QFont( "Arial", 12, QFont::Bold, false ) );
    }
//#endif  // SKIP

    // draw vertex name
#ifdef  SKIP
    glColor4d( 0.6, 0.8, 0.4, 1.0 );
    BGL_FORALL_VERTICES( vertex, *g, Graph ){

        string idstr = to_string( vertexID[ vertex ] );
        Coord2 coord = vertexGeo[ vertex ];
        renderText( coord.x() + shift, coord.y() + shift, 0.0, idstr.c_str(),
                    QFont( "Arial", 12, QFont::Bold, false ) );
    }

    glColor4d( 0.0, 0.0, 0.6, 1.0 );
    BGL_FORALL_VERTICES( vertex, *g, Graph ){

        string idstr = to_string( vertexID[ vertex ] );
        Coord2 coord = vertexSmooth[ vertex ];
        renderText( coord.x() + shift, coord.y() + shift, 0.0, idstr.c_str(),
                    QFont( "Arial", 12, QFont::Bold, false ) );
    }

    glColor4d( 0.6, 0.0, 0.0, 1.0 );
    BGL_FORALL_VERTICES( vertex, *g, Graph ){

        string idstr = to_string( vertexID[ vertex ] );
        Coord2 coord = vertexCoord[ vertex ];
        renderText( coord.x() + shift, coord.y() + shift, 0.0, idstr.c_str(),
                    QFont( "Arial", 12, QFont::Bold, false ) );
    }

#endif  // SKIP
    glColor4d( 0.6, 0.0, 0.0, 1.0 );
    BGL_FORALL_VERTICES( vertex, *g, Graph ){

        //string idstr = to_string( vertexSelectMag[ vertex ] );
        //string idstr = to_string( vertexName[ vertex ] );
        //string idstr = to_string( vertexScale[ vertex ] );
        //string idstr = to_string( vertexExternal[ vertex ].leaderWeight() );
        string idstr = to_string( vertexID[ vertex ] );
        //string idstr = to_string( vertexIndex[ vertex ] );
        //string idstr = to_string( vertexSelectTop[ vertex ] );
        //string idstr = to_string( vertexZone[ vertex ] );
        //string idstr = to_string( vertexGeodesic[ vertex ] );
        Coord2 coord = vertexCoord[ vertex ];
        renderText( coord.x() + shift, coord.y() + shift, 0.0, idstr.c_str(),
                    QFont( "Arial", 12, QFont::Bold, false ) );
    }

}


void Widget::_paintMagnification( void )
{
//    double radius = magnified_radius;
//    // double radius = MAGNIFICATION_RADIUS;
//    glColor4d( 0.3, 0.3, 0.3, 0.1 );
//
//    glEnable( GL_LINE_SMOOTH );
//    glEnable( GL_BLEND );
//    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
//    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
//
//    glMatrixMode( GL_MODELVIEW );
//    glPushMatrix();
//    glTranslated( _now_pointer_x, _now_pointer_y, 0.0 );
//    //cerr << "_now_x = " << _now_pointer_x << " _now_y = " << _now_pointer_y << endl;
//#ifdef  SKIP
//    glBegin( GL_LINE_LOOP );
//    glVertex2d( -radius/2, -radius/2 );
//    glVertex2d( -radius/2,  radius/2 );
//    glVertex2d(  radius/2,  radius/2 );
//    glVertex2d(  radius/2, -radius/2 );
//    glEnd();
//#endif  // SKIP
//    //_paintDisk( radius );
//    glPopMatrix();

#ifdef  SKIP
    // try Fermat's spiral
    int num = rng_samples;
    double golden_angle = M_PI * ( 3.0 - sqrt(5.0) );

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glTranslated( _now_pointer_x, _now_pointer_y, 0.0 );

    glBegin( GL_POINTS );
    for( int i = 0; i < num; i++ ){
        double theta = i * golden_angle;
        double r = radius * sqrt( i ) / sqrt( num );
        Coord2 coord(  r * cos( theta ), r * sin( theta ) );
        glVertex2dv( coord.element() );
    } 
    glEnd();
    _paintDisk( radius );
    glPopMatrix();

    glDisable( GL_TEXTURE_2D );
    glDisable( GL_ALPHA_TEST );
    glDisable( GL_BLEND );
#endif  // SKIP
}


void Widget::initTextures( const char * imagename )
{
    ostringstream ostr;
    Graph & g = _metro->g();

    VertexIndexMap vertexIndex = get( vertex_index, g );
    VertexTexNameMap vertexTexName = get( vertex_mytexname, g );
    VertexTextureIDMap vertexTextureID = get( vertex_mytextureid, g );
    VertexExternalMap vertexExternal = get( vertex_myexternal, g );
    VertexCoordMap vertexCoord = get( vertex_mycoord, g );

    static unsigned int     nTextures = 1;
    unsigned int *          textureID = NULL;

    // Parameter settings for the texture
    glPixelStorei       ( GL_UNPACK_ALIGNMENT, 4 );

#define INCREMENTAL_TEXTURE_ALLOCATION
#ifndef INCREMENTAL_TEXTURE_ALLOCATION
    // This may damage existing texture memory while I am not quite sure..
    for ( unsigned int k = 0; k < _metro->nStations(); ++k ) {
        if ( glIsTexture( k ) )
            glDeleteTextures( (GLsizei)1, &k );
    }
#endif  // INCREMENTAL_TEXTURE_ALLOCATION
    if ( textureID != NULL ) {
        cerr << "Texture IDs initialized." << endl;
        glDeleteTextures( (GLsizei)nTextures, textureID );
        textureID = NULL;
    }

    nTextures = 1 + _metro->nStations();
    cerr << " Number of thumbnail photos = " << nTextures << endl;
    textureID = new unsigned int [ nTextures ];
    glGenTextures( nTextures, textureID );

    unsigned int curID = 1;
    BGL_FORALL_VERTICES( vertex, g, Graph )
    {
        ostr.clear();
        ostr.str("");

        if( vertexIndex[ vertex ] >= _metro->nStations() ) break;
        ostr << imagename << vertexTexName[ vertex ] << "_1.png" << ends;
     //ostr << imagename << "modern11.png" << ends;

        cerr << "==============================" << endl;
        IplImage * img = cvLoadImage( ostr.str().c_str() );
        if ( !img ) { // If the input image file cannot be found, just exit.
            cerr << "############################################################" << endl;
            cerr << "Cannot load the file: " << ostr.str().c_str() << endl;
            cerr << "############################################################" << endl;
            return;
        }
        else {
            cerr << "Loading the thumbnail file : " << ostr.str().c_str() << endl;
        }

        // Just transform the image into its mirror image (upside down)
        if ( img->origin != 0 ) {
            cvFlip( img, img );
        }

        cerr << " Veretex ID = " << vertexIndex[ vertex ] << " cur texture ID = " << textureID[ curID ] << endl;
        glBindTexture( GL_TEXTURE_2D, textureID[ curID ] );
        vertexTextureID[ vertex ] = textureID[ curID ];
        //vertexExternal[ vertex ].width() = img->width;
        //vertexExternal[ vertex ].height() = img->height;
        vertexExternal[ vertex ].width() = 5.0;
        vertexExternal[ vertex ].height() = 5.0;

#ifdef DEBUG
        cerr << " k = " << vertexIndex[ vertex ] << endl
             << " Texture ID = " << textureID[ vertexIndex[ vertex ]  ]
             << " Coord = " << vertexCoord[ vertex ].x() << ", " << vertexCoord[ vertex ].y()
             << " size : " << vertexExternal[ vertex ].width() << " x " << vertexExternal[ vertex ].height() << endl;
#endif  // DEBUG

        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,
                      img->width, img->height,
                      0, GL_BGR_EXT, GL_UNSIGNED_BYTE, img->imageData );


        //const int block = DEFAULT_TEXTURE_BLOCK_SIZE;
        //vertexExternal[ vertex ].setFrame( Grid2( (int)(img->width/block), (int)(img->height/block) ) );

        cerr << "imgWidth= " << img->width << " imgHeight= " << img->height << endl;


#ifdef DEBUG
        cerr << "Image label box: width = " << vertexExternal[ vertex ].width()
             << " height = " << vertexExternal[ vertex ].height() << endl;
#endif  // DEBUG

        cvReleaseImage( &img );
        curID++;

    }

    delete [] textureID;
    textureID = NULL;

}


void Widget::_plotStations( void )
{
    Graph               & g             = _metro->g();
    VertexIDMap         vertexID        = get( vertex_myid, g );
    VertexGeoMap        vertexGeo       = get( vertex_mygeo, g );
    //VertexCoordMap      vertexCoord     = get( vertex_mycoord, g );

    // for Picking
    glInitNames();
    glPushName( -1 );

    glPointSize( 1.0 );
    BGL_FORALL_VERTICES( vd, g, Graph ) {
        unsigned int    idV = vertexID[ vd ];
        glLoadName( ( int )idV );
        Coord2 &        coord = vertexGeo[ vd ];
        //Coord2 &        coord = vertexCoord[ vd ];
        glBegin( GL_POINTS );
        glVertex2dv( coord.element() );
        glEnd();
        glLoadName( -1 );
    }

    glPopName();
    glInitNames();
}


void Widget::_plotLines( void )
{
    Graph               & g             = _metro->g();
    //VertexIDMap         vertexID        = get( vertex_myid, g );
    VertexCoordMap      vertexCoord     = get( vertex_mycoord, g );
    EdgeIDMap           edgeID          = get( edge_myid, g );

    // for Picking
    glInitNames();
    glPushName( -1 );

    glPointSize( 1.0 );
    BGL_FORALL_EDGES( ed, g, Graph ) {
        unsigned int    idE = edgeID[ ed ];
        VertexDescriptor vdS = source( ed, g );
        VertexDescriptor vdT = target( ed, g );
        glLoadName( ( int )idE );
        Coord2 &        coordS = vertexCoord[ vdS ];
        Coord2 &        coordT = vertexCoord[ vdT ];
        glBegin( GL_LINES );
        glVertex2dv( coordS.element() );
        glVertex2dv( coordT.element() );
        glEnd();
        glLoadName( -1 );
    }

    glPopName();
    glInitNames();
}


bool Widget::_handleVertex( int nHits, unsigned int * buf, const int button, const int modifier )
{
    unsigned int *      ptr             = NULL; //, names;
    float               minDepth        = 1000.0;
    int                 hitID           = NO_INDEX;
    Graph &             g               = _metro->g();

    VertexIDMap             vertexID            = get( vertex_myid, g );
    VertexSelectMagMap      vertexSelectMag     = get( vertex_myselectmag, g );

    ptr = buf;

    for ( int i = 0; i < nHits; ++i ) { // for each bit
        if ( ptr[ 0 ] != 1 ) {
            cerr << " Number of names for hit = " << ( int )ptr[ 0 ] << endl;
            assert( ptr[ 0 ] == 1 );
        }
        float curDepth = (float)ptr[ 1 ]/0xffffffff;
        int curID = ( int )ptr[ 3 ];
#ifdef DEBUG
        cerr << " i = " << i
             << " [0]: " << ptr[ 0 ]
             << " [1]: " << ptr[ 1 ]
             << " [2]: " << ptr[ 2 ]
             << " [3]: " << ptr[ 3 ] << endl;
#endif  // DEBUG
        if ( ( curDepth < minDepth ) && ( curID != NO_INDEX ) ) {
            minDepth = curDepth;
            hitID = ptr[ 3 ];
        }
        ptr += 4;
    }

    cerr << " hitVID = " << hitID << " depth = " << minDepth << endl;

    if ( hitID != NO_INDEX ) {

        BGL_FORALL_VERTICES( vertex, g, Graph ){

            // shortest path
            if( button == Qt::LeftButton && modifier == Qt::ControlModifier ) {
                if( vertexID[ vertex ] == hitID ) sourceVD = vertex;
            }
            else if( button == Qt::MiddleButton && modifier == Qt::ControlModifier ) {
                if( vertexID[ vertex ] == hitID ) targetVD = vertex;
            }

            // semantic zoom
            if( modifier == Qt::ShiftModifier ) {
                if( vertexID[ vertex ] == hitID ){
                    vertexSelectMag[ vertex ] = true;
                    focusVD = vertex;
                }
                else{
                    vertexSelectMag[ vertex ] = false;
                }
            }
        }
        //vd = vertex( ( unsigned int )hitID, g );
        return true;
    }

    return false;
}


bool Widget::_handleEdge( int nHits, unsigned int * buf )
{
    unsigned int *      ptr             = NULL; //, names;
    float               minDepth        = 1000.0;
    int                 hitID           = NO_INDEX;
    Graph &             g               = _metro->g();

    EdgeIDMap           edgeID          = get( edge_myid, g );
    EdgeSelectShiftMap  edgeSelectShift = get( edge_myselectshift, g );
    EdgeSelectCtrlMap   edgeSelectCtrl  = get( edge_myselectctrl, g );

    ptr = buf;

    for ( int i = 0; i < nHits; ++i ) { // for each bit
        if ( ptr[ 0 ] != 1 ) {
            cerr << " Number of names for hit = " << ( int )ptr[ 0 ] << endl;
            assert( ptr[ 0 ] == 1 );
        }
        float curDepth = (float)ptr[ 1 ]/0xffffffff;
        int curID = ( int )ptr[ 3 ];
#ifdef DEBUG
        cerr << " i = " << i
             << " [0]: " << ptr[ 0 ]
             << " [1]: " << ptr[ 1 ]
             << " [2]: " << ptr[ 2 ]
             << " [3]: " << ptr[ 3 ] << endl;
#endif  // DEBUG
        if ( ( curDepth < minDepth ) && ( curID != NO_INDEX ) ) {
            minDepth = curDepth;
            hitID = ptr[ 3 ];
        }
        ptr += 4;
    }

    cerr << " hitEID = " << hitID << " depth = " << minDepth << endl;

    if ( hitID != NO_INDEX ) {
        bool found = false;
        BGL_FORALL_EDGES( edge, g, Graph ){
            if ( edgeID[ edge ] == hitID ){
                edgeSelectShift[ edge ] = true;
            }
            else{
                edgeSelectShift[ edge ] = false;
            }
        }
        return true;
    }

    return false;
}

void retrieveCluster( int nHits, unsigned int * buf, vector< unsigned int > & ids )
{
    unsigned int * ptr = NULL; //, names;
    //float minDepth = 1000.0;
    //int hitID = NO_INDEX;

    cerr << "**** retrieveCluster ****" << endl;
    ids.clear();

    ptr = buf;

    for ( int i = 0; i < nHits; ++i ) { // for each bit
        if ( ptr[ 0 ] != 1 ) {
            cerr << " Number of names for hit = " << ( int )ptr[ 0 ] << endl;
            assert( ptr[ 0 ] == 1 );
        }
        //float curDepth = (float)ptr[ 1 ]/0xffffffff;
        int curID = ( int )ptr[ 3 ];
#ifdef DEBUG
        cerr << " i = " << i
             << " [0]: " << ptr[ 0 ]
             << " [1]: " << ptr[ 1 ]
             << " [2]: " << ptr[ 2 ]
             << " [3]: " << ptr[ 3 ] << endl;
#endif  // DEBUG
        if ( curID != NO_INDEX ) {
            vector< unsigned int >::iterator it = find( ids.begin(), ids.end(), curID );
            if ( it == ids.end() ) ids.push_back( curID );
        }
        ptr += 4;
    }

#ifdef  DEBUG
    for ( unsigned int k = 0; k < ids.size(); ++k ) {
        cerr << " set[ " << setw( 2 ) << k << " ] = " << setw( 3 ) << ids[ k ];
        if ( k % 2 == 1 ) cerr << endl;
    }
    if ( ids.size() % 2 == 1 ) cerr << endl;
#endif  // DEBUG

    glutPostRedisplay();
}

void Widget::selectMagnification( double x, double y, int button, int modifier )
{
//    Graph &                 g                   = _metro->g();
//    VertexIDMap             vertexID            = get( vertex_myid, g );
//    VertexGeoMap            vertexGeo           = get( vertex_mygeo, g );
//    VertexNameMap           vertexName          = get( vertex_myname, g );
//    VertexSelectMagMap      vertexSelectMag     = get( vertex_myselectmag, g );
//
//    unsigned int selectBuf[ BUFFER_SIZE ];
//    int nHits;
//    int viewport[ 4 ];
//
//    cerr << " pickVertexSet " << " x = " << x << " y = " << y << endl;
//
//    glGetIntegerv( GL_VIEWPORT, viewport );
//
//    // Picking begins here
//    glSelectBuffer( BUFFER_SIZE, selectBuf );
//    glRenderMode( GL_SELECT );
//
//    glInitNames();
//
//    glMatrixMode( GL_PROJECTION );
//    glPushMatrix(); // <====
//    glLoadIdentity();
//    // create small picking region near cursor location
//    //const double tolerance = 10.0;
//    const double tolerance = 2.0 * magnified_radius;
//    //const double tolerance = 2.0 * MAGNIFICATION_RADIUS;
//#ifdef DEBUG
//    cerr << "viewport = " << viewport[0] << "," << viewport[1] << ","
//         << viewport[2] << "," << viewport[3] << endl;
//    cerr << " x = " << x << " y = " << y << endl;
//#endif  // DEBUG
//    gluPickMatrix( (double)x, (double)(viewport[3]-y), tolerance, tolerance, viewport );
//
//    glOrtho( -width()/2, width()/2, -height()/2, height()/2, -2.0, 2.0 );
//
//    glMatrixMode( GL_MODELVIEW );
//    glPushMatrix();     // <====
//    glLoadIdentity();
//
//    _plotStations();
//
//    glMatrixMode( GL_PROJECTION );
//    glPopMatrix();      // <====
//    glMatrixMode( GL_MODELVIEW );
//    glPopMatrix();      // <====
//
//    glFlush();
//
//    vector< unsigned int > ids;
//    nHits = glRenderMode( GL_RENDER );
//    retrieveCluster( nHits, selectBuf, ids );
//
//    // cerr << "idSize = " << ids.size() << endl;
//
//    // clear selection
//    BGL_FORALL_VERTICES( vd, g, Graph ) {
//        vertexSelectMag[ vd ] = false;
//    }
//    for ( unsigned int k = 0; k < ids.size(); k++ ) {
//        BGL_FORALL_VERTICES( vd, g, Graph ) {
//            if ( vertexID[ vd ] == ids[ k ] ) {
//                Coord2 coord( x - width()/2, -y + height()/2 );
//                double dist = ( vertexGeo[ vd ] - coord ).norm();
//                //double dist = ( vertexCoord[ vd ] - coord ).norm();
//                // cerr << "coord = " << coord << "vertex = " << vertexCoord[vd] << endl;
//                if( dist < magnified_radius ) {
//                    vertexSelectMag[ vd ] = true;
//                    cerr << "pickID = " << vertexID[ vd ] << ", name = " << vertexName[ vd ] << endl;
//#ifdef  SKIP
//                    OutEdgeIterator e, e_end;
//                    for ( tie( e, e_end ) = out_edges( vd, g ); e != e_end; ++e ) {
//                        EdgeDescriptor ed = *e;
//                        VertexDescriptor vT = target( ed, g );
//                        vertexSelectMag[ vT ] = true;
//                    }
//#endif  // SKIP
//                }
//            }
//        }
//    }
//    // Picking ends here
//    _isDeformed = false;
//    glutPostRedisplay();
}

void Widget::selectTopology( void )
{
    Graph &                 g                   = _metro->g();
    VertexIDMap             vertexID            = get( vertex_myid, g );
    VertexCoordMap          vertexCoord         = get( vertex_mycoord, g );
    VertexSelectMagMap      vertexSelectMag     = get( vertex_myselectmag, g );
    VertexSelectTopMap      vertexSelectTop     = get( vertex_myselecttop, g );

    BGL_FORALL_VERTICES( vd, g, Graph ) {
        vertexSelectTop[ vd ] = false;
    }
    BGL_FORALL_VERTICES( vd, g, Graph ) {
        if( vertexSelectMag[ vd ] == true ) {
            focusVD = vd;
            calcGeodesicDistance();
        }
    }
}


void Widget::setMetroWeight( void )
{
    Graph &             g               = _metro->g();
    double              meanNodeSize    = _metro->meanVSize();
    double              length          = _metro->dBeta();
    setGraphWeight( _last_pointer_x, _last_pointer_y, meanNodeSize, length, g );
}

void Widget::setSimMetroWeight( void )
{
    Graph &             g               = _simmetro->g();
    double              meanNodeSize    = _simmetro->meanVSize();
    double              length          = _metro->dBeta();
    setGraphWeight( _last_pointer_x, _last_pointer_y, meanNodeSize, length, g );
}

void Widget::pickVertex( double x, double y, int button, int modifier )
{
    Graph &             g               = _metro->g();
    VertexIDMap         vertexID        = get( vertex_myid, g );
    VertexCoordMap      vertexCoord     = get( vertex_mycoord, g );

    unsigned int selectBuf[ BUFFER_SIZE ];
    int nHits;
    int viewport[ 4 ];

    //cerr << " pickVertex " << endl;

    glGetIntegerv( GL_VIEWPORT, viewport );

    // Picking begins here
    glSelectBuffer( BUFFER_SIZE, selectBuf );
    glRenderMode( GL_SELECT );

    glInitNames();

    glMatrixMode( GL_PROJECTION );
    glPushMatrix(); // <====
    glLoadIdentity();
    // create small picking region near cursor location
    //double const size tolerance = 1.0;
    //double const size tolerance = 5.0;
    //const double tolerance = 10.0;
    const double tolerance = 20.0;
#ifdef DEBUG
    cerr << "viewport = " << viewport[0] << "," << viewport[1] << ","
         << viewport[2] << "," << viewport[3] << endl;
    cerr << " x = " << x << " y = " << y << endl;
#endif  // DEBUG
    gluPickMatrix( (double)x, (double)(viewport[3]-y), tolerance, tolerance, viewport );

    glOrtho( -width()/2, width()/2, -height()/2, height()/2, -2.0, 2.0 );

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();     // <====
    glLoadIdentity();

    _plotStations();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();      // <====
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();      // <====

    glFlush();

    nHits = glRenderMode( GL_RENDER );

    switch ( button ) {
      case Qt::LeftButton:
          _handleVertex( nHits, selectBuf, button, modifier );
          break;
      case Qt::MiddleButton:
          _handleVertex( nHits, selectBuf, button, modifier );
          break;
      default:
          // do nothing
          break;
    }
    // Picking ends here

    glutPostRedisplay();

}

void Widget::pickEdge( double x, double y, int button, int modifier )
{
    Graph &             g               = _metro->g();
    //VertexIDMap         vertexID        = get( vertex_myid, g );
    //VertexCoordMap      vertexCoord     = get( vertex_mycoord, g );

    unsigned int selectBuf[ BUFFER_SIZE ];
    int nHits;
    int viewport[ 4 ];

    //cerr << " pickEdge " << endl;

    glGetIntegerv( GL_VIEWPORT, viewport );

    // Picking begins here
    glSelectBuffer( BUFFER_SIZE, selectBuf );
    glRenderMode( GL_SELECT );

    glInitNames();

    glMatrixMode( GL_PROJECTION );
    glPushMatrix(); // <====
    glLoadIdentity();
    // create small picking region near cursor location
    //double const size tolerance = 1.0;
    //double const size tolerance = 5.0;
    //const double tolerance = 10.0;
    const double tolerance = 20.0;
#ifdef DEBUG
    cerr << "viewport = " << viewport[0] << "," << viewport[1] << ","
         << viewport[2] << "," << viewport[3] << endl;
    cerr << " x = " << x << " y = " << y << endl;
#endif  // DEBUG
    gluPickMatrix( (double)x, (double)(viewport[3]-y), tolerance, tolerance, viewport );

    glOrtho( -width()/2, width()/2, -height()/2, height()/2, -2.0, 2.0 );

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();     // <====
    glLoadIdentity();

    _plotLines();

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();      // <====
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();      // <====


    glFlush();

    nHits = glRenderMode( GL_RENDER );

    switch ( button ) {
      case Qt::LeftButton:
          _handleEdge( nHits, selectBuf );
          break;
      default:
          // do nothing
          break;
    }
    // Picking ends here

    glutPostRedisplay();

}

bool Widget::calcShortestPath( void )
{
    Graph & g = _metro->g();
    VertexIDMap     vertexID        = get( vertex_myid, g );

    if( sourceVD != NULL && targetVD != NULL ){
        cerr << "Calculating shortest path between ID(" << vertexID[ sourceVD ] 
             << ") and ID(" << vertexID[ targetVD ] << ")" << endl;
        shortestPath( sourceVD, targetVD, g );
        return true;
    }
    
    return false;
}

bool Widget::calcGeodesicDistance( void )
{
    Graph & g = _metro->g();
    VertexIDMap     vertexID        = get( vertex_myid, g );

    if( focusVD != NULL ){
        // cerr << "Calculating geodesic distance of ID(" << vertexID[ focusVD ] << ")" << endl;
        geodesicDistance( focusVD, g );
        return true;
    }

    return false;
}
