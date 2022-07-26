/*
*	Stationmatching.cpp
*	@author Tobias Batik
*	@version 1.0  14/03/2021
*
*   Creating matching, and transforming guide
*   based on matching
*/

#include "Stationmatching.h"

void Stationmatching::init( void ) 
{
    UndirectedGraph &gGuide = _guide->g(); 
    _nVertices = 0; 
    BGL_FORALL_VERTICES( vertex, gGuide, UndirectedGraph) {
        _nVertices++; 
    }
    _nDifConstrs = 2; // translating x,y
    _nDifConstrs += 2; // Rotating x,y
    _nDifConstrs += 2; // scaling x,y  
    _nVars = 2 * _nVertices;
    _nConstrs = _nDifConstrs * _nVertices;
    
    _wTransX = _wTransY = 1.0;
    _wRotate = 0.0;
    _wScale = 1.0; 

    _initVars();
    _initCoefs(); 
    _initOutputs();
    _updateOutputs();
}

void Stationmatching::prepare( Metro * __metro, Guide * __guide,  double __half_width, double __half_height ) 
{
    _metro = __metro; 
    _guide = __guide; 
    _half_width = __half_width; 
    _half_height = __half_height; 

    _nVars = _nConstrs = 0; 
}

/*
*   input   none
*   output  none
*   
*   Distort Metro Map not uniformly based on matching 
*/
void Stationmatching::distortMap( void ) 
{
    cout << "distorting Map" << endl; 
    if (_stationDisc.size() <= 2) return; 
    
    // Find shortest Guide
    Coord2 minS; 
    Coord2 minG;
    double minDist = 1e+10;
    UndirectedGraph &gMetro = _metro->g();
    UndirectedGraph &gGuide = _guide->g(); 

    for ( int i = 0; i < (int) _stationDisc.size(); i++ ) {
        Coord2 s = *gMetro[_stationDisc[i]].smoothPtr; 
        Coord2 g = *gGuide[_guideDisc[i]].coordPtr; 
        Coord2 delta = s - g; 
        double d = magnitude(delta);
        if ( d < minDist ) {
            minS = s; 
            minG = g; 
            minDist = d; 
        }
    }
    Coord2 cp = minS;
    
    // Get average x,y scaling value
    double sSumX = 0; 
    double sSumY = 0; 
    double dSumX = 0; 
    double dSumY = 0; 
    for ( int i = 0; i < (int)_stationDisc.size(); i++ ) {
        Coord2 gp = *gGuide[_guideDisc[i]].coordPtr; 
        double dgx = abs(cp.x() - gp.x()); 
        double dgy = abs(cp.y() - gp.y()); 

        Coord2 sp = *gMetro[_stationDisc[i]].smoothPtr;
        double dsx = abs(cp.x() - sp.x()); 
        double dsy = abs(cp.y() - sp.y());

        if (dsx != 0 && dsx != 0) {
            double sx = dgx / dsx; 
            double sy = dgy / dsy; 

            sSumX += dsx * sx; 
            sSumY += dsy * sy; 
            dSumX += dsx; 
            dSumY += dsy;
        }
    }
    double scaleX = sSumX / dSumX; 
    double scaleY = sSumY / dSumY; 

    // Distort Metro Map
    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        Coord2 loc = *gMetro[ vertex ].smoothPtr;
        loc = loc - cp; 
        loc = Coord2(loc.x() * scaleX, loc.y() * scaleY);
        loc = loc + cp;  
        gMetro[vertex].smoothPtr->x() = loc.x();
        gMetro[vertex].smoothPtr->y() = loc.y();

        gMetro[vertex].coordPtr->x() = loc.x();
        gMetro[vertex].coordPtr->y() = loc.y();
    }
}

/*
*   Input   none
*   Output  none
*
*   calculates the transformation iterative 
*/
void Stationmatching::conjugateGradient( void ) 
{
    Eigen::MatrixXd A(_coef.rows(), _coef.cols()); 
    Eigen::VectorXd b; 
    Eigen::VectorXd Ap; 
    A = _coef.transpose() * _coef; 
    b = _coef.transpose() * _output;; 

    Eigen::VectorXd err = b - A * _var; 
    Eigen::VectorXd p = err; 
    double rsold = err.transpose() * err; 
    
    int numRounds = 1; 
    for( int i = 0; i < (int) _var.size(); i++ ) {  
        b = _coef.transpose() * _output;; 
        
        Ap = A * p; 
        double alpha = rsold / (double) (p.transpose() * Ap); 
        err = b - A * _var; // has been change otherwise always err = 0; 
        _var = _var + alpha * p;

        double rsnew = err.transpose() * err; 
        if ( sqrt( rsnew ) < 1e-10 ) {   
            cerr << "sqrterror(" << i << ") = " << sqrt( err.adjoint() * err ) << endl;
            break;
        }
        
        p = err + (rsnew / rsold)  * p;
        rsold = rsnew; 
        
        retrieve(); 
        _updateOutputs(); 
        numRounds++; 
    }
}

/*
*   Input   none
*   Output  none
*
*   retrieve new guid positions
*/
void Stationmatching::retrieve( void ) 
{
    UndirectedGraph &gGuide = _guide->g(); 
    int nRows = 0; 
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        double x = _var(nRows * 2 + 0); 
        double y = _var(nRows * 2 + 1);  
        gGuide[ vertex ].coordPtr->x() = x; 
        gGuide[ vertex ].coordPtr->y() = y;  
        nRows++;
    }
}

/*
*   Input   none
*   Output  none
*
*   calculate new Output 
*/
void Stationmatching::_updateOutputs( void ) 
{
    UndirectedGraph &gMetro = _metro->g();
    UndirectedGraph &gGuide = _guide->g(); 

    // Find average translation
    double transX = 0; 
    double transY = 0;  
    for (int i = 0; i < (int) _stationDisc.size(); i++ ) {
        Coord2 s = *gMetro[ _stationDisc[i] ].smoothPtr; 
        Coord2 n = *gGuide[ _guideDisc[i] ].coordPtr;
        double dx = s.x() - n.x(); 
        if (isnan(dx)) cout << "dx is nan" << endl; 
        double dy = s.y() - n.y();  
        transX += dx; 
        transY += dy; 
    }
    transX /= _stationDisc.size(); 
    transY /= _stationDisc.size(); 

    //Udpate output translate
    int nRows = 0;
    Coord2 center = getCenterOfGuideNodes();
    
    // Translate 
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        Coord2 gc = *gGuide[vertex].coordPtr;
        _output(nRows) = (double) ( gc.x() + transX ) * _wTransX ;
        nRows++;  
        _output(nRows) = (double) (gc.y() + transY) * _wTransY; 
        nRows++; 
    }

    // Calculate Average Angle
    double a_pos = 0;
    double a_neg = 0; 
    int num_pos = 0; 
    int num_neg = 0;  
    for (int i = 1; i < (int) _stationDisc.size(); i++ ){
        Coord2 s0 = *gMetro[ _stationDisc[i - 1] ].smoothPtr; 
        Coord2 n0 = *gGuide[ _guideDisc[i - 1 ] ].coordPtr; 
        Coord2 s1 = *gMetro[ _stationDisc[i] ].smoothPtr; 
        Coord2 n1 = *gGuide[ _guideDisc[i] ].coordPtr; 
        Coord2 s = s1 - s0; 
        Coord2 n = n1 - n0;
        if ( magnitude(s) != 0 && magnitude(n) != 0) {
            double dns = n.x() * s.x() + n.y() * s.y(); 
            double ai = acos( dns / ( magnitude(n) * magnitude(s)) );      

            double sig_i = n[0]*s[1] - n[1]*s[0];
            if (sig_i > 0) {
                a_pos += ai; 
                num_pos++; 
            } else {
                a_neg -= ai;
                num_neg++;  
            }
        }
    }
    if (num_pos > 0) a_pos /= num_pos; 
    if (num_neg > 0) a_neg /= num_neg; 
    double a = 0; 
    if (num_pos + num_neg != 0) a = a_pos * num_pos / (num_pos + num_neg) + a_neg * num_neg / (num_pos + num_neg); 

    // Update output rotate
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        Coord2 gc = *gGuide[vertex].coordPtr;
        Coord2 gd = gc - center; 
        _output(nRows) =(double) ( ( gd.x() * cos(a) - gd.y() * sin(a) ) + center.x()) * _wRotate;
        nRows++; 
        _output(nRows) = (double) ( ( gd.x() * sin(a) + gd.y() * cos(a) ) + center.y()) * _wRotate;
        nRows++;    
    } 

    // Scaling 
    double scale = 1.0; 
    if (_stationDisc.size() > 1 ) { 
        double scaleSum = 0; 
        double sMsum = 0; 
        for (int i = 1; i < (int) _stationDisc.size(); i++) {
            Coord2 g1 = *gGuide[_guideDisc[i]].coordPtr; 
            Coord2 m1 = *gMetro[_stationDisc[i]].smoothPtr; 
            Coord2 g0 = *gGuide[_guideDisc[i - 1]].coordPtr; 
            Coord2 m0 = *gMetro[_stationDisc[i - 1]].smoothPtr; 
            
            Coord2 dG = g1 - g0; 
            Coord2 dM = m1 - m0;
            double sM = magnitude( dG); 
            double sG = magnitude( dM); 
            double s = sG / sM; 
            if (isnan(s))  s = 1;
            scaleSum += sM * s;
            sMsum += sM;  
        }

        scale = (double) scaleSum / sMsum;
        if (isnan(scale)) scale = 1; 
    }

    // Update output Scaling 
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        Coord2 loc = *gGuide[vertex].coordPtr; 
        _output(nRows) = ( ( (loc.x() - center.x() ) * scale) + center.x() ) * _wScale; 
        nRows++; 
        _output(nRows) = ( ( (loc.y() - center.y() ) * scale) + center.y() ) * _wScale; 
        nRows++; 
    }
}

/*
*   Input   vector containing coordinates
*   Output  center of the points
*
*   calculates the center of a list of points  
*/
Coord2 Stationmatching::getAllCenters( vector<Coord2> points) {
    if (points.size() == 1 ) return points[0]; 
    vector<Coord2> newPoints; 
    for (int i = 1; i < (int) points.size(); i++ ) {
        Coord2 p0 = points[i-1]; 
        Coord2 p1 = points[i]; 
        Coord2 d = p1 - p0;
        Coord2 c = p0 + ( d / 2.0);
        newPoints.push_back(c);  
    }
    return getAllCenters(newPoints);
}

/*
*   Input   none
*   Output  center of all matched guide nodes
*
*   calculates the center of all guide nodes
*/
Coord2 Stationmatching::getCenterOfGuideNodes() {
    vector<Coord2> gPoints; 
    UndirectedGraph &gGuide = _guide->g(); 
    for (int i = 0; i < (int) _stationDisc.size(); i++) {
        Coord2 p = *gGuide[_guideDisc[i]].coordPtr;
        gPoints.push_back(p); 
    }
    return getAllCenters(gPoints);
}

/*
*   Input   none
*   Output  none
*
*   initializing variables 
*/
void Stationmatching::_initVars( void ) 
{
    UndirectedGraph& gGuide = _guide->g();
    _var.resize( _nVars );
    _var.fill(0);
    
    unsigned int nRows = 0;
    BGL_FORALL_VERTICES( vertex, gGuide, UndirectedGraph) {
        Coord2 np = *gGuide[vertex].coordPtr; 
        _var(nRows * 2 + 0) = np.x(); 
        _var(nRows * 2 + 1) = np.y();  
        nRows++; 
    }
}

/*
*   Input   none
*   Output  none
*
*   initializing coefficients 
*/
void Stationmatching::_initCoefs( void ) 
{
    _coef.resize(_nConstrs, _nVars ); 
    _coef.fill( 0 );

    int nRows = 0;
    //Translate 
    for (int i = 0; i < (int) _nVertices; i++) {
        _coef(nRows, i * 2) = _wTransX;
        nRows++;  
        _coef(nRows, i * 2 + 1) = _wTransY; 
        nRows++; 
    }
    //Rotate
    for (int i = 0; i < (int) _nVertices; i++) {
        _coef(nRows, i * 2) = _wRotate; 
        nRows++; 
        _coef(nRows, i * 2 + 1) = _wRotate;
        nRows++; 
    }
    //Scaling 
    for (int i = 0; i < (int) _nVertices; i++) {
        _coef(nRows, i * 2) = _wScale; 
        nRows++; 
        _coef(nRows, i * 2 + 1) = _wScale; 
        nRows++; 
    }
}

/*
*   Input   none
*   Output  none
*
*   initializing outputs 
*/
void Stationmatching::_initOutputs( void )
{
    _output.resize(_nConstrs); 
    _output.fill(0); 
}

/*
*   Input   coordinates of node that should be removed
*   Output  true if has been removed
*
*   removes a node from matching  
*/
bool Stationmatching::removeNode(Coord2 cord) 
{
    UndirectedGraph & gMetro = _metro->g();
    UndirectedGraph & gGuide = _guide->g();
    
    int index = 0; 
    double minDist = 1e+15;  

    for (int i = 0; i < (int) _stationDisc.size(); i++) {
        Coord2 gp = *gGuide[_guideDisc[i]].coordPtr; 
        Coord2 sp = *gMetro[_stationDisc[i]].coordPtr; 
        double dg = magnitude(gp - cord); 
        double ds = magnitude(sp - cord); 
        double d = min(dg,ds); 
        if (d < minDist) {
            minDist = d; 
            index = i;
        }
    }

    if (minDist < distTol) {
        _stationDisc.erase(_stationDisc.begin() + index); 
        _guideDisc.erase(_guideDisc.begin() + index);
        cout << "removed Guide" << endl; 
        return true;
    }

    return false; 
}

/*
*   Input   coordinates of the added node
*   Output  true if new matching added
*
*   adds new node to matching
*/
bool Stationmatching::addNewNode(Coord2 cord) 
{
    UndirectedGraph        & gMetro            = _metro->g();
    UndirectedGraph        & gGuide             = _guide->g();

    Coord2 nearestCord = Coord2(0,0); 

    if (_guideDisc.size() == _stationDisc.size()) {
        double minDist = 100000; 
        VertexDescriptor nearest = VertexDescriptor(); 
        BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
            Coord2 vertCord = *gGuide[vertex].coordPtr;
            Coord2 delta = cord - vertCord; 
            double dist = magnitude(delta); 
            if (dist < minDist) {
                minDist = dist; 
                nearest = vertex; 
                nearestCord = vertCord; 
            }
        }
        if (minDist < distTol) {
            _guideDisc.push_back(nearest);
            _guideCoord.push_back(nearestCord);
            int id = gGuide[nearest].id; 
            // cout << "Added new Guide Node — id: " << id << endl;
            return true; 
        } else {
            // cerr << "Could not find a Guide Node close, Dist: " << minDist << endl;
            return false; 
        }
    }
    else if (_guideDisc.size() > _stationDisc.size()) {
        // cout << "–> Adding new Station Disc" << endl;
        double minDist = 10e+10; 
        VertexDescriptor nearest = VertexDescriptor(); 

        BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
            Coord2 vertCord = *gMetro[vertex].coordPtr;
            Coord2 delta = cord - vertCord; 
            double dist = magnitude(delta); 
            if (dist < minDist) {
                minDist = dist; 
                nearest = vertex;
                nearestCord = vertCord; 
            }
        }
        if (minDist < distTol) {
            _stationDisc.push_back(nearest);
            _stationCoord.push_back(nearestCord); 
            int id = gMetro[nearest].id;
            // cout << "–> Added new Station Node — id " << id << " with name " << *gMetro[nearest].namePtr << "; distance to guide" << gMetro[nearest].smoothDistance << ";\n closest porint on edge " << gMetro[nearest].closestPointOnEdge.x() << ", " << gMetro[nearest].closestPointOnEdge.y() << "\n position: " << *gMetro[nearest].coordPtr << endl;
            cout << "\'" << *gMetro[nearest].namePtr << "\'" << ", ";
            // cout << "closest point" << gMetro[nearest].closOnMetro << ", distance " << (*gMetro[nearest].coordPtr - gMetro[nearest].closOnMetro).norm() << endl;
            return true;
        } else {
            // cerr << "Could not find a Metro Station close, Dist: " << minDist << endl;
            return false; 
        }
    }
    return false; 
}

/*
*   Input   Coordinate of vector
*   Output  magnitude
*
*   calculates magnitude of vector from 0,0 to v
*/
double Stationmatching::magnitude(Coord2 v) 
{
    double r = v.x() * v.x() + v.y() * v.y(); 
    r = sqrt(r); 
    return r; 
}

/* 
*   Input   center of cord. system and cord of point
*   Output  point in Polar cord. [r, angle]
* 
*   Coverts points from cat. cord. system to polar cord. system
*/ 
Coord2 Stationmatching::getAsPolar(Coord2 center, Coord2 point) {
    Coord2 p = point - center; 
    double r = sqrt(p.x() * p.x() + p.y() * p.y()); 
    double a = atan2(p.y(), p.x()); 
    return Coord2(r, a);
}

/*
*   Input   center of pol cord. system and point [r, angle]
*   Output  Coordinates in car. system
*   
*   Converts point from Pol.coord. system to Cart. coord. system
*/  
Coord2 Stationmatching::polarAsCart(Coord2 center, Coord2 point) {
    double x = point.x() * cos(point.y()); 
    double y = point.x() * sin(point.y()); 
    Coord2 p = Coord2(x,y); 
    return (p + center); 
}

Stationmatching::Stationmatching( void ) {}
Stationmatching::~Stationmatching( void ) {}