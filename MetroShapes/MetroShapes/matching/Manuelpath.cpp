//
// Created by tobias batik on 26.09.21.
//

#include "Manuelpath.h"

void Manuelpath::init(Metro * __metro, Guide * __guide) {
    _metro = __metro;
    _guide = __guide;
    _maxDistanceVertex = 5;
}

bool Manuelpath::addVertex( int id ) {
    UndirectedGraph        & gMetro        = _metro->g();
    bool successful = false;
    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        if ( gMetro[vertex].id == id ) {
            _path.push_back(vertex);
            _pathId.push_back(id);
            successful = true;
        }
    }

    for (int i = 0; i < _path.size(); i++) {
        VertexDescriptor vdI = _path[i];
        for (int j = 0; j < _path.size(); j++ ) {
            if ( i != j ) {
                VertexDescriptor vdJ = _path[j];
                if (edge(vdI, vdJ, gMetro).second) {
                    EdgeDescriptor ed = edge(vdI, vdJ, gMetro).first;
                    gMetro[ed].matchPath = true;
                    gMetro[vdI].autoPath = true;
                }
            }
        }
    }
    return successful;
}

bool Manuelpath::addVertex( Coord2 coord ) {
    UndirectedGraph        & gMetro        = _metro->g();
    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        auto stationCoord = *gMetro[vertex].coordPtr;
        auto delta = coord - stationCoord;
        if ( delta.norm() < _maxDistanceVertex ) {
            return addVertex(gMetro[vertex].id);
        }
    }
}

bool Manuelpath::addVertices( vector<int> ids ) {
    bool successful = true;
    for ( int id : ids ) {
        if (!addVertex(id))
            successful = false;
    }
    return successful;
}

void Manuelpath::alingGuide() {
    scaleGuide();

    auto guideCenter = _guide->getCenterCoord();
    auto pathCenter = getCenterPath();
    _guide->translate(pathCenter - guideCenter);

    // guideCenter = _guide->getCenterCoord();
    // pathCenter = getCenterPath();
    // _guide->translate(pathCenter - guideCenter);

}

Coord2 Manuelpath::getCenterPath() {
    UndirectedGraph        & gMetro        = _metro->g();

    Coord2 result = Coord2(0., 0.);
    double minX = 10e5;
    double minY = 10e5;
    double maxX = -10e5;
    double maxY = -10e5;

    for (auto vertex : _path) {
        Coord2 coord = *gMetro[vertex].coordPtr;
        minX = min(minX, coord.x());
        minY = min(minY, coord.y());
        maxX = max(maxX, coord.x());
        maxY = max(maxY, coord.y());
    }

    double x = (maxX - minX) * 0.5 + minX;
    double y = (maxY - minY) * 0.5 + minY;
    result.setX(x);
    result.setY(y);

    return result;
}

void Manuelpath::scaleGuide() {
    UndirectedGraph        & gMetro        = _metro->g();

    if (_path.size() > 1 ) {
        double minX = 10e5;
        double minY = 10e5;
        double maxX = -10e5;
        double maxY = -10e5;

        for (auto vertex : _path) {
            Coord2 coord = *gMetro[vertex].coordPtr;
            minX = min(minX, coord.x());
            minY = min(minY, coord.y());
            maxX = max(maxX, coord.x());
            maxY = max(maxY, coord.y());
        }
        double width = maxX - minX;
        double height = maxY - minY;

        double guideWidth = _guide->getWidth();
        double guideHeight = _guide->getHeight();

        double deltaWidth =  width / guideWidth;
        double deltaHeight = height / guideHeight;

        Coord2 guideCenter = _guide->getCenterCoord();
        _guide->translate(Coord2(-guideCenter.x(), -guideCenter.y()));
        _guide->scale(max(deltaWidth, deltaHeight));
        _guide->translate(guideCenter);

        // if (_deformUnuniformly && _path.size() > 5) {
        //    auto pathCenter = getCenterPath();
        //    _metro->translate(Coord2(pathCenter.x() * -1., pathCenter.y() * -1.));
        //    if (deltaWidth < deltaHeight)
        //        _metro->scale(Coord2(1., 1.0 + deltaWidth));
        //    else
        //        _metro->scale(Coord2(1.0 + deltaHeight, 1.0));
        //    _metro->translate(pathCenter);
        // }
    }
}

void Manuelpath::scaleMetroNonUniformly() {
    cout << "scaling Metro Non-Uniformly" << endl;
    UndirectedGraph        & gMetro        = _metro->g();

    if (_path.size() > 1 ) {
        double minX = 10e5;
        double minY = 10e5;
        double maxX = -10e5;
        double maxY = -10e5;

        for (auto vertex : _path) {
            Coord2 coord = *gMetro[vertex].coordPtr;
            minX = min(minX, coord.x());
            minY = min(minY, coord.y());
            maxX = max(maxX, coord.x());
            maxY = max(maxY, coord.y());
        }
        double width = maxX - minX;
        double height = maxY - minY;

        double guideWidth = _guide->getWidth();
        double guideHeight = _guide->getHeight();

        double deltaWidth =   guideWidth / width;
        double deltaHeight =  guideHeight / height;


        if (_deformUnuniformly and _path.size() > 1) {
           auto pathCenter = getCenterPath();
           _metro->translate(Coord2(pathCenter.x() * -1., pathCenter.y() * -1.));

           if (deltaWidth < deltaHeight) // The smaller one is done by the path scaling?
               _metro->scale(Coord2(1.0, deltaHeight));
           else
               _metro->scale(Coord2(deltaWidth, 1.0));
           _metro->translate(pathCenter);
        }
    }

}