//
// Created by tobias batik on 28.04.21.
//

#include "AutoMatching.h"

void AutoMatching::init(Metro *__metro, Guide *__guide, double __half_width, double __half_height) {
    _metro = __metro;
    _guide = __guide;
    _tolerance = M_PI_4 * .5;

    _nVertices = 0;

    UndirectedGraph        & gGuide        = _guide->g();
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        _guideVds.push_back(vertex);
    }
}

void AutoMatching::findPath() {
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();

    double bestScore = 0.0;
    double minError = 10e4;

    vector<VertexDescriptor> bestPath;

    int index = 0;
    cerr << "Starting the space djigstars" << endl;
    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        // if (gMetro[vertex].id == 56) {
        if (true) {
            if (out_degree( vertex, gMetro ) > 3 || true ) {
                auto spaceDjigstar = spaceDijgstar(vertex);
                auto path = spaceDjigstar.first;
                auto error = spaceDjigstar.second;

                cout << "\nfound path ";
                if (path.size() > 1 ) {
                    for (auto i : path) {
                        cout << ", " << gMetro[i].id;
                    }
                    cout << endl;

                    cout << "   error: " << error << endl;
                    if (error < minError ){
                        // if (bestPath.size() < path.size()) {
                        minError = error;
                        bestPath = path;
                    } else {
                        cout << "    not good enough error: " << error << "/" << minError << endl;
                    }
                }
            }
        }
        index++;
    }

    cout << "\nFound a Path with error: " << minError << endl;
    index = 0;
    cout << "   with ids: ";

    BGL_FORALL_EDGES(edge, gMetro, UndirectedGraph) {
        gMetro[edge].matchPath = false;
    }

    for (auto v_desc : bestPath ) {
        gMetro[v_desc].autoPath = true;
        if (index > 0 && edge(v_desc, bestPath[index - 1], gMetro).second) {
            cout << gMetro[v_desc].id << ", " << gMetro[bestPath[index - 1]].id << ", ";
            EdgeDescriptor ed = edge(v_desc, bestPath[index - 1], gMetro).first;
            gMetro[ed].matchPath = true;
        }
        index++;
    }
    cout << endl;

    _foundMetroPath = bestPath;

#ifdef TESTSIMILARITYFORALLPATHS
    // TEST Path similiarity messerument
    cout << "------------------- TEST Path similiarity messerument ------------------" << endl;
    // guideVD
    vector<VertexDescriptor> guideVD;
    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        guideVD.push_back(vertex);
    }
    vector<VertexDescriptor> metroVD;
    int indexCount = 0;
    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        if (indexCount == 87 || indexCount == 198 || indexCount == 175 ) metroVD.push_back(vertex);
        indexCount++;
    }
    double cost = computeDFD(metroVD, guideVD);
    double partCost = computePartDFD(metroVD, guideVD, 0).first;
    double partLenght = computePartDFD(metroVD, guideVD, 0).second;
    cout << "       similarity between the two paths: " << cost << endl;
    cout << "       part Similarity:                  " << partCost << ", " << partLenght << endl;
    cout << "--------------------------------------------------------------------------\n" << endl;
#endif
}

void AutoMatching::translateGuide() {
    cout << "Translating Guide ..." << endl;
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();

    Coord2 metroMin = getMinValues(_foundMetroPath, gMetro);
    cout << "metro Min value " << metroMin << endl;
    Coord2 guideMin = getMinValues(_guideVds, gGuide);
    Coord2 deltaMin = metroMin - guideMin;
    // _guide->translate(deltaMin);

    // scaling
    Coord2 center = _guide->getCenterCoord();
    Coord2 guideMax = getMaxValues(_guideVds, gGuide);
    guideMin = getMinValues(_guideVds, gGuide);
    Coord2 guideDelta = guideMax - guideMin;
    Coord2 metroMax = getMaxValues(_foundMetroPath, gMetro);
    Coord2 metroDelta = metroMax - metroMin;

    double sy = metroDelta.y() / guideDelta.y();
    double sx = metroDelta.x() / guideDelta.x();
    Coord2 centerGuide = _guide->getCenterCoord();
    Coord2 centerPath = getCenterOfPath(_foundMetroPath, gMetro);

    // scale to same size
    _guide->translate(Coord2(centerGuide.x() * -1., centerGuide.y() * -1.));
    if (sx >= sy) { // x difference is smaller
        _guide->scale(sx);
        _guide->translate(centerGuide);

        if (_deformUnuniformly) {
            _metro->translate(Coord2(centerPath.x() * -1., centerPath.y() * -1.));
            _metro->scale(Coord2(1., sx / sy));
            _metro->translate(centerPath);
        }
    } else {
        _guide->scale(sy);
        _guide->translate(centerGuide);

        if (_deformUnuniformly) {
            _metro->translate(Coord2(centerPath.x() * -1., centerPath.y() * -1.));
            _metro->scale(Coord2(sx * sy, 1.0));
            _metro->translate(centerPath);
        }
    }


    // align center coordinates
    centerPath = getCenterOfPath(_foundMetroPath, gMetro);
    centerGuide = _guide->getCenterCoord();
    Coord2 deltaCenter = Coord2(centerPath.x() - centerGuide.x(), centerPath.y() - centerGuide.y());
    _guide->translate(deltaCenter);

    centerPath = getCenterOfPath(_foundMetroPath, gMetro);
    centerGuide = _guide->getCenterCoord();
    deltaCenter = Coord2(centerPath.x() - centerGuide.x(), centerPath.y() - centerGuide.y());
    _guide->translate(deltaCenter);

}


pair<vector<VertexDescriptor>, double> AutoMatching::spaceDijgstar(VertexDescriptor vi) {
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();

    vector<VertexDescriptor> guideDis;

    BGL_FORALL_VERTICES(vertex, gGuide, UndirectedGraph) {
        guideDis.push_back(vertex);
    }

    map<VertexDescriptor , double> distances;
    map<VertexDescriptor , double> partDistances;
    map<VertexDescriptor , double> pathLength;
    map<VertexDescriptor , int> partMatchIndex;
    set<VertexDescriptor> vertexQue;
    set<VertexDescriptor> removedVertex;
    map<VertexDescriptor, VertexDescriptor> previousSelected;
    set<VertexDescriptor> notSetVertex;

    BGL_FORALL_VERTICES(vertex, gMetro, UndirectedGraph) {
        distances.insert(pair<VertexDescriptor, double>(vertex, MAX_COST));
        partDistances.insert(pair<VertexDescriptor, double>(vertex, MAX_COST));
        pathLength.insert(pair<VertexDescriptor, double>(vertex, 0.0));
        partMatchIndex.insert(pair<VertexDescriptor, int>(vertex, 0));
        notSetVertex.insert(vertex);
        previousSelected.insert(pair<VertexDescriptor, VertexDescriptor>(vertex, nullptr));
    }

    notSetVertex.erase(vi);
    vertexQue.insert(vi);

    int counter = 0;
    // main Alogirthm
    while(vertexQue.size() > 0) {
        pair<VertexDescriptor, double> selectedVertex = getShortestDistanceInQuee(vertexQue, partDistances);
        vertexQue.erase(selectedVertex.first); // for the first Vertex
        notSetVertex.erase(selectedVertex.first); // can not gurante that there is a "longer" way with additional edges
        auto metroPath_u = getPathVertex(selectedVertex.first, previousSelected);
        auto neighbors = getAllNeighboarsNotRemoved(selectedVertex.first, vertexQue ,notSetVertex);


        // Interate over all neighbors of u
        for (auto neighbor : neighbors) {
            auto metroPath_v = metroPath_u;
            metroPath_v.push_back(neighbor.first);
            // VertexDescriptor firstVert = metroPath_v[0]
            double costColors = 1;
#ifdef LESSSECTIONS
            double numberColors = pathNumberOfColors(metroPath_v);
            costColors = 1.0 + numberColors * _weightColor;
#endif

            auto partMatch = computePartDFD(metroPath_v, guideDis, partMatchIndex[selectedVertex.first]);
            double fd_v_part = partMatch.first * costColors;
            double fd_v_full = computeDFD(metroPath_v, guideDis) * costColors; // M[metroPath_v.size() - 1][guideDis.size() - 1];
            int matchIndex = partMatch.second;
            if (isnan(fd_v_part)) cerr << "NaN: spaceDjigstar: fd_v_part is " << fd_v_part << endl;


            if (fd_v_full < distances[selectedVertex.first]) {
                distances[neighbor.first] = fd_v_full ;
                partDistances[neighbor.first] = fd_v_part;
                pathLength[neighbor.first] = metroPath_v.size();
                previousSelected[neighbor.first] = selectedVertex.first;
                partMatchIndex[neighbor.first] = matchIndex;
                vertexQue.insert(neighbor.first);
            } else {
#ifdef DEBUG
                cerr << "fd_v_part: " << fd_v_part << endl;
#endif
            }

        }
        counter += 1;
    }

    // Find the best Path
    VertexDescriptor bestVert;
    double minCost = MAX_COST;
    double maxLenght = -1;
    double maxMatchIndex = 0;
    int in = 0;

    for (auto dist : distances) {
        if (abs(minCost) > abs(dist.second)) {
            maxLenght = pathLength[dist.first];
            bestVert = dist.first;
            minCost = dist.second;
            maxMatchIndex = partMatchIndex[dist.first];
        }
        in++;
    }

    auto path = getPathVertex(bestVert, previousSelected); //TODO

    path = getPathVertex(bestVert, previousSelected);
    // minCost = computeDFD(path, guideDis);
    int index = 0;
    for (auto v_desc : path ) {
        if ( index > 0 && edge(v_desc, path[index-1], gMetro).second) {
            EdgeDescriptor  ed = edge(v_desc, path[index-1], gMetro).first;
        } else {
            // cerr << "found an edge that does not exist; index " << index << endl;
        }
        index++;
    }
    return pair<vector<VertexDescriptor>, double>(path, minCost);
}

// TODO
int AutoMatching::pathNumberOfColors(vector<VertexDescriptor> path) {
    if (path.size() < 1 ) {
        return 10e10;
    } else if (path.size() == 1 ) {
        return 1;
    }
    int numberSegments = 1;
    UndirectedGraph &gMetro = _metro->g();
    VertexDescriptor vi = path[0];
    VertexDescriptor vj = path[1];
    auto currentEdge = edge(vi, vj, gMetro);
    if (currentEdge.second) {
        EdgeDescriptor edgeDis = currentEdge.first;
        vector< unsigned int > lastLineID = gMetro[currentEdge.first].lineID;
        for (unsigned int i = 1; i < path.size(); i++ ) {
            vi = vj;
            vj = path[i];
            currentEdge = edge(vi, vj, gMetro);
            auto currentEdge2 = edge(vj, vi, gMetro);
            if (currentEdge.second != currentEdge2.second) cerr << "They differer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            if (currentEdge.second) {
                auto currentLineID = gMetro[currentEdge.first].lineID;
                // if two vectors do not share a lineID numbersegments++
                bool share = false;
                for (auto idCur : currentLineID) {
                    for (auto idLast : lastLineID ) {
                        if (idCur == idLast ) {
                            share = true;
                        }
                    }
                }
                if (gMetro[currentEdge.first].isMetro) lastLineID = currentLineID;
                if (!share) numberSegments++;
            } else {
                // cerr << "This edge does not exist " << endl;
            }
        }
    } else {
        cerr << "can not calc number of colors because invalide edge" << endl;
    }
    // cout << "number segments = " << numberSegments << endl;
    return numberSegments / path.size();
}

Coord2 AutoMatching::getCenterOfPath(vector< VertexDescriptor> path, UndirectedGraph g) {
    Coord2 avgCord = Coord2(0,0);
    double count = 0;
    for (int i = 0; i < path.size(); i++) {
        Coord2 cordTem = *g[path[i]].coordPtr;
        avgCord += cordTem;
        count++;
    }
    // return avgCord / count;

    //Bounding Box Center
    double minX = 10e10;
    double minY = 10e10;
    double maxX = -10e10;
    double maxY = -10e10;

    for (int i = 0; i < path.size(); i++) {
        Coord2 c = *g[path[i]].coordPtr;
        if (c.x() < minX) {
            minX = c.x();
        } else if (c.x() > maxX) {
            maxX = c.x();
        }
        if (c.y() < minY) {
            minY = c.y();
        }else if (c.y() > maxY) {
            maxY = c.y();
        }
    }

    double width = maxX - minX;
    double height = maxY - minY;
    Coord2 center = Coord2(0., 0.);
    center.setX(width * 0.5 + minX);
    center.setY(height * 0.5 + minY);
    return center;
}

pair<double, int> AutoMatching::computePartDFD(vector<VertexDescriptor> metroPathVD, vector<VertexDescriptor> guideVD, int minMatch){
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();

    metroPathCd.erase(metroPathCd.begin(), metroPathCd.end());
    _metroPathVD.erase(_metroPathVD.begin(), _metroPathVD.end());
    guideCD.erase(guideCD.begin(), guideCD.end());
    metroVertPar.erase(metroVertPar.begin(), metroVertPar.end());
    guideVertPar.erase( guideVertPar.begin(), guideVertPar.end());

    for (auto vd : metroPathVD ) {
        Coord2 v = *gMetro[vd].coordPtr;
        metroPathCd.push_back(v);
        _metroPathVD.push_back(vd);
    }

    for (auto vd : guideVD ) {
        Coord2 g = *gGuide[vd].coordPtr;
        guideCD.push_back(g);
    }
    double metroLength = 0.0;
    for ( int i = 1; i < metroPathCd.size(); i++ ) {
        Coord2 delta = metroPathCd[i] - metroPathCd[i-1];
        metroLength += delta.norm();
    }
    double guideLength = 0.0;
    for (int i = 1; i < guideCD.size(); i++ ) {
        Coord2 delta = guideCD[i] - guideCD[i-1];
        guideLength += delta.norm();
    }
    metroVertPar.push_back(0);
    for (int i = 1; i < metroPathCd.size(); i++ ) {
        Coord2 delta = metroPathCd[i] - metroPathCd[i-1];
        double value = metroVertPar[i-1] + (delta.norm() / metroLength );
        metroVertPar.push_back(value);
    }
    guideVertPar.push_back(0);
    for (int i = 1; i < guideCD.size(); i++ ) {
        Coord2 delta = guideCD[i] - guideCD[i-1];
        double value = guideVertPar[i-1] + (delta.norm() / guideLength );

        guideVertPar.push_back(value);
    }
    // init M
    M.erase(M.begin(), M.end());
    for (int i = 0; i < metroPathCd.size(); i++ ) {
        vector<double> row;
        for (int j = 0; j < guideCD.size(); j++ ) {
            row.push_back(-10e+4);
        }
        M.push_back(row);
    }

    computeSegmentDFD(metroPathCd.size() - 1, guideCD.size() -1);

    int i = metroPathCd.size() - 1;

    if (minMatch == 0 ) minMatch = 1;

    double result = M[i][minMatch];
    int matchIndex = minMatch;

#ifdef  NORMALIZEPARTMATCHING
    result = M[i][minMatch] * (1 + 1) / (1.001 + guideVertPar[minMatch]);
    for (int j = minMatch; j < guideCD.size(); j++ ) {
        if (result >= abs(M[i][j]) * (1 + 1) / (1.001 + guideVertPar[j]) ) {
            result = M[i][j] * (1 + 1) / (1.001 + guideVertPar[j]);
            matchIndex = j;
        }
    }
#endif
#ifndef NORMALIZEPARTMATCHING
    for (int j = minMatch; j < guideCD.size(); j++) {
        if (result >= abs(M[i][j])) {
            result = M[i][j];
            matchIndex = j;
        }
    }
#endif

    if (isnan(result)) cerr << "computDFD result is: " << result << endl;

    #ifdef DEBUG
    cout << "------------- Print M ---------------------- " << endl;
    for (int i = 0; i < M.size(); i ++ ) {
        for (int j = 0; j < M[i].size(); j++ ) {
            cout << M[i][j] << ", ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------- " << endl;
    cout << "result part " << result << endl;
    #endif


    return pair<double, int>(result, matchIndex );
}


double AutoMatching::computeDFD(vector<VertexDescriptor> metroPathVD, vector<VertexDescriptor> guideVD) {
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();
    //vector<Coord2> metroPathCd, guideCD;
    metroPathCd.erase(metroPathCd.begin(), metroPathCd.end());
    guideCD.erase(guideCD.begin(), guideCD.end());
    metroVertPar.erase(metroVertPar.begin(), metroVertPar.end());
    guideVertPar.erase( guideVertPar.begin(), guideVertPar.end());

    for (auto vd : metroPathVD ) {
        Coord2 v = *gMetro[vd].coordPtr;
        metroPathCd.push_back(v);
    }
    for (auto vd : guideVD ) {
        Coord2 g = *gGuide[vd].coordPtr;
        guideCD.push_back(g);
    }
    double metroLength = 0.0;
    for ( int i = 1; i < metroPathCd.size(); i++ ) {
        Coord2 delta = metroPathCd[i] - metroPathCd[i-1];
        metroLength += delta.norm();
    }
    double guideLength = 0.0;
    for (int i = 1; i < guideCD.size(); i++ ) {
        Coord2 delta = guideCD[i] - guideCD[i-1];
        guideLength += delta.norm();
    }
    // vector<double> metroVertPar, guideVertPar;
    metroVertPar.push_back(0);
    for (int i = 1; i < metroPathCd.size(); i++ ) {
        Coord2 delta = metroPathCd[i] - metroPathCd[i-1];
        double value = metroVertPar[i-1] + (delta.norm() / metroLength );
        metroVertPar.push_back(value);
    }
    guideVertPar.push_back(0);
    for (int i = 1; i < guideCD.size(); i++ ) {
        Coord2 delta = guideCD[i] - guideCD[i-1];
        double value = guideVertPar[i-1] + (delta.norm() / guideLength );

        guideVertPar.push_back(value);
    }
    // init M
    M.erase(M.begin(), M.end());
    for (int i = 0; i < metroPathCd.size(); i++ ) {
        vector<double> row;
        for (int j = 0; j < guideCD.size(); j++ ) {
            row.push_back(-10e+04);
        }
        M.push_back(row);
    }

    double result = computeSegmentDFD(metroPathCd.size() - 1, guideCD.size() -1);
    return result;
}

double AutoMatching::computeSegmentDFD( int i, int j) {
    if ( M[i][j] > -10e+3) {
        return M[i][j];
    }
    else if ( i == 0 && j == 0) {
        double c_uivj = costSegmentDFD(i, j);
        M[i][j] = c_uivj;
    } else if ( i > 0 && j == 0) {
        double newVal = computeSegmentDFD(i-1, 0);
        double c_uiv1 = costSegmentDFD(i, 0);
        M[i][j] = newVal + c_uiv1;
    } else if ( i == 0 && j > 0 ) {
        double newVal = computeSegmentDFD(0, j -1);
        double c_u1vj = costSegmentDFD(0, j);
        M[i][j] = newVal + c_u1vj;
    } else if (i > 0 && j > 0 ) {
        double c1 = computeSegmentDFD(i - 1,  j);
        double c2 = computeSegmentDFD(i - 1, j - 1);
        double c3 = computeSegmentDFD(i, j - 1);
        double mini = min(c1, c2);  // TODO min stat max
        mini = min(mini, c3);   // TODO min stat max
        double c_uivj = costSegmentDFD(i, j);
        M[i][j] = mini + c_uivj; // TODO
        // M[i][j] = max(mini, c_uivj);
        // M[i][j] = c2 + c_uivj;
    } else {
        M[i][j] = -10e+4;
        cerr << "else No soultion? i: " << i << " metroPathLength: " << metroPathCd.size() << ", j: " << j << endl;

    }
    // if (isnan(M[i][j])) cerr << "NaN: compute segmentDFD: M_ij is " << M[i][j] << endl;
    return M[i][j];
}

double AutoMatching::costSegmentDFD(int i, int j) {
    // if (i == 0 || j == 0) return 0.01; // TODO no idea why 0.001 fixes problems
    if (i == 0 ) i++;
    if (j == 0 ) j++;

    auto metroV = metroPathCd[i];
    auto metroU = metroPathCd[i-1];
    Coord2 metroDelta = metroV - metroU;

    auto guideV = guideCD[j];
    auto guideU = guideCD[j-1];
    Coord2 guideDelta = guideV - guideU;

    double dist = metroVertPar[i] - metroVertPar[i-1];
    dist += guideVertPar[j] - guideVertPar[j-1];
    double  result;
    double norms = (metroDelta.norm() * guideDelta.norm());
    if (norms == 0) norms = 0.001;
    double ax = (metroDelta.x() * guideDelta.x() + metroDelta.y() * guideDelta.y()) / norms;

    double angle = acos( (min(max(ax,-1.0),1.0))); // Returns always a positive angle; calculate signum?
    result = angle * dist;

    UndirectedGraph        & gMetro        = _metro->g();
    VertexDescriptor vi = _metroPathVD[i];
    VertexDescriptor vj = _metroPathVD[i-1];

    auto fedge = edge(_metroPathVD[i], _metroPathVD[i-1], gMetro);
    if (fedge.second) {
        EdgeDescriptor  ed = fedge.first;
        if (gMetro[ed].isMetro) {
            result *= _costMetro;
        } else {
            result *= _costAdd;
        }
    } else {
        cerr << "could not find the edge when calc weight" << endl;
    }


        #ifdef DEBUG
    if (isnan(ax)) {
        cerr << "NaN: cost: ax is " << ax << endl;
        cerr << "   metroDelta " << metroDelta.x() << ", " << metroDelta.y() << endl;
        cerr << "   guideDelta " << guideDelta.x() << ", " << guideDelta.y() << endl;
    }
    if (isnan(angle)) cerr << "NaN: cost: angie is " << angle << endl;
    if (isnan(result)) cerr << "NaN: cost: result is " << result << endl;
        #endif

    return  result;
}

map<VertexDescriptor, EdgeDescriptor> AutoMatching::getAllNeighboarsNotRemoved(VertexDescriptor selectedVertex,
                                                                               set<VertexDescriptor> vertexQuee,
                                                                               set<VertexDescriptor> notSetVertex) {
    UndirectedGraph &gMetro = _metro->g();
    map<VertexDescriptor, EdgeDescriptor> result;

    OutEdgeIterator e, e_end;
    for (tie(e, e_end) = out_edges(selectedVertex, gMetro); e != e_end; ++e) {
        EdgeDescriptor ed = *e;
        VertexDescriptor vj = source(ed, gMetro);
        if (vj == selectedVertex) vj = target(ed, gMetro);
        if (notSetVertex.find(vj) != notSetVertex.end()) {
            result.insert(pair<VertexDescriptor, EdgeDescriptor>(vj, ed));
            // cout << gMetro[vj].id << endl;
        }
    }
    return result;
}

map<VertexDescriptor, EdgeDescriptor> AutoMatching::getAllNeighboarsNotRemovedOrStart(VertexDescriptor selectedVertex,
                                                                                      set<VertexDescriptor> notSetVertex,
                                                                                      VertexDescriptor vi) {
    UndirectedGraph &gMetro = _metro->g();
    map<VertexDescriptor, EdgeDescriptor> result;

    OutEdgeIterator e, e_end;
    for (tie(e, e_end) = out_edges(selectedVertex, gMetro); e != e_end; ++e) {
        EdgeDescriptor ed = *e;
        VertexDescriptor vj = source(ed, gMetro);
        if (vj == selectedVertex) vj = target(ed, gMetro);
        if ( notSetVertex.find(vj) != notSetVertex.end() || vj == vi ) {
            result.insert(pair<VertexDescriptor, EdgeDescriptor>(vj, ed));
            if (vj == vi) cout << "returned start: " << gMetro[vj].id << endl;
        }
    }
    return result;
}

vector<double> AutoMatching::convertToIPR(vector<VertexDescriptor> vds, UndirectedGraph graph) {
    return convertToIPR(vds, graph, false);
}

vector<double> AutoMatching::convertToIPR(vector<VertexDescriptor> vds, UndirectedGraph graph, bool final) {
    vector<double> result;
    auto coords = convertVertDisToCoord(vds, graph);

    if (coords.size() > 0 ) {
        auto u0 = coords[0];
        auto v0 = coords[1];
        auto uv0 = v0 - u0;
        double a_uv0 = atan2(uv0.y(), uv0.x());
        result.push_back(a_uv0);
        if (final) graph[vds[0]].inflectionPoint = true;
    }

    for (int i = 2; i < sizeof(coords) - 1; i++ ) {
        auto t = coords[i-2];
        auto u = coords[i-1];
        auto v = coords[i];
        auto w = coords[i+1];
        Coord2 tu = u - t;
        Coord2 uv = v - u;
        Coord2 vw = w - v;
        double a_tu = atan2(tu.y(), uv.x());            // TODO not every loop new
        double a_uv = atan2(uv.y(), uv.x());
        double a_vw = atan2(vw.y(), vw.x());

        double a_dif1 = a_uv - a_tu;
        double a_dif2 = a_vw - a_uv;
        bool sig1, sig2;
        a_dif1 > 0 ? sig1 = true : sig1 = false;
        a_dif2 > 0 ? sig2 = true : sig2 = false;

        // https://stackoverflow.com/questions/60030983/finding-inflection-points-of-a-set-of-points
        double winding2 = uv.x() * vw.y() - uv.y() * vw.x();                 //a1b2 - a2b1
        double winding1 = tu.x() * uv.y() - tu.y() * uv.x();

        if ( winding2 * winding1 < 0.0 ) {
            // if (sig1 != sig2) {                         // is ref point
            double at = ( a_vw + a_uv ) * .5;
            result.push_back(at);
            // VertexDescriptor vdi = vds[i];
            if(final) graph[vds[i]].inflectionPoint = true;
        }
    }
    if (coords.size() > 1 ) {
        int last = coords.size() - 1;
        auto u0 = coords[last-1];
        auto v0 = coords[last];
        auto uv0 = v0 - u0;
        double a_uv0 = atan2(uv0.y(), uv0.x());
        result.push_back(a_uv0);
        if (final) graph[vds[vds.size() - 1 ]].inflectionPoint = true;       //todo
    }

    return result;
}

int AutoMatching::longestPrefix(vector<double> mapPathIPR, vector<double> shapeIPR, double tolerance) {
    for (int i = 0; i < shapeIPR.size() && i < mapPathIPR.size(); i++) {
        if (abs(shapeIPR[i] - mapPathIPR[i]) > tolerance ) { // &&
            // abs(shapeIPR[i] - mapPathIPR[i] - M_PI) > tolerance ) {
            return i;
        } else {
            // cout << "was true: shape " << shapeIPR[i] << " map " << mapPathIPR[i] << endl;
            // cout << "    " << abs(shapeIPR[i] - mapPathIPR[i]) << endl;
            // cout << "    " << abs(shapeIPR[i] - mapPathIPR[i] - M_PI) << endl;
        }
    }
    return ( mapPathIPR.size() < shapeIPR.size() ) ? mapPathIPR.size() : shapeIPR.size();
}

vector<VertexDescriptor> AutoMatching::getPathVertex(VertexDescriptor endVert,
                                                     map<VertexDescriptor, VertexDescriptor> previousSelected) {
    vector<VertexDescriptor> result;
    VertexDescriptor vi = endVert;
    UndirectedGraph        & gMetro        = _metro->g();

    while(vi != nullptr) {
        result.insert(result.begin(), vi);
        vi = previousSelected[vi];
    }

    return result;
}

pair<double, double> AutoMatching::pathMatchingIPRScore(vector<VertexDescriptor> mapPathVD, vector<VertexDescriptor> shapeVD, double tolerance) {
    UndirectedGraph        & gMetro        = _metro->g();
    UndirectedGraph        & gGuide        = _guide->g();
    auto mapIPR = convertToIPR(mapPathVD, gMetro);
    auto shapeIPR = convertToIPR(shapeVD, gGuide);
    double score = 0.0;
    double errorSum = 0.0;

    for (int i = 0; i < shapeIPR.size() && i < mapIPR.size(); i++) {
        double error = abs(shapeIPR[i] - mapIPR[i]);
        if ( error < tolerance ) {
            score += 1.0;// / ((double) shapeIPR.size());
            errorSum += error;
        }
    }
    return pair<double, double>(score, errorSum);
}

bool AutoMatching::pathSharePrefix(vector<double> mapPathIPR, vector<double> shapeIPR, double tolerance ) {
    if (mapPathIPR.size() == 0 ) return true;

    int i_shape = 0;
    int i_map = 0;

    // while(i_shape < shapeIPR.size() && i_map < mapPathIPR.size() ) {

    // }


    if (shapeIPR.size() >= mapPathIPR.size() ) {
        vector<double> diffs;
        for (int i = 0; i < shapeIPR.size(); i++ ) {
            double diff = abs(mapPathIPR[i] - shapeIPR[i]);
            if( diff < tolerance ) {
                diffs.push_back(diff);
            } else {
                break;
            }
        }

        if (diffs.size() > 0) {
            return true;
        } else {
            return false;
        }
    } else {    // shape path can not be longer?
        cout << "shape path only ";
        for (auto i : shapeIPR) {
            cout << i << ", ";
        }
        cout << endl;
        return false;
    }
    return true;
}

vector<Coord2> AutoMatching::convertVertDisToCoord(vector<VertexDescriptor> vds, UndirectedGraph graph) {
    vector<Coord2> result;
    for (VertexDescriptor vd : vds) {
        Coord2* vc = graph[vd].geoPtr;
        result.push_back(*vc);
    }

    return result;
}

Coord2 AutoMatching::getMinValues(vector<VertexDescriptor> vds, UndirectedGraph graph) {
    double minX = 10e10;
    double minY = 10e10;
    for (auto vertex : vds ) {
        Coord2 v = *graph[vertex].coordPtr;
        if (v.x() < minX) {
            minX = v.x();
        }
        if (v.y() < minY) {
            minY = v.y();
        }
    }

    return Coord2(minX, minY);
}

Coord2 AutoMatching::getMaxValues(vector<VertexDescriptor> vds, UndirectedGraph graph) {
    double maxX = -10e10;
    double maxY = -10e10;
    for (auto vertex : vds ) {
        Coord2 v = *graph[vertex].coordPtr;
        if (v.x() > maxX) {
            maxX = v.x();
        }
        if (v.y() > maxY) {
            maxY = v.y();
        }
    }

    return Coord2(maxX, maxY);
}
/*
 * returns vert dis, min Dist
 */
pair<VertexDescriptor, double> AutoMatching::getShortestDistanceInQuee(set<VertexDescriptor> vertexQuee, map<VertexDescriptor, double> distances) {

    VertexDescriptor minVertex;
    double minDist = 10e10;

    for (auto element : vertexQuee ) {
        // double dist = distances[element].first;
        double dist = distances[element];
        if (abs(dist) < abs(minDist)) {
            minDist = dist;
            minVertex = element;
        }
    }
    return pair<VertexDescriptor, double>(minVertex, minDist);

}

map<VertexDescriptor, EdgeDescriptor> AutoMatching::getAllNeighborsInQuee(VertexDescriptor selectedVertex,
                                                                          set<VertexDescriptor> vertexQuee) {
    UndirectedGraph        & gMetro        = _metro->g();
    map<VertexDescriptor, EdgeDescriptor> result;

    OutEdgeIterator e, e_end;
    for ( tie( e, e_end ) = out_edges( selectedVertex, gMetro ); e != e_end; ++e ) {
        EdgeDescriptor ed = *e;
        VertexDescriptor vj = source(ed, gMetro);
        if (vj == selectedVertex) vj = target(ed, gMetro);
        if (vertexQuee.find(vj) != vertexQuee.end()) {
            result.insert(pair<VertexDescriptor, EdgeDescriptor>(vj, ed));
        }
    }
    return result;
}

void AutoMatching::pathFromMetroLine(int l) {
    cerr << "select Line! " << l << endl;
    UndirectedGraph        & gMetro        = _metro->g();
    _foundMetroPath.erase(_foundMetroPath.begin(), _foundMetroPath.end());
    BGL_FORALL_EDGES(edge, gMetro, UndirectedGraph) {
        vector<unsigned int> lineid = gMetro[edge].lineID;
        bool contains = false;
        for (int i = 0; i < lineid.size() && !contains; i++) {
            if (lineid[i] == l) contains = true;
        }
        if (contains) {
            _foundMetroPath.push_back(source(edge, gMetro));
            _foundMetroPath.push_back(target(edge, gMetro));
            gMetro[edge].matchPath = true;
        } else {
            gMetro[edge].matchPath = false;
        }
    }
}

