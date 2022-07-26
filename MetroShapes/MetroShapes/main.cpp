
#include <string>
#include <iostream>


#include "MainWindow.h"
#include "core/MetroShapes.h"

int main (int argc, char **argv) 
{    
    QApplication app(argc, argv); 
    Ui::MainWindow window; 

    MetroShapes metroShapes; 
    metroShapes.init( argc, argv );
    // metroShapes.computeMatching();  //TODO Should be called by button
    Base *b_ptr = &metroShapes;
	window.init( &metroShapes.getMetro(), b_ptr, &metroShapes._smooth, &metroShapes._mixedlayout ,&metroShapes.getGuid(), &metroShapes.getMatching(), &metroShapes.getAutoMatching());
    window.show();
    return app.exec();
}