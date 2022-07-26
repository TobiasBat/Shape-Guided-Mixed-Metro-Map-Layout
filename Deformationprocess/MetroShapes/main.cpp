
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
    if (argc == 5) {
        cout << "Save As Png" << endl;
        string outputName = "../data/output/frame_";
        outputName += argv[3];
        outputName += "_";
        outputName += argv[4];
        outputName += ".png";
        window.savePNG(3840, 2160, outputName);
        // app.exec();
        app.closeAllWindows();
        app.quit();
        return 0;
    }
    return app.exec();
}