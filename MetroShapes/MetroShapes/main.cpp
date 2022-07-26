
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
    Base *b_ptr = &metroShapes;
	window.init( &metroShapes.getMetro(), b_ptr, &metroShapes._smooth, &metroShapes._mixedlayout ,&metroShapes.getGuid(), &metroShapes.getMatching(), &metroShapes.getAutoMatching());
    window.show();

    if (argc == 5 ) {
        string argv_4 = argv[4];
        if (argv_4.compare("1") == 0) {
            string outputname = "../../../output-optimisation/";
            outputname += argv[3];
            outputname += "/output_image.svg";
            window.saveSVG(1200, 1200, outputname);
            app.closeAllWindows();
            app.quit();
            return 0;
        }
    }

    return app.exec();
}