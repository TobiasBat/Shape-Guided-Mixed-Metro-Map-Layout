# MetroShapes: Embedding User-Defined Shapes into Metro Map Layouts

Thesis: [https://www.cg.tuwien.ac.at/research/publications/2021/Batik-2021/](https://www.cg.tuwien.ac.at/research/publications/2021/Batik-2021/)



### Prerequisites 

To compile the code, the following libraries have to be installed on the computer: 

+ cmake 3.10+
+ qt 5.13
+ Boost Graph Library 1.71
+ Eigen3 3.3 



### Compilation

To compile the program, clone the repository into a new directory and execute the following commands: 

```
$ git clone MetroShapes
$ cd MetroShapes
$ mkdir build
$ cd build
$ cmake ../
$ make
$ cd bin
$ ./MetroShapes
```



### Using

The program can be started by the following command: 

```
$ ./MetroShapes <metro-network>.txt <guide-name>.txt
```

To switsch between automatic scenario and manual one change in Smooth.h
```c++
    bool            _exclusivlyPathMode = true;
```

To export the GraphML Files, go to the edit menu or press cmd + O. 
To create an svg press cmd + s

Testcases for the Vis Paper 

```
Auto: 
./MetroShapes taipei-new.txt heart-guide.txt taipei_hearth
./MetroShapes montreal-metro.txt montreal_guide.txt montreal
./MetroShapes berlin-sbahn+UBahn-Planar-multiroutes.txt bear_head_guide-closed.txt
./MetroShapes paris-metro-newColor-Planar-2-fixTempl.txt eye-guide-smaller.txt paris-eye
./MetroShapes tokyo-all.txt heart-guide.txt tokyo-all


Manuelle
./MetroShapes lisbon_center_manuell_planar.txt lisabon2-guide.txt lisabon
./MetroShapes singapore-MRT-2021-planar-circleExtension.txt circle-guide.txt singapore-compleCircle
./MetroShapes moscow-simplified2-metro-planrized.txt circle-guide.txt moscow
./MetroShapes berlin-sbahn+UBahn-Planar-multiroutes.txt oval-guide.txt berlinSUOval

Additional one for survey
Auto: 
./MetroShapes paris-metro-newColor-Planar-2-fixTempl.txt cloud-guide.txt paris-cloud
./MetroShapes singapore-MRT-2021-planar-circleExtension.txt heart-guide.txt singapore-compleCircleHeart
./MetroShapes berlin-sbahn-Planar.txt  bear_head_guide-closed.txt berlinS-Bear

Manuelle:
./MetroShapes paris-metro-newColor-Planar-2-fixTempl.txt circle-guide.txt paris-circle
```




### Interface 

A new matching between a metro station and a guide vertex can be defined by first clicking onto a guide vertex, and then on a metro station. 
By clicking onto the button "Match" the guide is transformed. 
The buttons "Smooth" and "Mixedlayout" are used to generate the mixed/smooth layout. 
By clicking the button "Recalculate" the system transforms the guide, calculates a smooth layout, and after a mixed layout. 
On the left sidebar, it is possible to define which metro line is affected by the guide path. 



### Tested environment

The framework has been developed and tested under MacOs system. 


### Export Frames
requires adaptation of code
```python
import subprocess

for i in range(0,10):
    print('calling ', i)
    subprocess.call(['/Users/tobiasbatik/Documents/TU/8-Semester/05_MetroProject/_code/MetroShapes/build/bin/MetroShapes', 'paris-metro-newColor.txt', 'circle-guide.txt', 'paris_circle', str(i)])

```

----



The code developed in this project builds upon the code of Hsiang-Yun Wu.


