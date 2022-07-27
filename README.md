

## Requirements

+ cmake 3.10+
+ qt 5.13
+ Boost Graph Library 1.71
+ Eigen3 3.3 



## Build 

```bash
./build.sh
```



## Run

Run with the default test case:

```bash
./execute.sh
```



Parameters can be added to run the program with other metro networks: 

```bash
 ./execute.sh <test-case-name> <metronetwork>.txt <guide-name>.txt <compute automatically>
```



To run the Automatic Testcases with the default parameters:

```bash
./execute.sh montreal montreal-metro.txt montreal_guide.txt 1
./execute.sh taipei taipei-new.txt heart-guide.txt 1
./execute.sh paris-eye paris-metro-newColor-Planar-2-fixTempl.txt eye-guide-smaller.txt 1
```



To adapt the parameters and reproduce the manual testcases change the fourth parameter from 1 to 0. 

```bash
./execute.sh lisabon lisbon_center_manuell_planar.txt lisabon2-guide.txt 0 
./execute.sh berlin berlin-sbahn+UBahn-Planar-multiroutes.txt oval-guide.txt berlinSUOval 0
```

To manually select a path, click on the stations you want to add to the user-defined path. After the Optimization stage is finished, close the window.





