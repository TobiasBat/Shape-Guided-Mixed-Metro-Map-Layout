

## Requirements





## Build 

```bash
./build.sh
```



## Run

__Default Test case__

```bash
./execute.sh
```



To run the Automatic Testcases with the default parameters:

```bash
./execute.sh montreal montreal-metro.txt montreal_guide.txt 1
./execute.sh taipei taipei-new.txt heart-guide.txt 1
./execute.sh berlinBear berlin-sbahn+UBahn-Planar-multiroutes.txt bear_head_guide-closed.txt 1
./execute.sh paris-eye paris-metro-newColor-Planar-2-fixTempl.txt eye-guide-smaller.txt 1
./execute.sh tokyo tokyo-all.txt heart-guide.txt 1
```



To adapt the parameters and reproducte the manual testcases change the fourth parameter from 1 to 0. 

```bash
./execute.sh lisabon lisbon_center_manuell_planar.txt lisabon2-guide.txt 0 
./execute.sh berlin berlin-sbahn+UBahn-Planar-multiroutes.txt oval-guide.txt berlinSUOval 0
./execute.sh moscow moscow-simplified2-metro-planrized.txt circle-guide.txt 0
./execute.sh singapore-compleCircle singapore-MRT-2021-planar-circleExtension.txt circle-guide.txt 0

```

To manually select a path, click on the stations you want to add to the user-defined path. After the Optimization stage is finished, close the window.



Additional Test Cases: 

```bash
// Automatic
./execute.sh paris-cloud paris-metro-newColor-Planar-2-fixTempl.txt cloud-guide.txt 0
./execute.sh singapore-compleCircleHeart singapore-MRT-2021-planar-circleExtension.txt heart-guide.txt 0
./execute.sh berlinS-Bear berlin-sbahn-Planar.txt  bear_head_guide-closed.txt 0


// Manualy 
./execute.sh paris-circle paris-metro-newColor-Planar-2-fixTempl.txt circle-guide.txt paris-circle 0
```

