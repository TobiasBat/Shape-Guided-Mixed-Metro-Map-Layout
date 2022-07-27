

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
```



To adapt the parameters and reproducte the manual testcases change the fourth parameter from 1 to 0. 

```bash
./execute.sh lisabon lisbon_center_manuell_planar.txt lisabon2-guide.txt 0 
./execute.sh berlin berlin-sbahn+UBahn-Planar-multiroutes.txt oval-guide.txt berlinSUOval 0

./execute.sh paris-cloud paris-metro-newColor-Planar-2-fixTempl.txt cloud-guide.txt 0

```

After the Optimization stage is finished, close the window.



