To adapt the settings of the optimization process execute the program with the fourth parameter as 0. 

For example: 

```bash
./execute.sh taipei taipei-new.txt heart-guide.txt 0
```



**Default Parameter:** 

```
// Smooth
w_c: 4,00
w_l: 1,00
w_a: 2,00
w_p: 0.16
crossing Weight: 3,16

// Mixed 
w_o: 2,00
w_p: 0,10
w_c: 10,00
Crossing Weight: 14,00

// Route Matching 
cost Metro: 1,00
cost Add Edges: 3,00
weight Color: 0,300

MIN_DISTANCE: 7.0
```

For the more complex examples, `MIN_DISTANCE` has to be adapted. Change `MIN_DISTANCE` to `\#define MIN_DISTANCE (3.0)` in `Base/src/Common.h`. 



Parameter-adaptations to demonstrate the effect of the parameters:

Opt_1: 

```
// Smooth
w_c: 10,00
```



Opt_2: 

```
// Smooth
w_a: 0.5
```



Opt_3: 

```
// Smooth
w_a: 5.0
```



Opt_4: 

```
// Mixed
w_o: 5
w_c: 5
```



Opt_5:

```
// Mixed
w_o: 10
w_c: 2
```

