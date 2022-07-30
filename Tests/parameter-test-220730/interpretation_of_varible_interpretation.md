# Interpretation of different Parameters



To adapt the settings of the optimization process execute the program with the fourth parameter as 0. 

For example: 

```bash
./execute.sh taipei taipei-new.txt heart-guide.txt 0
./execute.sh moscow moscow-simplified2-metro-planrized.txt circle-guide.txt 0 
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



Opt_6: 

```
// Mixed
w_o: 10
```



Opt_7: 

```
// Mixed: 
w_c: 2
```



----



## Interpretation



**Smooth Layout Varibles (Testcase 1-3)** In the smooth figure of the Taipei system with the default parameters, you can see that the metro lines are straightened and stations are evenly spaced apart. Options 2 illustrated the effect of $w_a$ clearly. By decreasing the value, the system prioritizes the other objective functions and metro lines are less straight. Resulting in a less clear layout after the mixed stage. By increasing $w_a$ (Opt_3) metro lines are even more straightened at the cost of the other objective function. Resulting in a layout that does not preserve the geographic shape of the metro network. For the more complex Moscow metro system, the effects of the different parameters are even greater. Indicating the potential conflict between the different objective functions and that the system is not able to find a reasonable solution. For example, if $w_a$ is decreased the approach does not produce good results. 



**Mixed Layout Varibles (Testcases 4-7)** For the Taipei-test data, the variation of the mixed weights $w_o$ and $w_c$ have relatively little effect on the results of the final layout. Indicating that the process is relatively stable and the system can find a good solution for each objective function. For the Moscow-test case, the system fails to rotate a greater amount of edges to an octolinear angle independent from the parameters. Especially in the dense center of the map. Therefore when the weight $w_c$ is decreased (Opt_5), the system prioritizes the colinearity of edges and pushes stations away from the guide that should be located close to the guide shape. 

We compare the different results of the smaller, less complex metro system compared to the larger more complex Moscow system it can be seen that in the complex systems the weights have to be balanced carefully. While in simple networks the system produces good results for a variety of parameters and the system finds a solution that satisfys the different objective functions. 