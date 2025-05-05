# TSP

TSP solver for the Operations Research 2 course.
By: Candau Jaime, Dario Alessandro, Zerbinati Riccardo.

## Build

Relase build using CMake:

```
mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCPLEX_ROOT_DIR=/path/to/cplex/root .. && make
```

## Performance Profile 

Run performance profile on a .csv file 

```
python src/perfprof.py -D , -T 3600 -S 2 -M 20 perf_prof.csv pp.pdf -P “Name_Graph” -X "ratio"
```