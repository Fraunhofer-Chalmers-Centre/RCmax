## The console interface
RCmax.exe is an executable to run the LR algorithm on an input or generated instance. (built with clang for windows)

RCmax and RCmax.out are the corresponding binaries for linux systems, built with gcc and clang, respectively. 

Invoke as: RCmax.exe path_to_instance.dat time_limit_in_seconds flag1 flag2...

Potential flags are (sa. sensitivity analysis in paper):

        -nV             not use variable fixing
        -nS             use few subgradient iterations
        -nD             not use deflected subgradient
        -nL             not use large neighbourhood search
        -nE             not use ergodic sequence heuristic
        -o path_to_result_file
                        direct the result to the output_file
        -v              print progress
        -g m n t        generate a new instance with m machines and n tasks and write to the instance path, where t = u,j,m, (unrelated, correlated jobs, correlated machines)


The instances are assumed to be comma separated matrices of the processing times p_ij, e.g., an instance with 3 machines and 5 tasks can specified as:

```
66,91,15,56,86
81,75,18,55,78
81,86,14,55,76
```

use the -g flag to generate instances on this format.

## Building from source:
On linux build using cmake (we observed slightly better performance using clang, remove first line if you prefer gcc)
```
export CXX=/usr/bin/clang++
cmake -DCMAKE_BUILD_TYPE=Release
make
```
On windows use visual studio and open as a cmake project, again prefer using the `clang compiler for windows addon`

## Citing RCmax 

If you find RCmax useful in your work, we kindly request that you cite the following paper:

```
@article{ablad2021rcmax,
     author = {Edvin {\AA}blad and Domenico Spensieri and Ann-Brith Str{\"o}mberg},
     title = {Exact makespan minimization of unrelated parallel machines},
     journal = {Open Journal of Mathematical Optimization},
     eid = {4},
     publisher = {Universit\'e de Montpellier},
     volume = {2},
     year = {2021},
     doi = {10.5802/ojmo.4},
     language = {en},
     url = {https://ojmo.centre-mersenne.org/articles/10.5802/ojmo.4/}
}
```

