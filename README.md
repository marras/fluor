Fluor
=====

Diffusion and fluorescence simulation package. Works on Linux, Mac OS X, and probably Windows (at least used to, some time ago)

Installation
------------
Download & unzip the package. then run
$ make

Running simulation
---------------
```
$ cp fluor simulation_directory/
$ cd simulation_directory/
```

Edit configuration file: grompp.mdp

Simulate diffusion and fluorescence
```
$ ./fluor -f
```

Calculate the correlation function
```
$ ./fluor -b
```

The output appears in corr.dat .