# SciFiTestbeamAndSimulation
To compile the code, create a build directory and run cmake:
```
mkdir build
cd build
cmake ..
```

Then you can run
```
make <the program I want to compile>
```
where the possible programs are present in the main folder:
```
produceCorrectedFile : produces a new file that contains the pedestal and gain corrected data
clusterAnalysis : searches for clusters in a corrected data file
```
To run the compiled program from the build folder, run
```
./bin/the_program
```

Available options for `produceCorrectedFile` are:
```
--file2correct, -f : test beam data file containing the data run
--umax, -u",  : uplink number to stop at [1, ... , 8], default_value(4)
--umin, -l", : "uplink number to start at [1, ... , 8], default_value(3)
```

Available options for `clusterAnalysis` are:
```
--file, -f : corrected test beam data file(s), wild cards like *.root are also supported
--simulation, -s : add this optiion if you are running on a simulated file
--clusteralg, -c : clustering algorithm: b for Boole or m for Maxs ,default_value("b")
--tag, -t : tag that is added to the output file name, default_value("")
```
