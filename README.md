This code estimates fit parameters from the data collected from the TIFR stack. 

Compile using:
``g++ -o analysis ino_analysis.cpp main.cpp `root-config --glibs --cflags ``

Run using:
`./analysis.out`

A sample file is given for testing: `INORUN_20160819_053749.ire`

The output file is `out.fit`