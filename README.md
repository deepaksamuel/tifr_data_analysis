This code estimates fit parameters from the data collected from the TIFR stack. 

Compile using:
``g++ -o analysis ino_analysis.cpp main.cpp `root-config --glibs --cflags ``

Run using:
`./analysis.out`

A sample file is given for testing: `INORUN_20160819_053749.ire`

The default output file is `out.fit`

This code can also convert TIFR raw data from ROOT format (.ire) to ASCII.
Only the Event number, timestamp and Hits are pushed to the output
For this, in main.cpp, disable the line:
```// #include "ino_analysis.h"```
and enable the line
```#include "convert_ascii.h" ```
Then disable the following line:
```  //ino_analysis *a = new ino_analysis(T);```
and enable the following line:
```convert_ascii *a = new convert_ascii(T); ```

The compile command is now:
``g++ -o analysis convert_ascii.cpp main.cpp `root-config --glibs --cflags ``
The output is by default written in out.csv.
The first column in eventID, second is timestamp.
Only the stripIDs fired in an event are then written in sequence. The stripIDs are generated as follows
- x side layer 0 strip id will be from 0 -31
- x side layer 1 strip id will be from 32 -63 
- max strip id on xside will be 383 for 12 layers and 32 strips
- y side layer 0 strip id will be from 384 -415
- y side layer 1 strip id will be from 416 -447 

For instance, the line
100, 2, 32, 417 
indicates that 
- the eventID was 100
- timestamp was 2 
- xside layer 1 strip 0 was fired
- yside layer 1 strip 1 was fired

101, 2, 32, 417, 0
indicates that 
- the eventID was 101
- timestamp was 2 
- xside layer 1 strip 0 was fired
- yside layer 1 strip 1 was fired
- xside layer 0 strip 0 was fired
