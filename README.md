# TimeDeadMapTools
Code for production and check of ITS efficiency maps.

## How to run

You need:
+ To be inside the CERN network (to reach the bookkeeping and query the API)
+ An up-to-date (22 march 2024) O2 environment


In you working directory you must have:
+ `rundeadmap.py`
+ `DeadMapQA.C`
+ `Logger.h`

(do `chmod +x rundeadmap.py` if not already executable).

To prduce the map a single coomand line is needed:

```
./rundeadmp.py <run_number>
```

Notice that the O2 workflow can take several minutes. Better to use a `tmux` session.

## The output and the log files

If everything goes well, your working directory will be populated with the following files:
```
<your working dir>
├── log.log
├── 544180.log
└── output/
    └── 544180/
        ├── full_ctf_list.dat
        ├── alien_ctf_epn291.dat
        ├── its_time_deadmap.root
        ├── mft_time_deadmap.root
        ├── o2-deadmapbuilder.err
        ├── o2-deadmapbuilder.log
        ├── orbits.png
        └──ITSQA/
           ├── DeadMapQA.log
           ├── DeadMapQA.png
           ├── DeadMapQA.root
           └── root.log
```

The log of the `rundeadmap.py` script will be in the working directory (`<run_number>.log`), while every other output file will be in `output/<run_number>/`. Please note:
+ You don't have to create any directory in advance, they will be created by the script.
+ If you re-process a run, the full directory `<run_number>/` will be deleted and overwritten!

Description of the output:
+ `its_time_deadmap.root` and `mft_time_deadmap.root` : **the time-dependent maps**. The script checks from bookkeeping that MFT was in the run. If not, only ITS workflow is run. 
+ `full_ctf_list.dat`: the list of CTFs on grid for that run
+ `alien_ctf_epn<..>.dat`: the selected CTFs, those analyzed by the workflow
+ `o2-deadmapbuilder.log/err` std out and std err of the O2 workflow
+ `orbits.png` a sketch of the orbit gap in betweemn the map steps, read from the O2 workflow logs
+ `ITSQA/` results of the root macro processing and checking the ITS object. *This directory deos not exist if the root macro fails*.
    + `DeadMapQA.log` the output printed on this file by the macro
    + `DeadMapQA.png` a summary of the quality: that's the first picture to look at, the assess the object quality
    + `DeadMapQA.png` it contains the same info of `.png`, each pad in root format
    + `root.log` that contains std err and std out of the command `root -b DeadMapQA.C`


## What to check?

The main script log (`<run_number>.log`) contains a summary of the full process. It checks the O2 workflow logs, the orbits gap in the map, the map duration compared to run duration and the presence of BAD quality spotted by the root QA macro. Therefore:

+ Just `grep` the strings "WARNING", "ERROR" or "FATAL" in the main script log! If nothing is found, most probably everything was good. 
+ For commissioning purpose, please also check that "*rundeadmap.py reached the end.*" is the last message of the main script log. This guarantees that unexpected crashes happened.

If bad quality of the QA output is reported, the details are in the last lines of `DeadMapQA.log` and in the `.png` as well. 

## More details...

to be completed 

## Extra/temporary

We are interested in finding runs where the EOR significantly exceeds the compressed data on disk because e.g. run crash. A signature for this is a missing EOX in the RCT ccdb object. To find them, look for the following message in the `root.log` file:

```
[WARN]  | Missing/invalid EOX -> use EOR
```