# TimeDeadMapTools
Code for production and check of ITS efficiency maps.

## How to run

You need:
+ A bookkeeping token
+ An up-to-date (22 march 2024) O2 environment

Write your token in a file named `token.dat`.

In you working directory you must have:
+ `rundeadmap.py`
+ `mylogger.py`
+ `DeadMapQA.C`
+ `Logger.h`
+ `token.dat`

(do `chmod +x rundeadmap.py` if not already executable).


To prduce the map a single coomand line is needed:

```
./rundeadmp.py <run_number>
```

Notice that the O2 workflow can take several minutes. Better to use a `tmux` session.

## The output and the log files

If everything goes well, your working directory will be populated with the following files:
```
<working dir>
├── log.log
└── output/
     └── 544180/
         ├── main.log
	 ├── period.txt
	 ├── run.json
         ├── full_ctf_list.dat
         ├── alien_ctf_epn291.dat
         ├── its_time_deadmap.root
         ├── mft_time_deadmap.root
         ├── o2-deadmapbuilder.err
         ├── o2-deadmapbuilder.log
         ├── orbits.png
         └──ITSQA/
             ├── DeadMapQA.log
             ├── DeadMapQA1.png
             ├── DeadMapQA2.png
             ├── DeadMapQA3.png
             ├── DeadMapQA4.png
             ├── DeadMapQA5.png
             └── root.log
```

Please note:
+ You don't have to create any directory in advance, they will be created by the script.
+ If you re-process a run, the full directory `<run_number>/` will be deleted and overwritten!

Description of the output:
+ `main.log` is the log of the main script `rundeadmap.py`. When the script is being executed, the log is written in the working directory and named `<run_number>.log`. It's moved in the output folder as the very last action of the script.
+ `period.txt`. It contains a single string with the LHC period of the run.
+ `run.json`. The run info read from bookkeeping.
+ `its_time_deadmap.root` and `mft_time_deadmap.root` : **the time-dependent maps**. The script checks from bookkeeping that MFT was in the run. If not, only ITS workflow is run. 
+ `full_ctf_list.dat`: the list of CTFs on grid for that run
+ `alien_ctf_epn<..>.dat`: the selected CTFs, those analyzed by the workflow
+ `o2-deadmapbuilder.log/err` std out and std err of the O2 workflow
+ `orbits.png` a sketch of the orbit gap in betweemn the map steps, read from the O2 workflow logs
+ `ITSQA/` results of the root macro processing and checking the ITS object. *This directory deos not exist if the root macro fails*.
    + `DeadMapQA.log` the output printed on this file by the macro
    + `DeadMapQA1.png` a summary of the quality: that's the first picture to look at, to assess the object quality
    + `DeadMapQA2.png` the complete history of the lanes status vs orbit
    + `DeadMapQA3.png` the average dead time, stave by stave
    + `DeadMapQA4.png` the average dead time, lane by lane
    + `DeadMapQA5.png` with the averaged time evolution of dead time of each layer
    + `root.log` contains std err and std out of the command `root -b DeadMapQA.C`


## What to check?

The main script log contains a summary of the full process. It checks the O2 workflow logs, the orbits gap in the map, the map duration compared to run duration and the presence of BAD quality spotted by the root QA macro. Therefore:

+ Just `grep` the strings "WARN", "ERROR" or "FATAL" in the main script logs! If nothing is found, most probably everything is good.

None of the main script logs should remain in the working dir (`<run_number>.log`). If this happens, the script crashed or failed in an unexpected way.


If bad quality of the QA output is reported, the details are in the last lines of `DeadMapQA.log` and in the `.png` as well. 

## More details...

to be completed 
