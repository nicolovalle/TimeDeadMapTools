
# TimeDeadMapTools

This repository contains the code to generate and verify ITS efficiency maps.

## Prerequisites

Before you begin, ensure you have the following:
- A valid bookkeeping token
- An up-to-date O2 environment (as of March 22, 2024)

Place your token in a file named `token.dat`.

The following files should be present in your working directory:
- `rundeadmap.py`
- `mylogger.py`
- `DeadMapQA.C`
- `Logger.h`
- `token.dat`

Make the `rundeadmap.py` executable, if necessary:
```bash
chmod +x rundeadmap.py
```

## Running the Script

To generate the map, jyst use the following command:
```bash
./rundeadmap.py <run_number>
```
Please note that the O2 workflow can take several minutes to complete. It's recommended to use a `tmux` session to avoid interruptions.

## Output and Log Files

Upon successful execution, your working directory will contain the following:

```
<working dir>
├── log.log
└── output/
     └── <run_number>/
         ├── main.log
         ├── period.txt
         ├── run.json
         ├── full_ctf_list.dat
         ├── alien_ctf_epn<epn>.dat
         ├── its_time_deadmap.root
         ├── mft_time_deadmap.root
         ├── o2-deadmapbuilder.err
         ├── o2-deadmapbuilder.log
         ├── orbits.png
         └── ITSQA/
             ├── DeadMapQA.log
             ├── DeadMapQA1.png
             ├── DeadMapQA2.png
             ├── DeadMapQA3.png
             ├── DeadMapQA4.png
             ├── DeadMapQA5.png
             └── root.log
```

**Important notes:**
- You do not need to create any directories manually; the script handles that for you.
- If you reprocess a run, the directory `<run_number>/` will be deleted and overwritten.

### Output File Descriptions
- `main.log`: Log of the main script `rundeadmap.py`, which moves to the output folder after completion.
- `period.txt`: Contains the LHC period of the run.
- `run.json`: Information about the run, fetched from the bookkeeping system.
- `its_time_deadmap.root` and `mft_time_deadmap.root`: The generated time-dependent maps.
- `full_ctf_list.dat`: A list of CTFs available on the grid for the run.
- `alien_ctf_epn<epn>.dat`: The selected CTFs analyzed by the workflow.
- `o2-deadmapbuilder.log/err`: Standard output and error logs from the O2 workflow.
- `orbits.png`: A visualization of the orbit gap, as read from the workflow logs.
- `ITSQA/`: Output from the quality assessment (QA) checks of the ITS object. This directory is not created if the root macro fails.

### QA Output Files:
- `DeadMapQA.log`: The log output from the QA macro.
- `DeadMapQA1.png`: A summary of the ITS object quality.
- `DeadMapQA2.png`: The lane status history versus orbits.
- `DeadMapQA3.png`: The average dead time, stave by stave.
- `DeadMapQA4.png`: The average dead time, lane by lane.
- `DeadMapQA5.png`: The average time evolution of dead time for each layer.
- `root.log`: Standard output and error logs from the command `root -b DeadMapQA.C`.

## What to check

The `main.log` file provides a summary of the process, including checks for the O2 workflow logs, orbit gaps in the map, and run duration versus map duration. It also flags any bad quality detected by the QA macro. 

You can check for issues by searching for "WARN", "ERROR", or "FATAL" in the main script logs. If no such strings are found, everything likely ran successfully.

If a QA issue is flagged, check the final lines of `DeadMapQA.log` and inspect the `.png` files for more details.

### Automatic checks

The following automatic checks are implemente in the `DeadMapQA.C` macro.

- **Avg dead time IB**:
  - `GOOD` if the average dead time of IB after the first 10 seconds is below 3%
  - `MEDIUM` if it is below 10%
  - `BAD` otherwise
- **Avg dead time OB**:
  - `GOOD` if the averga dead time of OB after the first 10 seconds is below 5%
  - `MEDIUM` if it is below 10%
  - `BAD` otherwise
- **Chip interval**:
  - `GOOD` if the chip IDs in the time-evolving part of the map are properly grouped into lanes
  - `FATAL` otherwise
- **Fully dead IB**:
  - `GOOD` if the number of IB chips marked as dead in every step of the map is lower than 9
  - `MEDIUM` if such number is less than 10% of the IB chips (i.e. less than 44)
  - `BAD` otherwise
- **Fully dead OB**:
  - `GOOD` if the number of OB lanes with at least one chip which is always dead is lower than 68 (roughly 2% of the lanes)
  - `BAD` otherwise
- **Default object**:
  - `FATAL` if the object is the default one
  - Not declared otherwise
- **Map size**:
  - If the object is the default one:
    - `GOOD` if both the static and time-evolving maps are empty
    - `FATAL` otherwise
  - If the object is not the default one:
    - `GOOD` if both the statis and time-evolving maps have entries
    - `BAD` if the static map is empty
    - `FATAL` if the time-evolving map is empty
- **Null orbit**:
  - `GOOD` if orbit = 0 is not among the map keys
  - `MEDIUM` if the map contains orbit = 0 and all the chips at that step are marked as bad
  - `BAD` if the map contains orbit = 0 with alive chips, or if the map contains also negative orbits
- **Orbit gaps**:
  - `BAD` if there is at least one gap in between steps above 330k orbits (this is the "un-anchorable" threshold in the digitizer) or if more than 25% of the steps have gap above 380 TFs = 12160 orbits
  - `MEDIUM` if there is at least one gap above 760 TFs or at least three gaps above 380 TFs. The BAD condition is evaluated first.
  - `GOOD` otherwise
- **Orbit range**:
  - `GOOD` if the time interval covered by the map is within 5 seconds of the run duration as specified by the RCT ccdb object
  - `MEDIUM` if the map's time interval exceeds the run duration by more than 5 seconds or is shorter by no more than 30 seconds
  - `BAD` otherwise (i.e. the difference between the run duration and the map duration is larger than 30 seconds).
- **Un-anchorable fraction**:
  - This is the fraction of orbits within the total map range that are not anchorable to a map element by the digitizer. The digitizer rejects the map if the nearest orbit is more than 330k orbits away from the requested one.
     - `GOOD` if the fraction is below 2% (consider than a single gap above 330k orbits is not expected, and it trigger a BAD "Orbit gaps" check).
     - `MEDIUM` if the fraction is below 5%
     - `BAD` otherwise
