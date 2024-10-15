
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
