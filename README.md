
# TimeDeadMapTools

This repository contains the code to verify the ITS efficiency maps.
The code to produce maps out of CTFs is available in the `main` branch of the repo.

## Prerequisites

Before you begin, ensure you have the following:
- An up-to-date O2 environment (as of March 22, 2024)
- `DeadMapTREE.C`
- `DeadMapQA.py`
- `MakeCanvas.py`
- `mylogger.py`
- `Logger.h`




## Running the Scripts

`DeadMapQA.py` analyzes the map content. It uses `MakeCanvas.py` as library to prduce plots. The input is a specific ROOT file with simple trees, containing the dead map and some metadata. Such file is produced with `DeadMapTREE.C`.

Therefore, having the deadmap `.root` object on local disk, to run the analysis on it, the chain is:

```bash
root -b DeadMapTREE.C("dmap_file.root",run_number)
python3 DeadMapQA.py
```

The `.C` macro creates the file called `DeadMapTREE.root`, taken as default input by `DeadMapQA.py`. One can call the `.py` script with different options (run `python3 DeadMapQA.py --help` to get instructions).

The output `.png` files are saved into the `./canvas/` directory and the logs of the python script are saved in the local directory as `QApy.log`. 


### QA Output Files:
- `QApy.log`: The log output from the QA macro.
- `canvas/full_canvas1.png`: A summary of the ITS object quality with many pads
- `canvas/full_canvas2.png`: The time-evolution of the lane dead fraction
- `canvas/full_canvas2_zoomX.png`: The time-evolution of the lane dead fraction in each time region where the dead fraction of IB or OB was above 8%
- `canvas/full_canvas4.png`: The percentage of dead time lane by lane and stave by stave
- `canvas/full_canvas5.png`: The trends over time of the dead fraction and recovery rate for each layer



### Checks

The following checks are implemented in the `DeadMapQA.py` script

- **Invalid orbit**:
  - An invalid orbit is a map key (orbit) which is more that 330k orbits distant from the run duration 
  - `UNKNOWN` if the CTP orbit reset could not be fetched by ccdb (so that the map range wrt run duration is unknown)
  - `GOOD` if there are no invalid orbits
  - `MEDIUM` otherwise

Invalid orbits are filtered out before diaplaying the map statistics and evaluating the following checks. Information on this is still printed on `full_canvas1.png`.

- **Avg dead time IB**:
  - `GOOD` if the average dead time of IB after the first 10 seconds is below 3%
  - `MEDIUM` if it is below 10%
  - `BAD` otherwise
- **Avg dead time OB**:
  - `GOOD` if the average dead time of OB after the first 10 seconds is below 5%
  - `MEDIUM` if it is below 10%
  - `BAD` otherwise
- **Fully dead IB**:
  - `GOOD` if the number of IB chips marked as dead in every step of the map is lower than 9
  - `MEDIUM` if such number is less than 10% of the IB chips (i.e., less than 44)
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
    - `GOOD` if both the static and time-evolving maps are filled
    - `BAD` if the static map is empty
    - `FATAL` if the time-evolving map has less than 2 entries
- **Orbit gaps**:
  - `BAD` if there is at least one gap in between steps larger than 330k orbits (this is the "un-anchorable" threshold in the digitizer) or if more than 25% of the steps have a gap larger than 380 TFs = 12160 orbits
  - `MEDIUM` if there is at least one gap larger than 760 TFs or at least three gaps larger than 380 TFs. The BAD condition is evaluated with priority.
  - `GOOD` otherwise
- **Orbit range**:
  - `UNKNOWN` if the CTP orbit reset could not be fetched by ccdb (so that the map range wrt run duration is unknown)
  - `GOOD` The map's orbit range is fully contained within the run's orbit range (as stored in the RCT object), and both map boundaries lie within 3.2 seconds of the corresponding run edge
  - `MEDIUM` if the condition for GOOD is not satisfied, but each map boundary is still within 330k orbits of the run edges.
  - `BAD` otherwise (i.e. at least one map boundary lies more than 330k orbits away from the run edges).
 
