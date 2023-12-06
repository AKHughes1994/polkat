--- 

### What is this?

This is a modifed version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat), that has been designed to make polarization calibration as hands-off as possible. I direct the reader to the original oxkat documentation for a more detailed description of the pipeline. Here I will highlight some key differences.


---
##### Overarching Changes

* This uses a mix of containers, namely `oxkat-0.41.sif` and `polkat-0.0.1.sif`. The latter is only on the idia cluster, and, as a result, if you want to use polkat on other clusters (or locally) you have to download the container and put it in the appropriate directories.

* There are a lot of similarites between the two containers, however, `polkat-0.0.1.sif` now uses [`casa`](https://casa.nrao.edu/) version 6.x. and I have the intention of (eventually) changing all of the [`cubical`](https://github.com/ratt-ru/CubiCal) scripts into [`quartical`](https://github.com/ratt-ru/CubiCal) script. 


---
##### Changes of not with respect to oxkat

* The ms file is averaged to 1024 channels in the INFO step (oxkat does this in 1GC. The reason for this has to do with getting the numeric indexes of the (sub-)fields that are used in each calibration run. INFO also identifies the polarization calibator automatically. Currently, MeerKAT defaults to specifying the scan intent as "UNKOWN" for the polarization calibrator.

* 1GC is the biggest change to the calibration routine:
  1. We are now solving for additional leakage (Df) terms using the (unpolarized) primary calibrator, Cross-hand delay (KCROSS) terms using the polarization calibrator, and Cross-hand phase (Xf) terms using the polarization calibrator.
  2. The Stokes I models for the secondary calibrator, and the polarization calibrator are solved for using interative imaging and self-calibration (Step name is `RECAL` and uses the `1GC_casa_refinement.py` script).
  3. Ampltiude solutions are transferred from the primary calibrator (not the secondary), any temporal evolution is mapped onto the source through an additional gaintype "T" calibration table that solves for a single temporal evoltion for both the X and Y. This is to avoid transfering polarized signal from the (polarized) seconday onto the source due to differences betweenteen the X and Y gain solutions from the secondary. 
  4. The calibration routine is now broken into three steps: (largely) parallel-hand XX,YY with `1GC_casa_refcal.py`, secondary/polarization model refinemnet with `1GC_casa_refinement.py`, and (mostly) cross-hand with `1GC_casa_polcal.py`
 
* At the end of 1GC, by defualt, you will image the seconday in Stokes IQUV. You can use this to check the calibration against the [cataloged values](https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1452146701/L-band+gain+calibrators)
