--- 

### What is this?

This is a modifed version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat), that has been designed to make polarization calibration as hands-off as possible. I direct the reader to the original oxkat documentation for a more detailed description of the pipeline. Here I will highlight some key differences: 

---
##### Overarching Changes

* This uses a mix of containers, namely `oxkat-0.41.sif` and `polkat-0.0.1.sif`. The latter is only on the idia cluster, and, as a result, if you want to use polkat on other clusters (or locally) you have to download the container and put it in the appropriate directories.

* There are a lot of similarites between the two containers, however, `polkat-0.0.1.sif` now uses [`casa`](https://casa.nrao.edu/) version 6.x. and I have the intention of (eventually) changing all of the [`cubical`](https://github.com/ratt-ru/CubiCal) scripts into [`quartical`](https://github.com/ratt-ru/CubiCal) script. 


---
##### Changes of not with respect to oxkat

* The ms file is averaged in the ``INFO.py''
