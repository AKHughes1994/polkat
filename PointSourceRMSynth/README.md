### Here is a 1D RM Synthesis Script that is adopted from the CIRADA

* Original repo can be found [here](https://github.com/CIRADA-Tools/RM-Tools). The  only change is that the noise is scaled up using the RMS rather than the MAD (which is the default). From past expeirence with this code, the RMS is normally a factor of ~2 higher, and will make results more convincing. 

* For those unfamiliar with RM Synthesis it the classic paper to read is [(Brentjens & de Bruyn 2005)](https://arxiv.org/abs/astro-ph/0507349), which details the methods used in the attached code. Note that this method coherently adds the polarization vectors from each spectral channel coherently. Therefore, the detection threshold is set by the integrated (i.e., band-averaged S/N) and not the signal to noise of the detections in each spectral channels. Exploratory studies have shown that for band-averaged S/N > 7 RM synthesis is reliable [(e.g., Macquart et al. 2012)](https://arxiv.org/abs/1203.2706). Don't use RM synthesis for 5sigma detections (you likely just will not be able to get reliable RM measures at these S/N)

* Interestingly, RM Synthesis is a direct (1-dimensional) analog to synthesis imaging with radio interferometers. The polarization intensity at different frequencies, can be fourier transformed into a "Faraday Dispersion Function" (FDF) at different "Faraday Depths". Due to incomplete sampling the intial FDF has artifacts is dirty, and therefore, is CLEANed. For many sources the CLEAN FDF will have a single pronounced peak, where the value of the peak is related to the de-rotated Polarization Intensity (i.e., absent any bandwidth depolarization), the ratio of the real and imaginary components at the Peak is the polariation angle, and the Faraday depth that corresponds to the peak is the rotation measure.

Below is a somewhat detailed example.

---

The first step requires a txt file in the appropriate format. There is an example in `datafiles`. The format is: freq (Hz), I_flux (Jy), Q_flux (Jy), U_flux (Jy), I_rms (Jy), Q_rms (Jy), U_rms (Jy). Note, that `polkat/tools` contains a script called `extract_rm_synth.py`. Following a sucessful polkat run, within the home directory run, for example,  

```
singularity exec /idia/software/containers/polkat-0.0.1.sif casa --nologger -c tools/extract_rm_synth.py -r "circle[[17:27:43.41,-16.12.20.37],15arcsec]" -p "SwiftJ1727_rmsynth" -i "IMAGES/img_1696160474_sdp_l0_1024ch.ms_SwiftJ1727_pcalmask"
```

and the necessary file is produced.

From here RM Synthesis is broken into two steps

First the initial creation (or "synthesis") of the dirty FDF:

```
python3 ../RMtools_1D/do_RMsynth_1D.py ../datafiles/SwiftJ1727_rmsynth_60218.48831018519.txt -S -v -p -o 4
```

Here the `-S` means that you save all of the dirty files, `-v` is the verbosity flag which outputs the values to the terminal, `-p` toggles on the plots, `-o` is the order of the first to the Stokes I spectra (it doesn't have a large effect as long as the spectra is "wacky").  There are other useful inputs that I am not usinghere; in the past I've used `-d` which manually sets the widths of the Faraday depth bins (units of rad/m^2), and `-l` which specifies the maximum (+/-) Faraday depth. This is analagous to changing the cellsize and image dimensions during synthesis imaging. 

Following make the dirty FDF you clean it: 

```
python3 ../RMtools_1D/do_RMclean_1D.py ../datafiles/SwiftJ1727_rmsynth_60218.48831018519.txt -v -p
```

From here you will have your final data products. Below I've included an example output which is a CLEAN FDF. 

![image](https://github.com/AKHughes1994/XKATPol/assets/49698839/914c8ded-cdfa-4fe7-8e7a-78665b8701a5)





