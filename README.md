# RADISH: Rabi Distance Search

Welcome to our repository for RADISH.

These simple MATLAB scripts will allow you to process WASABI data into $\delta\omega$ and $B_1$ maps within seconds.

Please view **radish_demo.m** for demonstration on how to set up your parameters and run RADISH.

A manuscript has been submitted to MRM, pending reviews.

RADISH is based on the original work by Schuenke et al. (2017)

* P. Schuenke, J. Windschuh, V. Roeloffs, M. E. Ladd, P. Bachert, M. Zaiss. *Simultaneous mapping of water shift and B1(WASABI)—Application to field-Inhomogeneity correction of CEST MRI data.* Magnetic Resonance in Medicine. 2017;77(2):571–580.

If you have any issues running RADISH, or if you have any feedback/suggestions, please contact Mara Quach at trinhq@student.unimelb.edu.au

## Trouble-shooting
1. If radish.m runs successfully without a parfor loop, but fails with one, consider changing the parallel environment from "processes" to "threads". 
