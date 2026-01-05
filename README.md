# Migrating a Hippocampal CA1 Pyramidal Glutamatergic (GLU) Neuron from NEURON to COMSOL Multiphysics

This MATLAB code imports a hippocampal CA1 pyramidal glutamatergic (GLU) neuron [1,2] from the NEURON [3] framework into COMSOL Multiphysics. 

## Usage

1. Run the MATLAB script.
2. The code will generate a COMSOL model that can be further modified to:
   - Increase the number of neurons
   - Add background noise
   - Incorporate electrodes to record extracellular signals from the specified neuron

## Output

Running the script creates a COMSOL Multiphysics environment with the neuron embedded within a surrounding grey matter domain, ready for simulation and further modification. As shown in the video below. 

https://github.com/user-attachments/assets/348f6c75-e861-410f-8e3e-3a1f3d14b638

## References
[1] Gold, C., Henze, D. A., Koch, C. & Buzsáki, G. On the origin of the extracellular action potential waveform: A modeling
study. 95, 3113–3128, DOI: 10.1152/jn.00979.2005. Publisher: American Physiological Society.
[2] Gold, C., Henze, D. A. & Koch, C. Using extracellular action potential recordings to constrain compartmental models. 23,
39–58, DOI: 10.1007/s10827-006-0018-2.
[3] ModelDBRepository/84589. Original-date: 2019-05-29T21:06:31Z.

