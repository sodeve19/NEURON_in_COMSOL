## **Information**
Authors: Rana Jaylani, Sian McConchie

Requirements:
* MATLAB with COMSOL LiveLink
* COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)

### **Electrode Types**

* Model\_Mircowire

Creates a COMSOL model with a microwire electrode

* Model\_Shank

Creates a COMSOL model with a shank of 8 electrodes

* Model\_HighdensityShank

Creates a COMSOL model with a high density shank of 188 electrodes

* Model\_Mircowire\_Array

Creates a 9x9 array of mircowire electrodes



### **Neuron Import**
* COMSOLformatting

Transforms Gold et al. format into 2 csvs suitable for COMSOL integration 

* Import\_neuron

Imports neuron at specified location and angle into a COMSOL model
Uses:

* segment\_currents\_originalGold .csv
* segment\_morphology\_originalGold.csv



### **Spike Grid**

* Import\_spikegrid\_points

Code that imports a specified grid of point probes into a COMSOL model

* Plot\_spikegrid\_points

Takes output csv from COMSOL and plots the spike grid on the neuron



### **Biological Noise**

* Import\_BioNoise 

Imports required density of neurons far from electrode site as points

