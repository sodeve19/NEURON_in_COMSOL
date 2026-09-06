# Migrating a Hippocampal CA1 Pyramidal Glutamatergic (GLU) Neuron from NEURON to COMSOL Multiphysics
This repository provides MATLAB scripts for transferring a hippocampal CA1 pyramidal glutamatergic (GLU) neuron model from the NEURON simulation environment into COMSOL Multiphysics.
The neuron morphology and compartmental model data are based on the work of Gold et al. [1,2]. The MATLAB scripts preprocess the model data into a format that can be imported into COMSOL, where the neuron can be combined with extracellular electrode configurations and a surrounding grey-matter domain.

The repository also includes scripts for generating electrode geometries, importing background biological noise, and creating grids of point probes for recording extracellular signals.


## Software requirements

1. MATLAB with COMSOL Livelink
3. COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)

## Repository Structure

### **Electrode Models**
The following MATLAB scripts generate the base COMSOL models containing different electrode configurations:

* Model\_Mircowire
Creates a COMSOL model with a microwire electrode

* Model\_Shank
Creates a COMSOL model with a shank of 8 electrodes

* Model\_HighdensityShank
Creates a COMSOL model with a high density shank of 188 electrodes

* Model\_Mircowire\_Array
Creates a 9x9 array of mircowire electrodes

These scripts are used to establish the electrode geometry and surrounding computational domain before importing the neuron.

### **Neuron Import**
The neuron-import workflow converts the compartmental neuron data from the Gold et al. model into a format that can be used to construct the neuron within COMSOL.

* COMSOLformatting
Preprocesses the original Gold et al. model data and generates two CSV files required by the neuron-import workflow:
- segment_morphology_originalGold.csv
- segment_currents_originalGold.csv

* Import\_neuron
Imports the processed neuron morphology and electrical data into an existing COMSOL model.
The neuron can be positioned at a specified location and orientation within the COMSOL geometry. The script uses:
- segment\_currents\_originalGold .csv
- segment\_morphology\_originalGold.csv

### **Spike Grid**

* Import\_spikegrid\_points
Code that imports a specified grid of point probes into a COMSOL model

* Plot\_spikegrid\_points
Takes output csv from COMSOL and plots the spike grid on the neuron

### **Biological Noise**

* Import\_BioNoise 
Imports required density of neurons far from electrode site as points

## Workflow

A typical simulation can be generated using the following workflow:

1. Create the base electrode model: Run one of the Model_* MATLAB scripts to generate the desired COMSOL electrode configuration and surrounding domain.
2. Prepare the neuron data: Run COMSOLformatting to convert the original Gold et al. model data into the CSV files required by COMSOL.
3. Import the neuron: Run Import_neuron to construct the CA1 pyramidal neuron within the COMSOL model at the desired position and orientation.
4. Add additional model components: The resulting COMSOL model can then be modified to include:
- additional neurons
- background biological noise
- spike-grid point probes
5. Run the COMSOL simulation

## Output

Running the workflow produces a COMSOL Multiphysics model containing the imported CA1 pyramidal neuron within a surrounding grey-matter/extracellular domain.

Depending on the selected model configuration, the simulation can also contain electrode geometries, additional background neural activity, and spatially distributed point probes for measuring extracellular potentials.

An example of the resulting model and workflow is shown below:
https://github.com/user-attachments/assets/348f6c75-e861-410f-8e3e-3a1f3d14b638

## References
[1] Gold, C., Henze, D. A., Koch, C. & Buzsáki, G. On the origin of the extracellular action potential waveform: A modeling
study. 95, 3113–3128, DOI: 10.1152/jn.00979.2005. Publisher: American Physiological Society.

[2] Gold, C., Henze, D. A. & Koch, C. Using extracellular action potential recordings to constrain compartmental models. 23,
39–58, DOI: 10.1007/s10827-006-0018-2.

[3] ModelDBRepository/84589. Original-date: 2019-05-29T21:06:31Z.

