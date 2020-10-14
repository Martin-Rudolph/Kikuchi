# An experimental code for the partial reconstruction of a Kikuchi sphere from recorded EBSD patterns
_Matlab 2018b_


**Peter Fischer, [Martin Rudolph](mailto:m.s.rudolph@outlook.com), Andreas Leineweber, David Rafaja, Institute of Materials Science, Gustav-Zeuner-Straße 5, TU Bergakademie Freiberg**

## 1 Installation instructions
Install a valid version of Matlab 2018b or higher. Download the folder MergeEBSD and run within this folder the script MergeEBSD.m.

## 2 Run the script
For a first example run no changes are required, just press the run button or F5.
However, several options and parameters can be changed in the **Input** section.
Brief parameter descriptions are given in the script as comment.
It is recommended to use as reconstruction type _overlayed_ (cut works not properly) and the others may result in blurry images in case of imperfect data.

### 2.1 Simulated patterns
The following images were simulated (DynamicS, Esprit, Bruker) for a specimen were the crystal coordinate system matches the specimen coordinate system (a||S1, b||S2, c||S3). Pattern center (PC) was for all images [0.5, 0.25] and the relative distande between screen and specimen (DD) was 0.6 (see images below PC-X, PC-Y, DD).

First image (reference phi = 0°):

![Simulated_0](/Pictures/Test_0.PNG)

Due to the specific orientation and the 4-fold roataion axis (c-axis) the same image is obtained for phi = 90°, 180° and 270°

Second image (reference phi = 30°)

![Simulated_30](/Pictures/Test_30.PNG)

Third image (reference phi = 60°)

![Simulated_60](/Pictures/Test_60.PNG)

For each pattern the reconstruction was carried out as depicted below:

![Simulated_30](/Pictures/ProjectionSimulated.PNG)

where S is the specimen coordinate (red) system, L is the laboratory coordinate system (green) and the red point marks the pattern center.

Superimposing the single images lead to the reconstruction of the Kikuchi sphere:

![Simulated_30](/Pictures/SphereSimulated.PNG)


from which the stereographic projection can be obtained:

![Simulated_30](/Pictures/StereoSimulated.PNG)

For comparison the stereographic projection from _DynamicS_ is given below:

![Simulated_30](/Pictures/SimulationSP.bmp)


### 2.2 Measured patterns
Comment the second PatterList images to show the reconstruction of the measured patterns instead of the simulated.


![Simulated_30](/Pictures/StereoMeasured.PNG)


## 3 Required data
The big advantage of our methode is that no indexing is required for the reconstruction. 
Therefore, some experimental parameters must be recorded beside the patterns: the distance between the phosphor screen and the specimen (relative to the diameter of the phosphor screen), the pattern center (relative to the diameter of the taken image) and the angles phi and psi. Latter is usually 70°.
Furthermore, the patterns must be recorded at the same specimen position (same crystallite). 
New scanning electron microscopes may provide an eucentric rotation.
Otherwise, the specimen position has to be searched after rotation around the specimen normal (phi) and the pattern center has to be calibrated again.

The size of the reconstructed pattern (effective size of phosphor screen) can be increased by the variation of the angles phi and psi and the shift of the pattern center (_e.g._ change in working distance).

**Increasing the effective size of the phosphor screen, allow better indexing and/or interpretation of complicated patterns due to pattern superpositions, pseudo-symmetry, etc.**


