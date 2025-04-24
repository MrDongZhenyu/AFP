# Analytic Fourier Ptychotomography (AFP) 

This repository contains the codes and example data for Analytic Fourier Ptychotomography (AFP), a technique used for volumetric refractive index imaging. It includes both the simulation and reconstruction codes. 

**Note:** Part of the codes are temporarily encrypted and will be fully released upon publication of the paper. The provided demos still work as expected.

Paper link: 

arXiv: https://arxiv.org/abs/2504.16247

Project page: https://mrdongzhenyu.github.io/AFP-Web/

Data source: https://osf.io/f7tqa/

Top-level folder structure:

```bash
├── Data                               # Directory for raw AFP experimental data
├── sampleFunctions                    # Sample functions used in the simulation
├── subfunctionAFP                     # Directory for AFP subfunctions used in the reconstruction process
├── Main_AFP_experiment.m              # Main Matlab script for experimental data reconstruction
├── Main_AFP_simulation.m              # Main Matlab script for AFP simulation
├── setSystemParameters.m              # Script to set system parameters, modify this for your customized experimental system
└── README.md                          # Project documentation (this file)
```

## Table of Contents
- [Introduction](#introduction)
- [Dependency](#dependency)
- [Data source](#data-source)
- [Usage](#usage)
- [Contact](#contact)
- [BiBTeX](#bibtex)
- [Copyright](#copyright)

## Introduction
AFP is an analytic technique used for refractive index tomography of weakly scattering samples. It retrieves the aberration function of an imaging system and extends the Fourier spectrum into the darkfield regime, both analytically, under the assumption of a finite sample thickness (FST) prior. By incorporating this prior, AFP transforms the inverse scattering problem into a linear problem, enabling more robust and efficient reconstruction of the refractive index distribution, without any paramter tuning or iterative optimizations. 

## Dependency
The codes are written in Matlab, and tested with versions 2021b-2024b. 

Packages required: Image Processing Toolbox & Symbolic Math Toolbox.

If GPU computing is enabled, the Parallel Computing Toolbox is also needed.

## Data source
Please download data and put them in the 'data' folder. GitHub 'data' folder currently gives a few example experimental data including beads, embryos, root and pathology. Some of the data are cropped because of GitHub data size limit, please download the entire data for best performance.

## Usage
### 1. AFP simulation
Run `Main_AFP_simulation.m`, choose `'sampleName'` as "bead" or "lymphn node", use `'approxMethod'` as "Born" or "Rytov", and set `'useAbeCorrection'` and `'useDarkfield'` to choose whether to do aberration correction and darkfield extention respectively. 

You can define your own simulation sample functions and put in the 'sampleFunctions' folder, change the sample thickness (`'samThick'`), the background refractive index (`'n_media'`), the aberration function (`'zernikeGT'`) and the illumination pattern (`'kIllu'`) as you want. The program also visulaizes the comparison with the ground truth after AFP reconstruction. 

### 2. AFP experiment
1. Define your system paremters in `setSystemParameters.m` first. By default, it includes parameters for a Laser system, an LED array, and a UV LED setup already. 

2. Run `Main_AFP_experiment.m`. In default, we have sample choices of "bead", "embryo", "root" and "pathology". Set `'useAbeCorrection'` and `'useDarkfield'` to choose whether to do aberration correction and darkfield extention respectively. `'useAbsorption'` is only set to `true` for the pathology sample. The illumination `'kIllu'` can be replaced by the designed/calibrated values for your own system.  

    Note: set  `'EnableROI'` to `ture` if you want to select a specific ROI to do the reconstruction with your own data.
    
### 3. Basic Parameters
1. `'pixelsizeZ'`: Can be set to be different from `'pixelSizeXY'`, i.e. an uneven grid size is allowed.
2. `'highresPad'`: To generate a high-resolution image, upsampling is typically requried due to the requirement of Nyquist sampling. `'highresPad'` tells the program the upsampling ratio.
3. `'approxMethod'`: We recommand to use "Rytov" approximation in the embryo reconstruction, and "Born" in the rest samples.
4. `'samThick'`: Sample thickness in pixels, depending on the sample used, the thicker the sample, the smaller region the FST prior can extend. 


### 4. Important options for AFP's sub-functions
We note that parameter tuning is unnecessary as AFP is an analytic method. However, there are several options that you can choose in performing the reconstruction. 

#### `recFieldKK.m` function
Two important optional arguments for `recFieldKK` are `'CTF'` and `'norm'`.
1. `'CTF'`: Incorporate (the absolute value of) the coherent transfer function for noise reduction. This will force the regions that are not covered by the CTF in reconstructed spectrum to be zero.
2. `'norm'`: Whether to use normalization in the reconstruction. We assume that the amplitude of zero-frequency in the sample's spectrum does not change with repect to the tilted illumination. Thus, when `'norm'` is set to `true`, we force the zero-frequency component in reconstructed spectrums to have the same amplitude (This is essentially correcting for illumination intensity variation).

#### `findAbeFromOverlap3D` function
Three important optional arguments for `findAbeFromOverlap3D` are `'weighted'`, `'similarity tolerence'` and `'vis_samplingmap'`.
1. `'weighted'`: When it is `true`, the program focuses more on places with larger weights and tends to ignore places with smaller weights, we have four consideration detailed in this function. When it is `false`, the function does not emphasize on any particular place and treats the signals equally. It is set to `true` by default
2. `'similarity tolerence'`: This parameter is used to calculate the overlap region using the FST prior. This tolerence value determines the half-width of a function (such as gaussian or sinc) when it drops to the `'tolerence'` percentage of its peak value. We use a default value of $\frac{\sqrt{2}}{2}$.
3. `'vis_samplingmap'`: Whether to visualize the aberration sampling map. 

#### `recFieldFromKnown3D` function
Six important optional arguments for `recFieldFromKnown3D` are `'drift'`, `'reg'`, `'threshold'`, `'intensity correction'`, `'similarity tolerence'` and `'use data intensity'`
1. `'drift'`: Whether to consider calibration error in the illumination angle. When set to `true`, this function provides more robust results. However, this results in assuming a smaller known spectrum, which would in turn require a larger overlap ratio.
2. `'reg'`: Regularization term to consider the effect of illumination angle misalignment, a larger number is needed for larger misalignments. It is recommended to choose a number in between `[1, 10]` for experiments. 
3. `'threshold'`: When the overlap ratio of the data is large, we can use this argument to improve the SNR of the reconstruction. The program reconstructs some regions of the spectrum multiple times using different measurements and then averages over these independent reconstructions. It is recommended to choose a threshold number in between `[0, 0.4]`. A larger threshold number leads to a longer reconstruction time.
4. `'intensity correction'`: When set to `true`, the function automatically compensates for illumination intensity differences.
5. `'similarity tolerence'`: As the darkfield measurements have disjoint support in 3D, this number defines the tolerence when the spectrum at a slightly different $k'_z = k_z + \delta k_z$
is treated to be approximately the same as the spectrum at $k_z$. If given as input, this is used to generate a sample-thickness-dependent variable (for the FST prior) computed within this function. In default, it is set to 0.3.
6. `'use data intensity'`: When set to `true`, the function automatically use the amplitude of the darkfield measurement to maintain the pixel-wise energy for the analytically solved darkfield spectrum.

### 5. How to run the code for your own data?
1. Prepare the data in the same format as the example data (.mat file)
    
    Data File (.mat) variables:

    * `'imStack'`: 3D data stack [NumOfPixels along x, NumOfPixels along y, numImages], format: "uint8", "uint16", "double" or "single"
    * `'ExposureTime'`: exposure time for each image (normalization purpose)
    * `'NA_Illu'`: illumination angle list (numImages-by-2).
        * Note 1: illumination NA = sqrt(NA_Illu(idx,1)^2 + NA_Illu(idx,1)^2)
        * Note 2: `'NA_Illu'` is dimensionless and the values belong to `[0, 1]`
2. Set the sample name and data path in Step 1 of the code
3. Set the parameters in `setSystemParameters.m` (design you own options)
4. Add calculation of `'kIllu'` in Step 3 in your sample name option

    Code example:
    ``` matlab
    kIllu = NA_Illu / NA * maxCTF;   % unit: 1/um
    kIllu(1:nNAmatching,:) = kIllu(1:nNAmatching,:)*0.98;   % for NA matching
    ```
5. Run the code `Main_AFP_experiment.m`

    * Note: For datasets with large field-of-view (FOV), it is recommended to divide the FOV into smaller patches of 256 x 256 pixels (as an example). After processing each patch, you can then stitch the AFP reconstructions back together to form a larger FOV. This approach helps mitigate field-dependent aberrations, as demonstrated in our Root and Pathology results.

## Contact
For any questions or comments about this code, please contact [zdong@caltech.edu](mailto:zdong@caltech.edu) or [rcao@alumni.caltech.edu](mailto:rcao@alumni.caltech.edu)

## BiBTex
```
@misc{dong2025analyticfourierptychotomographyvolumetric,
      title={Analytic Fourier ptychotomography for volumetric refractive index imaging}, 
      author={Zhenyu Dong and Haowen Zhou and Ruizhi Cao and Oumeng Zhang and Shi Zhao and Panlang Lyu and Reinaldo Alcalde and Changhuei Yang},
      year={2025},
      eprint={2504.16247},
      archivePrefix={arXiv},
      primaryClass={physics.optics},
      url={https://arxiv.org/abs/2504.16247}, 
}
```

## Copyright
Copyright (C) 2025 @ Biophotonics Laboratory, Caltech.

The code is licensed under the GPL-v3.0 license.

