# MiniBatch-QNP-NonlinearReco
## MATLAB Implementation

[Tao Hong](https://hongtao-argmin.github.io), [Thanh-an Pham](https://scholar.google.com/citations?user=gBO51K4AAAAJ&hl=fr), [Irad Yavneh](https://irad.cs.technion.ac.il), and [Michael Unser](https://scholar.google.com/citations?user=nKVDcQoAAAAJ&hl=en), ``[A Mini-Batch Quasi-Newton Proximal Method for Constrained Total-Variation Nonlinear Image Reconstruction](https://arxiv.org/abs/2307.02043)'', to appear in SIAM Journal on Imaging Sciences, 2026.


Our algorithms are evaluated on 3D optical diffraction tomography, and the experimental setup heavily relies on the [Lippmann-Schwinger Project](https://github.com/ThanhAnPham/Lippmann-Schwinger). If you use our code or reproduce our experimental settings, please also cite the relevant publications associated with the [Lippmann-Schwinger Project](https://github.com/ThanhAnPham/Lippmann-Schwinger).

In addition, our implementation is built upon the [GlobalBioIm library](https://biomedical-imaging-group.github.io/GlobalBioIm/). If you encounter any difficulties in setting up the experiments, we can share the specific version we used for reproducibility purposes. All intellectual credit remains with the original authors.

To use the full package, please follow these steps:

1.	Download the GlobalBioIm library and add it to your MATLAB path.
2.	Download this package and add it to your MATLAB path.
3.	Navigate to /simulation/main_simulate3D.m to generate simulated data.
Generating the Green’s function can be time-consuming. For a quick test, we provide pre-generated data for a 64^3 volume here:
[Download](https://drive.google.com/drive/folders/1M2cwVWwaagiqX_QiKiyW07uaq3Kpy5M9?usp=sharing)
4.	Run ODT_Reconstruction_Born.m (demo for Born-model-based reconstruction)
    or
    run ODT_Reconstruction.m (demo for Lippmann–Schwinger-model-based reconstruction).

Feel free to shoot me an email (tao.hong@austin.utexas.edu) if you would like to learn more about our work.

