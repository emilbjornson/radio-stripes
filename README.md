MMSE-Optimal Sequential Processing for Cell-Free Massive MIMO With Radio Stripes
==================

This is a code package is related to the following scientific article:

Zakir Hussain Shaik, Emil Björnson, Erik G. Larsson, “[MMSE-Optimal Sequential Processing for Cell-Free Massive MIMO With Radio Stripes](https://arxiv.org/pdf/2012.13928.pdf),” IEEE Transactions on Communications, To appear.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. We encourage you to also perform reproducible research!


## Abstract of Article

Cell-free massive multiple-input-multiple-output (mMIMO) is an emerging technology for beyond 5G with its promising features such as higher spectral efficiency and superior spatial diversity as compared to conventional multiple-input-multiple-output (MIMO) technology. The main working principle of cell-free mMIMO is that many distributed access points (APs) cooperate simultaneously to serve all the users within the network without creating cell boundaries. This paper considers the uplink of a cell-free mMIMO system utilizing the radio stripe network architecture with a sequential fronthaul between the APs. A novel uplink sequential processing algorithm is developed, which is proved to be optimal in both the maximum spectral efficiency (SE) and the minimum mean square error (MSE) sense. A detailed quantitative analysis of the fronthaul requirement or signaling of the proposed algorithm and its comparison with competing sub-optimal algorithms is provided. Key conclusions and implications are summarized in the form of corollaries. Based on the analytical and numerical simulation results, we conclude that the proposed scheme can significantly reduce the fronthaul signaling, without compromising the communication performance.

## Content of Code Package

The article contains 6 simulation figures, called Figure 4a, 4b, 5a, 5b, 6, and 7. Figure X is plotted using the Matlab script generateplotsFigureX.m and, in most cases, one first has to generate the simulation data using the corresponding script called ScriptGenerateDataforFigureX.m
The package also contains several Matlab functions that are called using the aforementioned scripts.

See each file for further documentation.

## Acknowledgements

E.~Björnson was supported by the Grant 2019-05068 from the Swedish Research Council.

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
