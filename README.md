## Open the Black-box: Interpreting the Machine-learned Low-dimensional manifold

<img src="./Images/Highlight.png" width=300 />

### 1. Highlight

- This work explained the physical insights discovered by a data-driven dimensionality reduction approach (i.e., Active Subspace).
- We analytically derived a low-dimensional representation of the original system, thanks to the insights revealed by the Active Subspace approach.
- We achieved even faster instability risk calculations by eliminating costly data-driven subspace identification processes.

This work was firstly presented in [37th International Symposium on Combustion](http://www.combustionsymposia.org/2018/loadpage/program/program), and was later published in the journal:

Guo S., Silva C. F., Bauerheim M., Ghani A., Polifke W., [Evaluating the impact of uncertainty in ﬂame impulse response model on thermoacoustic instability prediction: A dimensionality reduction approach](https://www.sciencedirect.com/science/article/abs/pii/S1540748918304383). *Proceedings of the Combustion Institute*, 2019, 37(4), pp. 5299-5306.

### 2. Motivation

Accurately quantify the uncertainty in combustion instability prediction is crucial for the industry to design reliable gas turbine systems.

Previously, by employing a data-driven dimensionality reduction approach (Active Subspace), we discovered an active direction inside of the high-dimensional parameter space, such that along this direction, the output of interest (combustion instability index) varies the most. Due to the black-box nature of the data-driven approach, we are hidden from the physical interpretation of the identified active direction. Without this knowledge, we cannot be certain when the data-driven approach will succeed or fail.

### 3. Methodology

To explain the physical interpretation conveyed by the derived low-dimensional manifold, we worked with the domain experts and start with manipulating the physical equations that govern the combustion instability phenomenon:

- We linearized the governing equations and performed scale analysis to eliminate unimportant terms/quantities.

- We employed the phasor plot of the flame transfer function to visualize results. This visualization practice is not common in the community but helped deliver rich insights in the framework of the current study. 

<img src="./Images/Phasor.png" width=500 height=300 />

### 4. Results

- We analytically recovered the results yielded by the Active Subspace approach.

<img src="./Images/Results.png" width=500 />

- We are able to visualize this active direction, and we can analytically derive the coefficients of this 1D projection.

- We achieved 80 times faster instability risk calculations by eliminating costly data-driven subspace identification processes.

- The insights revealed by the data-driven approach and our analytical derivations allow us to quickly screen the sensitivity of individual input parameters.

<img src="./Images/Sensitivity.png" width=300 />


### 5. Folder structure

**1. Presentation**: the slides presented in [37th International Symposium on Combustion](http://www.combustionsymposia.org/2018/loadpage/program/program) conference.

**2. MatlabScripts**: MATLAB source code and data to reproduce the results. The code and data are organized in individual folders corresponding to different sections in the paper. 