# README
## POLARIS: Path of Least Action Recursive Survey

This repository contains the software implementation for our paper **POLARIS: path of least action analysis on energy landscapes** (`Seitz E. and Frank J., doi: 10.1021/acs.jcim.9b01108`). It contains tools to apply the discussed methods to energy landscapes.

### Installation
1.  Install latest version of Anaconda: https://www.anaconda.com/download/

2.  Use the `.yml` file found within the repository to install the environment for your operating system, via either:
         
	    conda env create -f env_linux.yml
	    
	    conda env create -f env_mac.yml
	 
    If any problems emerge, try to continue instead with the following commands:
    
    2a.  Create environment:

	    conda create --no-default-packages -n polaris python=3
	
    2b.  Open environment first before installing packages inside:

	    conda activate polaris
	
    2c.  Install packages in environment:

	    pip install --upgrade pip
  	    pip install PyQt5
  	    pip install matplotlib
 	    conda install -c anaconda psutil

3. If environment not already open, it must be sourced each time before program is run:

	   conda activate polaris
	 
4. Run program from POLARIS_GUI.py directory via:

	   python POLARIS_GUI.py
	
5. When done using program, always exit environment via:

	   conda deactivate

### Attribution:
Please cite `E. Seitz and J. Frank (2019); doi: 10.1021/acs.jcim.9b01108` if you find this method useful in your research. If needed, this GitHub repository can additionally be cited via [![DOI](https://zenodo.org/badge/157617482.svg)](https://zenodo.org/badge/latestdoi/157617482)

### License:
Copyright 2018-2020 Evan Seitz

For further details, please see the LICENSE file.
