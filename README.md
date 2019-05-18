# POLARIS
algorithm for path of least action analysis on energy landscapes

# INSTALLATION
first, install latest version of Anaconda: https://www.anaconda.com/download/

create environment:

	conda create --no-default-packages -n POLARIS python=3.7
	
open environment:

	source activate POLARIS
	
install packages in environment:

	pip install --upgrade pip
  	pip install PyQt5
  	pip install matplotlib
 	conda install -c anaconda psutil
	 
run program from polaris_FE.py (frontend) directory via:

	python polaris_FE.py
	
when done using program, exit environment via:

	source deactivate
