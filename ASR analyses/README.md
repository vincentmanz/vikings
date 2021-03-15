# Exploratory script using Ancestral State Reconstruction


The first time follow the procedure


# Python environment
Install vistual environment 
	python3 -m pip install virtualenv

### create env
	cd ~
	python3 -m venv vikings
	source ./vikings/bin/activate
	cd PATH/TO/VIKING/PROJECT
	pip install --upgrade pip 
	pip install -r requirement_python.txt 


# R environment
 sudo apt install r-base-core
 Rscript requirement_R.txt

**If the python environment and the R module are installed run: **
	source ./vikings/bin/activate

Then run the script
	./vikings_ASR.py -t 13 -f ./csv/



To go out of the python environment (not necessary)
	deactivate
