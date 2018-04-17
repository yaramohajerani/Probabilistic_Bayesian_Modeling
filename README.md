# Probabilistic_Bayesian_Modeling
Bayesian Modeling Scripts in Python using Stan and Edward

Note compiled Stan files (`.pkl`) are ignored by git so if you
clone this repository you have to compile the models in
'compiled_models.dir' in your python directory.

Also data/ is set to be ignored in `.gitignore`. However, the code is
setup such that all the data is outside of the github folder right now.

In the same location where you see the folder `Stan_Bayesian_Modeling` make
another folder called `poc_data`, within which there should be two subdirectories
called `indata.dir` and `outdata.dir` for the input and output data, respectively.
Place the `pom_flux` folder from Team Drive into `indata.dir` before running the script.

Please make sure you have the latest version of PyStan. The earlier
versions gave me segmentation fault errors for some chains.
