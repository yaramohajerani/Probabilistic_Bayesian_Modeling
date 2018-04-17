# Probabilistic_Bayesian_Modeling
Bayesian Modeling Scripts in Python using Stan, Edward, and PyMC3

Note compiled Stan files (`.pkl`) are ignored by git so if you
clone this repository you have to compile the models in
'compiled_models.dir' in your python directory.

Also data/ is set to be ignored in `.gitignore`. However, the code is
setup such that all the data is outside of the github folder as of now.

In the same location where you see the folder `Probabilistic_Bayesian_Modeling` make
another folder called `poc_data`, within which there should be two subdirectories
called `indata.dir` and `outdata.dir` for the input and output data, respectively.
Place the `pom_flux` folder from Team Drive into `indata.dir` before running the script.

Please make sure you have the latest version of PyStan. The earlier
versions gave me segmentation fault errors for some chains.

If Edward gives you problems with the `computation` module in pandas, update `dask`.

If PyMC3 gives you problems with sampling the posterior, update `Numpy` to the latest version.
