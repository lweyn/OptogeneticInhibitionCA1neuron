# OptogeneticInhibtion of a CA1 pyramidal neuron
This folder contains the code associated with the simulation study that is reported in 

Weyn, L., Tarnaud, T., Schoeters, R., De Becker, X., Joseph, W., Raedt, R., Tanghe, E. (2025) *Computational analysis of optogenetic inhibition of CA1 neurons using a data-efficient and interpretable potassium and chloride conducting opsin model*

The code was developed in the anaconda environment described in the environmnet.yml file using NEURON version 8.2.0 [1].  

## Model files
- Mods_Gentilette should contain kcc2.mod and nakpump.mod files described in Gentiletti et al. (2022), available at modeldb.science/267499 [2,3]
- Mods_Tomko and CA1PYR_TK21.hoc contain the model files for the CA1 pyramidal neuron described in Tomko et al. (2021) [4] 
- Mods_Opto contains the .mod files of the opsin models and other tools to enable optogenetic intervention in the model.

## Notebooks
- OpsinModelFit.ipynb contains an example of how the opsin model can be fit to an experimental dataset.
- runSim.ipynb contains an example of how the simulations can be run. 

## References
[1] N. T. Carnevale and M. L. Hines, The NEURON Book. Cambridge University Press, 1 (2006).

[2] Gentiletti D, de Curtis M, Gnatkovsky V, Suffczynski P.  Focal seizures are organized by feedback between neural activity and ion concentration changes eLife. 11 (2022).

[3] McDougal RA, Morse TM, Carnevale T, Marenco L, Wang R, Migliore M, Miller PL, Shepherd GM, Hines ML. Twenty years of ModelDB and beyond: building essential modeling tools for the future of neuroscience. J Comput Neurosci. 42(1):1-10.(2017)

[4] Tomko, M., Benuskova, L. & Jedlicka, P. A new reduced-morphology model for CA1 pyramidal cells and its validation and comparison with other models using HippoUnit. Sci Rep 11, 7615 (2021).


