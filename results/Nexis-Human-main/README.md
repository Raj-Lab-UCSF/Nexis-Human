# Nexis-Human

Nexis set up to fit to individual human data (gene or global)

1. Go to 'Running_Nexis_human.ipynb.' This is the notebook where you will prepare data, optimize Nexis parameters, and run the model on every patient.

   - The model itself is in 'Nexis_human_model.py' and necessary functions are in 'Nexis_human_functions.py.'

2. Upload the connectivity matrix found in 'Connectivity matrix original.csv.'

3. Select initial vector method. 'Binary' sets a region of your selection to 1 and every other regions to 0. Default is entorhinal cortex. 'Baseline' sets the initial vector to t=0 from the EBM dervied regional tau interpolations for the whole ADNI3 cohort.

   If you select 'binary,' the model will scale the vector using parameter Gamma so that it matches the empircal EBM data. If you select 'baseline,' Gamma will automatcially be set to 1.

   'Regional tau time series.csv' contains the cohort level interpolation of tau in every region across time, dervied from EBM.

4. Load DK regional volumes for volume correction.

5. Upload gene data and enter gene of interest. Make sure data is in DK order first (the code is set up to reorder it appropriately). It should be a matrix of nROI x number of genes.

6. Define inputs to Nexis. U_global is set to a matrix of zeros such that gene mediated diffusion/aggregation is not modeled. U_gene is set to the gene data matrix you just uploaded.

7. Upload patient data from 'Cross-sectional stage and regional tau.csv,' which contains the most likely stage assignment from EBM and regional tau for every ADNI3 subject. The stage and corresponding tau vector for each subject is what Nexis will fit to.

8. Initially run the model on only one subject to make sure it is working before looping through all subject.

   a. Optimize parameters. You can set it to optimize either Nexis-global or Nexis-gene. Test both.

   b. Run Nexis with fitted parameters. Again, you can set it to run either model (make sure its the same version as the parameters you just optimized).

   c. Check results. Calculate R. Total tau plot gives a sense of how well the model is fitting to the individual's data.

9. Loop through all patients. This loop is set up to perform step 8 for every subject for both Nexis-global and Nexis-gene. Output is a data frame with with R^2 and the optimized parameters for both versions of the model for every subject.
