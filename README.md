# FBA_for_drug_target_identification

## **Project description**

## **Prerequisite and Installation** 

Make sure you have the correct [requirements](https://github.com/marionvdw21/FBA_for_drug_target_identification/blob/main/requirements.txt) and download the 'make_model' and ''

## **Usage**

Using metabolic modeling for drug target prediction requires : 
1) Adapting and validating a model from a generic model
2) Using this model to attempt drug target prediction.

Hence, the scripts are divided into two folders, which can be used separately.

### Make the model 
This folder is used to make and validate a metabolic model, starting from a generic model such as Recon models. 
The 'mainModel.py' script is the main script which calls the other functions in order to go from a generic model to a model adapted to a specific cell line. 
It contains many step, some of them being longer to execute than others, however each step can be carried out separately by calling the associated function. 

The user-provided inputs are the following: 

- myModel : object of class cobra.core.Model. It is the basic model which will be modified, for example Recon3. 
  It can be generated by downloading an SBML file in the directory (from http://bigg.ucsd.edu/models from example), and reading it using cobra : 
```
import cobra
myModel = cobra.io.read_sblm_model('Recon3.xml')
```

- 'expression.txt' : a text file containing the gene names and their associated expression in the cell line of interest cells (log2 read counts).

- 'minimal_media.txt' : a text file containing the reaction ids of the import and export reactions that are allowed in order to create the desired media.

- 'essentiality_Recon.txt' : a text file containing the experimental data of the cell line's gene ranked by essentiality, and their assciated rank score.


For each of the user-provided inputs, an example is given in the 'make_model' folder. 


### Use the model for drug target prediction

Once the model is made similar to the cell line with which the CRISPR screens have been done, it can be used to try to predict drug targets. 
For each drug, double gene deletions of every CRISPR hit of that drug VS all the other genes of the model will be performed in every pairs possible. 
The genes that are suspected to be interacting with CRISPR hits of that drug are considered candidate targets for that drug. An assesment of wether or not the true target was in the list of candidate target for that drug will be made, an a p-value will be associated to it. 

Many functions are involved in the above process. Unfortunately, there is no main function calling all of the others yet, partly because some of these functions take a very long time to execute, so it is better to divide them up. However, there is a pdf explaining how all of these functions call each other and how they should be used in the 'target_prediction' folder, under the name 'function_explanation'. 


## **License, Contact, Information** 

### Acknowledgments 
The cell line expression data and essentiality data ('essentiality.txt' and 'expression.txt') are provided by : Bertomeu, T., Coulombe-Huntington, J., Chatr-aryamontri, A., Bourdages, K. G., Coyaud, E., Raught, B., Xia, Y., &amp; Tyers, M. (2018). A high-resolution genome-wide CRISPR/Cas9 viability screen reveals structural features and contextual diversity of the Human Cell-Essential Proteome. Molecular and Cellular Biology, 38(1). https://doi.org/10.1128/mcb.00302-17 

We gratefully acknowledge support from McGill's Summer Undergraduate Research in Engineering programm. 

### Contact Information

email : marion.vandewynckele-bossut@mail.mcgill.ca


