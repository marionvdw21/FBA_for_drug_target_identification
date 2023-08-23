# FBA_for_drug_target_identification

## **Project description**

## **Prerequisite and Installation** 

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

## **License, Contact, Information** 

Chemogenomics data sets from : Jasmin Coulombe-Huntington,Thierry Bertomeu, Caroline Huard, Corinne St-Denis, Andrew Chatr-Aryamontri, Daniel St-Cyr, Mike Tyers. In preparation.


