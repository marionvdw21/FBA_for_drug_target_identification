def make_and_validate_model():

    import cobra.io as co

    # Note : each step can be performed individually using the associated scripts and function. 

    # THE USER-PROVIDED VARIABLES REQUIRED : expression file and generic model from which to start 
    ExpressionFile = 'expression.txt'
    myModel = co.read_sbml_model('Recon3.xml')
    mediaFile = 'minimal_media2.txt'
    experimental_essentiality = 'essentiality_Recon.txt' 
    # if the experimental essentiality file has gene names in HGNC, it can be passed through the 'HGNC_to_Recon_converter' function (in 'HGNC_to_Recon.py' file) to be converted
    


    # --------------------------------------------------------------------
    # ---------- STEP 1: convert expression file to Recon ---------------
    # -----------------------------------------------------------------------
    import HGNC_to_Recon
    HGNC_to_Recon.HGNC_to_Recon_converter(ExpressionFile, myModel)
    import orderExpression
    orderExpression.orderExpression(expressionFile=ExpressionFile) 
    # !!!  need to check the resulting file, sometimes 'expression_Recon_ordered.txt' and 'expression_Recon_deordered.txt' need modifications. 



    # --------------------------------------------------------------------------------
    # --------- STEP 2 : use this converted file to create a new adapted model-------
    #---------------------------------------------------------------------------------
    import convertMyModel
    adapted_model = convertMyModel.convertMyModel('expression_Recon_ordered2.txt', myModel)[0]
    # returns the model adapted according to the expression file. 



    # ---------------------------------------------------------------------------------
    # ---------------------------- STEP 3 : convert the media ------------------------
    #---------------------------------------------------------------------------------
    # the media is made closer to the experimental media by regulating the import and exchange reactions 
    # that are allowed by the model. 
    import convertMyMedia
    adapted_model2 = convertMyMedia.convertMyMedia(adapted_model, mediaFile)



    # --------------------------------------------------------------------------------------------------
    # ----- STEP 4 : perform single gene deletion on every gene on the model ------------------------
    #-----------------------------------------------------------------------------------------------
    import single_gene_deletion
    single_gene_deletion.singleGeneDeletion(adapted_model2)
    # this step can take a long time to run, so it is recommended to use an HPC system to execute it 
    # !!! import adapted_model2 to directory of system used



    # ---------------------------------------------------------------------------------------------------
    # ------------ STEP 5 : gene essentiality prediction to validate the model -----------------------------
    # --------------------------------------------------------------------------------------------------
    # the genes are ordered in terms of essentiality (using the 'allSingleGeneDel.txt' file), and outputed 
    # in the 'EssPred.txt' file. 
    import Essentiality_prediction
    Essentiality_prediction.PredictEss('allSingleGeneDel.txt', 'EssPred.txt')
    # the 'EssPred.txt' file is then used with the experimental predicted essentiality to asses the accuracy of the model
    import compare_essentiality_prediction
    fileList = ['EssPred2.txt']
    compare_essentiality_prediction.fisherTest(fileList, experimental_essentiality, 200, -4, 'two-sided')



    #--------------------------------------------------------------------------------------------------------------
    # ------------ STEP 6 : optional step to create a file to visualize the reactions with non-zero fluxes ---------
    # ---------------------------------------------------------------------------------------------------------------  
    # the platform VMH (virtual metabolic human) can be used to create layovers on top of the Recon3 map 
    # to visulized which pathways are used once our modified model is optimized. The VMH_layover script 
    # outputs a file that can be used at the VHM website (https://www.vmh.life/#reconmap). 
    # Note : an account at VMH needs to be created in order to use this feature 
    import VHM_layover
    VHM_layover.makeLayover(adapted_model2)



    return


if __name__ == '__main__': 
    make_and_validate_model()
