def getInteractions(singleGeneDeletionFile, directory_hits_recon, directory_interactions, lethal_cutoff, rescue_cutoff): 
    
    '''
    U S A G E : 
    The use of this file is to get every gene with which the crispr hits of each drug of interest is interaction. An interaction is considered as such if 
    the double knockout produces a significant change in the expected growth, which would be the product of the duble knockouts (rescue or lethal interaction). 
    For each drug in the 'directory_hits_recon' folder, a file will be created (and placed in the 'directory_interaction' folder) with the list of crispr hits, if any, that interacted
    with the some other genes in the model. 
    An example of usage is provided at the end. 

    I N P U T S : 
    --> singleGeneDeletionFile : the name of the file which contains the single gene deletion data for every gene in the model of interest. 
        It could be outputed from the 'single_gene_deletion.py' script. It should be a tab-delimited text file with the BiGG IDs in first column 
        in the format '['498_AT1', '498_AT2', '498_AT3']', and the second column containg the flux of the biomass maintenance reaction (which is the default solution when the model is optimzed).
    --> directory_hits_recon : the directory that leads to the folder which contains the txt file of each drug with its associated CRISPR hits, in BiGG id. 
        It could be outputed from the 'screen_converter.py' script. 
    --> directory_interactions: the directory which leads to the folder in which all the outputed interactions files will be placed. This folder should be empty. 
    --> lethal_cutoff : the cutoff at which an interaction between two genes is considered lethal. In that this context, 'lethal' is not used in the normal sense of 
        the term (in which deletion of the two genes would kill the cell), but rather means that the deletion of the two genes slows down the growth of the cell, 
        which suggest an interaction between the two genes. The number 'lethal_cutoff' represent the ratio of the double-gene-deletion growth to the single-gene-deletion growth. 
    --> rescue-cutoff : the cutoff at which an interaction between two genes is considered rescue. In that this context, 'rescue' is not used in the normal sense of 
        the term (in which deletion of the two genes saves the cell from death), but rather means that the deletion of the two genes increases the growth of the cell, 
        which suggest an interaction between the two genes. The number 'rescue-cutoff' represent the ratio of the double-gene-deletion growth to the single-gene-deletion growth. 
    
    O U T P U T S : 
    --> 'drugXXX_interactions.txt'
        For each drug, a text file with the predicted interaction between the crispr hits of that drug and other genes in the model will be created, and placed in the folder associated to
        'directory_interaction'. 

    N O T E S : 
    --> KO2 is the solution (flux through biomass maintenance reaction) of the model that underwent double gene knockout 
        These value come from the '{geneName}_VS_all.txt' file outputed from 'double_gene_deletion.py'. 
    --> KO1_1 is the solution (flux through biomass maintenance reaction) of the model that underwent single gene knockout, the single gene being the crispr hit. 
        These value come from the singleGeneDeletionFile. 
    --> KO1_2 is the solution (flux through biomass maintenance reaction) of the model that underwent single gene knockout, the single gene being any one of the gene in the model.
        These value come from the '{geneName}_VS_all.txt' file outputed from 'double_gene_deletion.py'. 
    --> no_sol_del is the solution (flux through biomass maintenance reaction) of the model that did not undergo any gene knockout. 

    --> The text files containing the results of the double knockouts for each crispr hit should be in the same directory as this function, 
        or else precise the specific directory in line 69.
        These text files should be named '{geneName}_VS_all.txt', and can be ouputed from the 'double_gene_deletion.py' script. 
    '''
    
    import os 

    # create dictionary to get KO1_2
    singleGeneDict = {}
    with open(singleGeneDeletionFile, 'r') as file : 
        for lines in file : 
            singleGeneDict[lines.split('\t')[0]] = lines.split('\t')[1]

    # Loop 1 : screens 
    countFile = 0
    for fileName in os.listdir(directory_hits_recon): 
        countFile += 1
        print(countFile)
        filePath = os.path.join(directory_hits_recon, fileName)
        newLines = []
        geneList = []
        with open(filePath, 'r') as file : 
            for lines in file : 
                geneList.append(lines.split('\t')[0])
        
        # Loop 2 : genes in screens 
        for geneName in geneList :
            geneFile = '%(geneName)s_VS_all.txt' %{'geneName' : geneName}
            header = ''
            header += 'Gene_name_1' + '\t' + 'Gene_name_2' + '\t' + 'biomass_flux' + '\t' + 'KO2:KO1_1' + '\t' + 'KO2:KO1_2' + '\t' 
            with open(geneFile, 'r') as file : 
                count = -1
                for lines in file : 
                    count += 1
                    myList_L = []
                    myList_R = []
                    if count == 0 : 
                        sol_no_del = float(lines.split('\t')[1])
                    if count == 2 : 
                        KO1_1 = float(lines.split('\t')[1])
                        newLines.append(header)
                    if count > 4 : 
                        KO2 = float(lines.split('\t')[1])
                        KO1_2 = float(singleGeneDict[lines.split('\t')[0]])

                        # ~ ~ ~ ~ ~ l e t h a l ~ ~ ~ i n t e r a c t i o n s ~ ~ ~ ~ 
                        if abs(KO1_1) > 0.01 : # protection against hit genes that interact with everyone 
                            if abs(KO1_2) == 0.0: # protection against 0 division 
                                KO1_2 = 0.00001
                            minimum = min(KO1_1, KO1_2)
                            # the 'if' statement below could be done another way, such as comparing double knowckout growth to the expected grwoth predicted by both the single knockouts (which is the more intuitive method) : 
                            # if KO2/((KO1_1 * KO1_2)/sol_no_del) < lethal_cutoff : 
                            # however the line above yielded to many lethal interaction so we switched for the alternative below. 
                            if KO2/minimum < lethal_cutoff :
                                myList_L.append(geneName) # gene name no1
                                myList_L.append(lines.split('\t')[0]) # gene name no2
                                myList_L.append(lines.split('\t')[1]) # flux thru biomass rxn
                                myList_L.append(str(KO2/KO1_1)) # KO2 / KO1_1
                                myList_L.append(str(KO2/KO1_2)) # KO2 / KO1_2
                                myList_L.append('lethal')

                        # ~ ~ ~ ~ ~ r e s c u e ~ ~ ~ i n t e r a c t i o n s ~ ~ ~ ~ 
                        if KO2 > 0.1: 
                            # the following lines aims to correct for the biomass production that end up being extremely small numbers that are not 
                            # not that different from each other (such as 5.536e-20 and -1.798e-23), but are taken as such if not corrected. 
                            if KO2 < 0.001 : 
                                KO2 = 0.001
                            if KO1_1 < 0.001 : 
                                KO2 = 0.001
                            if KO1_2 < 0.001 : 
                                KO1_2 = 0.001
                            # again, the 'if' statement below could be done in the more intuitive way of comparing double KO growth with the product of the two single KO growth, such as : 
                            # if KO2/((KO1_1*KO1_2)/sol_no_del) > rescue_cutoff:
                            # however the line above yielded to many rescue interaction so we switched for the alternative below. 
                            if KO2/KO1_1 > rescue_cutoff or KO2/KO1_2 > rescue_cutoff:
                                myList_R.append(geneName) #gene name no1
                                myList_R.append(lines.split('\t')[0]) #gene name no2
                                myList_R.append(lines.split('\t')[1]) #flux through biomass rxn
                                myList_R.append(str(KO2/KO1_1)) # KO2 / KO1_1
                                myList_R.append(str(KO2/KO1_2)) # KO2 / KO1_2
                                myList_R.append('rescue')

                    if len(myList_L) != 0 : 
                        newLines.append(myList_L)
                    if len(myList_R) != 0 : 
                        newLines.append(myList_R)

        outFileName = fileName[:-15] + '_interactions.txt'
        outFilePath = os.path.join(directory_interactions, outFileName)
        with open(outFilePath, 'w') as file : 
            for i in range(len(newLines)): 
                for j in range(len(newLines[i])): 
                    file.write(str(newLines[i][j]))
                    file.write('\t')
                file.write('\n')
    return 

if __name__ == '__main__': 
    # example of usage : 
    singleGeneDeletionFile = 'allSingleGeneDel.txt'
    directory_hits_Recon = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    directory_interactions = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'
    getInteractions(singleGeneDeletionFile, directory_hits_Recon, directory_interactions, 0.97, 1.002)