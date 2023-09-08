def getHits(directory, directory2 = 'none'): 

    '''
    U S A G E : 
    This function takes a folder full of CRISPR screens for various drugs, and outputs a text file containing
    all the CRISPR hits (without the redunduncies). 
    An example of usage is provided at the end. 

    I N P U T S : 
    --> directory : the directory of the folder containing all the CRISPR screens. The files in the folder should have the genes in BiGG id, 
        and only the genes present in the model should be in the CRISPR screens. An example of such a file is provided as 'screen_name_hits_Recon.txt'. Only the first column, 
        containing the gene's BiGG ID, matters.  
        Such files can be outputed from the 'screen_converter.py' script. 
    --> directory2 (optional) : if some double gene deletions have already been made on some of the CRISPR hits, the function avoids doing them again. 
        The directory2 is the directory in which those files of already-done double gene deletion would be stored. 

    O U T P U T S : 
    --> 'list_of_hits.txt' : a tab-delimited text file with each line containing the name of a CRISPR hit gene that should undergo double-gene-deletion. 
        An example of such a file is provided under 'list_of_hits_example.txt'. 

    '''

    import os 
    # creating a set of all the CRISPR hits across all the screens in 'directory'
    allHITS = []
    for fileName in os.listdir(directory): 
        filepath = os.path.join(directory, fileName)
        with open(filepath, 'r') as file : 
            for lines in file: 
                allHITS.append(lines.split('\t')[0])
    allHIT = set(allHITS)
    allHIT = list(allHIT)


    # the following block of code avoids doing double gene deletion on genes that are already done 
    if directory2 != 'none': 
        doneGenes = []
        for fileName in os.listdir(directory2): 
            # the following line should be changed depending on the name of given to the double gene deletion files in your directory
            if fileName[:14] == 'HITSdoubleKO6_' and fileName[-15:] == '_AT1_VS_all.txt': 
                genelenght =len(fileName) - 28
                geneName = fileName[14:17+genelenght]
                doneGenes.append(geneName)
        print(len(doneGenes))
        toDoGenes = []
        for i in range(len(allHIT)) : 
            if allHIT[i] not in doneGenes: 
                string = ''
                # the 'something' could be anything, as long as there is no '\n' after the gene's BiggID
                string += str(allHIT[i]) + '\t' + 'something'
                toDoGenes.append(string)
    
    # writing the list of genes that needs to be done (double gene deletion) in a file
    with open('list_of_hits.txt', 'w') as file : 
        for i in range(len(toDoGenes)): 
            file.write(toDoGenes[i])
            file.write('\n')
       
    return 

if __name__ == '__main__': 
    # example of usage
    directory = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    directory2 = r'C:\Users\mario\OneDrive\Documents\SURE'
    getHits(directory, directory2)
