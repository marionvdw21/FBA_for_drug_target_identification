def screen_converter(directory, directory_hits, directory_hits_recon , myModel): 

    '''
    U S A G E : 
    This function takes a folde with CRISPR screens and extracts the hits and converts them to BiGG ID. It works in two parts : 
    1 ) FindHITS : Within each CRISPR screens, the hits (those genes with an FDR < 0.05) are found. For each file in 'directory', a file 
        containing only the crispr hits will be outputed in 'directory_hits'
    2 ) ConvertHITS : For each file, the hits are converted from HGNC symbols to BiGG ids. 
        For each file in 'directory_hits', a file translated to BiGG ids will be outputed in 'directory_hits_recon'. 
    An example of usage is provided at the end. 

    I N P U T S : 
    --> 'directory' : the directory of the folder containing all the raw crispr screens. The files in this folder should be tab delimited txt files 
        with HGNC symbols at the first colomn, and the FDR values at the fourth column. An example is provided as '1_methyl_nicotinamide.txt'.
        The folder corresponding to this directory should be full of all the crispr screens of the drugs of interest. 
    --> 'directory_hits' : the directory in which the files of the CRISPR screen filtered for only the hits will be outputed. 
        The folder corresponding to this directory should be empty. 
    --> 'directory_hits_recon' : the directory in which the translated files (from HGNC to Recon) will be outputed.
        The folder corresponding to this directory should be empty. 

    O U T P U T S : 
    --> 'drugXXX_hits' : for each drug, a corresponding file will be outputed in the folder of 'directory_hits', with only the genes that are hits. 
    --> 'drugXXX_hits_Recon' : for each drug, a corresponding file of hits translated to BiGG id from HGNC will be outputed in the 
        folder of 'directory_hits_recon'. 
    --> eligible_screen : the number of eligible drug crispr screens, which means they contain crispr hits and that the crispr hits are gene in the model. 


    '''

    import HGNC_to_Recon2
    import os
    def findHITS():
        for fileName in os.listdir(directory): 
            newLines = []
            filePath = (os.path.join(directory, fileName))
            outFileName = fileName[:-4] + '_hits.txt'
            outFilePath = os.path.join(directory_hits, outFileName)

            with open(filePath, 'r') as file : 
                next(file)
                for lines in file : 
                    if float(lines.split('\t')[3]) < 0.05: 
                        newLines.append(lines)
            with open(outFilePath, 'w') as file : 
                for i in range(len(newLines)): 
                    file.write(newLines[i])
        return 
    findHITS()


    def convertHITS(): 
        import os 
        for fileName in os.listdir(directory_hits): 
            filePath = os.path.join(directory_hits, fileName)
            check = 0
            with open(filePath, 'r') as file : 
                if len(file.readlines()) > 2 : 
                    check = 1
                if check == 1: 
                    outFileName = fileName[:-4] + '_Recon.txt'
                    outFilePath = os.path.join(directory_hits_recon, outFileName)
                    HGNC_to_Recon2.HGNC_to_Recon_converter(filePath, outFilePath, myModel) 

        # check how many screns out of all the folder are eligible to double gene deletion 
        # (i.e. have CRISPR hits that are in our model)
        eligible_screen = 0
        for fileName in os.listdir(directory_hits_recon): 
            filePath = os.path.join(directory_hits_recon, fileName)
            with open(filePath, 'r') as file : 
                if len(file.readlines()) > 0 : 
                    eligible_screen += 1
        return eligible_screen
    eligible_screen = convertHITS()

    return eligible_screen


if __name__ == '__main__': 
    # Example of usage 
    import cobra.io as co
    directory = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS'
    directory_hits = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS'
    directory_hits_recon = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    screen_converter(directory, directory_hits, directory_hits_recon, myModel = co.read_sbml_model('Recon3.xml'))