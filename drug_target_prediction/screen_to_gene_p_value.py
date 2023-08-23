def screen_to_gene_p_value(target_table_Recon, directory_interactions, Ho_tot, candidate_gene_file = 'allSingleGeneDel.txt'): 

    '''
    U S A G E : 
    This function takes the results from 'getWorkingScreens.py' and associates a p-value to it, to verify if the screens 
    that worked are not a result of chance. 

    I N P U T : 
    --> target_table_Recon : as an input of the 'getWorkingScreen.py' script. See that script for more description on that input. 
    --> directory_interactions :  as an input of the 'getWorkingScreen.py' script. See that script for more description on that input. 
    --> Ho_tot : the 'n' variable for the p-value calculation. The larger the 'n', the more precise the p-value. 
    --> candidate_gene_file : a text file containing the pool of genes from which the random picker will pick during the p-value calculation. 
        The default is every gene in the model, but a more accurate pool of gene would be those more highly expressed gene, as they have more 
        chance of being more essential and being drug targets. 

    O U T P U T S : 
    --> a p-value : probability of getting the results of 'getWorkingScreens.py' by chance.
    '''

    import random
    import os 

    def is_float(string): 
        try : 
            float(string)
            return True
        except: 
            return False

    import getWorkingScreens
    ListScreenNames = getWorkingScreens.getWorkingScreens(target_table_Recon, directory_interactions)[0]
    number_of_sucess = getWorkingScreens.getWorkingScreens(target_table_Recon, directory_interactions)[1]
    screenDict = getWorkingScreens.getWorkingScreens(target_table_Recon, directory_interactions)[2]




    Ho_true = 0
    # making a list of all the 'true' targets
    targetContent = list(screenDict.values())
    

    for i in range(Ho_tot): 
        print(i)
        success = 0

        j = -1
        for screenName in ListScreenNames: 
            j += 1
            interGenes = []
            singleSucess = 0
            filePath = os.path.join(directory_interactions, screenName)
            with open(filePath, 'r') as file : # getting the number and name of interactions genes 
                for lines in file : 
                    if is_float(lines.split('\t')[0][0]) == True : 
                        interGenes.append(lines.split('\t')[1])
            interGenes = set(interGenes)
            n = len(interGenes)


            with open(candidate_gene_file, 'r') as file : # can change file
                randint = []
                content = []
                for lines in file : 
                    content.append(lines.split('\t')[0])

            for i in range(n): 
                picker = random.randint(2, 1100) # change depending on file --> picker = random.randint(2, 1849)
                while picker in randint : # to avoid duplicates 
                    picker = random.randint(2, 1100) # change depending on file --> picker = random.randint(2, 1849)
                randint.append(picker)
                randomGene = content[picker]
                #randomGene = allToatDict[randomGene] # change depending on file 
                
                if randomGene in targetContent[j]: 
                    singleSucess += 1

            if singleSucess > 0:
                success += 1 
        

        if success >= number_of_sucess : # this is the true alternative
            Ho_true  += 1
    
    print('Ho_true: ', Ho_true)
    p_value = Ho_true/Ho_tot

    print(p_value)
    return p_value 


if __name__ == '__main__': 
    # example of usage
    target_table_Recon = 'target_table_Recon.txt'
    directory_interactions = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'
    Ho_tot = 1
    screen_to_gene_p_value(target_table_Recon, directory_interactions, Ho_tot, candidate_gene_file = 'allSingleGeneDel.txt')

    