
def inverseSignature(directory_interactions, directory_hits_recon): 

    '''
    U S A G E : 
    The function 'get_working_screens.py' allows us to obtain the number of screens for which the 'true target' was predicted using double gene deletion. 
    This ('inverseSignature.py') function allows us to determine if the genes that were predicted to be drug targets by 'get_working_screens.py' were specfific to their
    drug, or if those genes are important genes interacting with every drug (i.e. found in every CRISPR screen), in which case the results from 'get_working_screens.py'
    would be less significant. 
    Hence, this function builds a signature for every candidate gene, which consists of all the screens with which this gene is interacting. 
    Then, the number of true gene-screen interaction is recorded. A 'true interaction' is considered as such if the candidate gene is in the list of 
    CRISPR hit of the drug with which it is though to be interacting.  
    A example of usage is provided at the end. 

    I N P U T : 
    --> directory_interactions : the directory of the folder containing all the interaction files ('drugXXX_interaction.txt'). The genes in those files
        are the candidate target genes. Those files can be issues from the 'get_interactions.py' script. 
    --> directory_hits_recon : the directory of the folder containing all the hits (in BiGG id) of the CRISPR screens. Those files should be tab-delimited text files, 
        and can be issues from the 'screen_converter.py' script. It is those files that will help us determine if an interaction between a gene and a screen is considered true or not. 

    O U T P U T : 
    --> sucess_count : the number of gene-screen interaction considered to be true. 
    --> allcandidateDict : a dictionary containing the signature of every candidate gene. It is in the form : 
        {candidategene : {'drug with which it is though to be interacting':['other gene with which it is though to be interacting']}}
        {candidategene1 : {'drug1': ['gene56'], 'drug67': ['gene159']}, candidategene2 : {'drug171': ['gene76']}, candidategene3 :  {'drug187': ['gene86'], 'drug3': ['gene1880']}, etc.}

    N O T E : 
    Some gene are associated to a great number of screens, which makes them non-specific to the screen with which they are truly interacting, and add noise
    to the data. Therefore, only the gene interacting with less than half of the total number of eligible screens are considered to calculate this p-value. 
    This value can be changed, simply change the number in line 81 (0.5*x). 

    N O T E : 
    The function 'inverseSignature' holds two functions: 
    1 ) 'inverseSignature1' : which builds the dictionnary 'allcandidateDict' presented above
    2 ) 'inverseSignature2' : which uses the dictionnary to asses the number of 'true' gene-screen interactions. 
    '''

    import os 

    def is_float(string): 
        try :
            float(string)
            return True 
        except: 
            return False
    
    def inverseSignature1():
        # this function builds a dictionnary of the target genes, and the screen with which they interact 
        from converter5 import allToat
        allToatDict = allToat()
        allcandidates = []

        for fileName in os.listdir(directory_interactions): 
            filePath = os.path.join(directory_interactions, fileName)
            with open(filePath, 'r') as file : 
                for lines in file : 
                    if is_float(lines.split('\t')[0][0]) == True : 
                        candidate = allToatDict[lines.split('\t')[1]]
                        allcandidates.append(candidate)
        allcandidates = set(allcandidates)
        allcandidates = list(allcandidates)



        allcandidateDict = {}
        x = len([name for name in os.listdir(directory_interactions) if os.path.isfile(os.path.join(directory_interactions, name))])
        for i in range(len(allcandidates)):
            candidateGene = allcandidates[i]
            geneScreenDict = {}
            for fileName in os.listdir(directory_interactions): 
                screenName = fileName[:-19]
                filePath = os.path.join(directory_interactions, fileName)
                interGene = []
                with open(filePath, 'r') as file : 
                    for lines in file : 
                        if is_float(lines.split('\t')[0][0]) == True : 
                            if allToatDict[lines.split('\t')[1]] == candidateGene: 
                                interGene.append(lines.split('\t')[0])
                if len(interGene) != 0 : 
                    geneScreenDict[screenName] = interGene
            if len(geneScreenDict) < (0.5*x): # can change this 
                allcandidateDict[candidateGene] = geneScreenDict

        #allcandidateDict looks something like this 
        #KEYS : VALUES 
        #1119_AT1 : {'Dihydrosphingosine_8uM': ['8879_AT1'], 'FTY720_3uM': ['56848_AT1']}
        #1503_AT1 : {'Aphidicolin_A_0_4uM': ['1635_AT1'], 'Dorsomorphin_5uM': ['50484_AT1'], 'Wisent_FBS_Performance_098-150_lot_185705': ['1635_AT1']}
        #1537_AT1 : {'AMTT_10uM': ['51805_AT1', '51004_AT1', '27235_AT1', '10229_AT1'], 'Aphidicolin_A_0_4uM': ['51805_AT1'], 'FK-506_30uM': ['10654_AT1'], 'Hippuristanol_0_12uM': ['51805_AT1', '10229_AT1', '51004_AT1'], 'HMS-I2_10_uM': ['10229_AT1'], 'MK2206_4uM': ['10229_AT1'], 'Nutlin-3A_1_6_uM': ['51004_AT1'], 'Phenformin_20_uM': ['10654_AT1'], 'Pimelic_High': ['4598_AT1'], 'QNZ_0_01uM': ['10654_AT1'], 'Quercetin_20uM': ['51004_AT1'], 'Thimerosal_0_85uM': ['10229_AT1'], 'Vemurafenib_6_6uM': ['10229_AT1']}

        return allcandidateDict


    def inverseSignature2(): 

        allcandidateDict = inverseSignature1()
        candidateGene = list(allcandidateDict.keys())
        geneScreenDict = list(allcandidateDict.values())

        sucess_count = 0
        count = 0

        for i in range(len(candidateGene)): # iterating throught genes 
            associated_screen_name = list(geneScreenDict[i].keys()) 
            n = len(associated_screen_name)
            for j in range(len(associated_screen_name)): # iterating through gene-screen relation
                fileName = associated_screen_name[j] + '_hits_Recon.txt'
                filePath = os.path.join(directory_hits_recon, fileName)
                count += 1
                with open(filePath, 'r') as file : 
                    realHits = []
                    for lines in file : 
                        realHits.append(lines.split('\t')[0])
                if candidateGene[i] in realHits: 
                    sucess_count += 1


        print('sucess_count : ', sucess_count)
        print('count:', count)



        return sucess_count, allcandidateDict
    
    
    return inverseSignature2()

if __name__ == '__main__': 
    # example of usage 
    directory_interactions = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'
    directory_hits_recon = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    inverseSignature(directory_interactions, directory_hits_recon)
