def getWorkingScreens(myTargets_Recon, directory_interactions):

    '''
    U S A G E : 
    This function goes through each 'drugXXX_interaction.txt' file of the 'interaction' folder and asses if the 
    real target is in the candidate genes (meaning the genes in the 'drugXXX_interaction.txt' files). 
    An example usage is provided at the end. 

    I N P U T : 
    --> myTargets_Recon : a tab-delimited text file containing the associations between the drug screens names and 
        their associated gene targets (in BiGG id). The Bigg ids should be in the first column, and the name of the screen
        (without any suffix such as '_hits' or 'hits_targets') in the second column. 
        An example is provided as 'myTargets_Recon.txt'. 
    --> directory_interactions : the directory of the folder containing all the interaction files ('drugXXX_interaction.txt'). The genes in those files
        are the candidate target genes. Those files can be issues from the 'get_interactions.py' script. 

    O U T P U T : 
    --> eligible_screens_to_p_value : a list containing the name of the screens that fit the following criteria : 
        - the screen name exists both in the myTargets_Recon and in the 'interaction' folder
        - the screen has one or more hits
        - the hit(s) of the screens are in the model. 
    --> totsucess : the number of screens for which the 'true target' (the one listed in 'myTargets_Recon') was
        found in the candidate genes of that screen (i.e. found in 'drugXXX_interactions.txt')
    --> screenDict : a dictionnary built from 'myTargets_Recon', containing the name of the screen as keys and 
        the associated 'true targets' as values. 

    '''

    def is_float(string): 
        try : 
            float(string)
            return True
        except : 
            return False

    # make a dictionary of every drug screen and it's associated targets 
    screenList = []
    with open(myTargets_Recon, 'r') as file : 
        for lines in file :
            if lines.split('\t')[1] not in screenList: 
                screenList.append(lines.split('\t')[1])
    screenDict = {}
    for i in range(len(screenList)): 
        myScreen = screenList[i]
        myTargets = []
        with open(myTargets_Recon, 'r') as file : 
            for lines in file :
                if lines.split('\t')[1] == myScreen : 
                    myTargets.append(lines.split('\t')[0])
                screenDict[myScreen[:-1]] = myTargets

    # now we find the target. In each of the eligible screen having a target in our model, we check if the target is 
    # in the candidate targets from the 'drugXXX_interaction.txt' files. 
    import os 

    # to transform '['36_AT1', '36_AT2', '36_AT3']' into '36_AT1'
    from converter5 import allToat
    allToatDict = allToat()
    

    # construct list of all screens with known recon targets 
    reconScreens = list(screenDict.keys())
    

    # go through all screens' interaction file in interaction folder, see if target in candidate genes 
    eligible_screens_to_p_value = []
    totsucess =0
    for fileName in os.listdir(directory_interactions):
        sucess = 0
        #change to fileName[:-x] depending on name of file. screenName should be the name of the screen only, without the '_interactions.txt'
        screenName = fileName[:-19] 
        if screenName in reconScreens: 
            filePath = os.path.join(directory_interactions, fileName)
            eligible_screens_to_p_value.append(fileName)
            with open(filePath, 'r') as file : 
                for lines in file : 
                    if is_float(lines.split('\t')[0][0]) == True: 
                        candidateGene = allToatDict[lines.split('\t')[1]]
                        targets = screenDict[screenName]
                        if candidateGene in targets: 
                            sucess += 1
        if sucess > 0 : 
            totsucess += 1
      
    return eligible_screens_to_p_value, totsucess, screenDict


if __name__ == '__main__':
    myTargets_Recon =  'myTargets_Recon.txt'
    directory_interactions = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'
    getWorkingScreens(myTargets_Recon, directory_interactions)