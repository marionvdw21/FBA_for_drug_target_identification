def convert_target_file(fileName, model):

    '''
    U S A G E : 
    This function converts a file containing drugs screen names and their associated targets to another file (named 'myTargets_Recon_example.txt') 
    that can be read by the script 'getWorkingScreens.py'. 

    I N P U T S : 
    --> fileName : the name of the tab delimited text file containing the name of the drug screens and their associated targets.
        The targets should be the name of the gene, in HGNC symbols. It should be in the following format : 
        __________________________________________________________
        |   screen name     drug name   target|target|target      |
        |   screen name     drug name   target                    |
        |   screen name     drug name   target|target|target      |
        |   ...                                                   |
        |   screen name     drug name   target|target             | 
        |_________________________________________________________|

        An example is provided as 'target_table_example.txt'

    --> model : the metabolic model of class cobra.core.Model. 

    O U T P U T S : 
    --> 'myTargets_Recon.txt' : a tab delimited txt file that can be used in the 'getWorkingScreens.py' script. 
        The fist column contains molecular targets in BiGG id, and the second column contains the name of their associated screen. 
        An example of such as file is provided under 'myTargets_Recon_example.txt'. 
    '''

    import HGNC_to_Recon
    screenDict = {}
    with open(fileName, 'r') as file : 
        for lines in file :
            screenDict[lines.split('\t')[0]] = lines.split('\t')[2][:-1].split('|')
    screenName = list(screenDict.keys())
    screenTargets = list(screenDict.values())
            
    with open('myTargets.txt', 'w') as file : 
        for i in range(len(screenName)): 
            for j in range(len(screenTargets[i])): 
                file.write(screenTargets[i][j])
                file.write('\t')
                file.write(screenName[i])
                file.write('\n')

    HGNC_to_Recon.HGNC_to_Recon_converter('myTargets.txt', model)

    return 


if __name__ == '__main__': 
    # example of usage
    import cobra.io as co
    convert_target_file('target_table.txt', co.read_sbml_model('Recon3.xml'))
