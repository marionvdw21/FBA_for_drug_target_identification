def screen_to_gene_p_value(directory_interactions, directory_hits_recon, Ho_tot): 
    
    '''
    U S A G E : 
    This function takes the results from 'get_inverse_signature.py' and associates a p-value to it, to verify if the results of the 
    'get_inverse_signature.py' script are not attributable to chance. 
    An example of usage is provided at the end. 

    I N P U T S :  
    --> directory_interactions : used as an input of the 'get_inverse_signature.py' script. See that script for more description on that input.
    --> directory_hits_recon : used as an input of the 'get_inverse_signature.py' script. See that script for more description on that input. 
    --> Ho_tot : the 'n' variable for the p-value calculation. The larger the 'n', the more precise the p-value. 

    O U T P U T S : 
    --> a p-value : probability of getting the results of 'get_inverse_signature.py' by chance

    N O T E : 
    Some gene are associated to a great number of screens, which makes them non-specific to the screen with which they are truly interacting, and add noise
    to the data. Therefore, only the gene interacting with less than half of the total number of eligible screens are considered to calculate this p-value. 
    This value can be changed, simply change the number in the script 'get_inverse_signature.py' (0.5*x). 
    '''
    
    import get_inverse_signature
    import os 
    import random

    threshold = get_inverse_signature.inverseSignature(directory_interactions, directory_hits_recon)[0]
    allcandidateDict = get_inverse_signature.inverseSignature(directory_interactions, directory_hits_recon)[1]

    
    Ha = 0


    # creating my pool of screens from which I will pick
    screenPool = []
    for fileName in os.listdir(directory_hits_recon): 
        screenPool.append(fileName)


    # process of picking randomly from my pool of screen. 
    allcandidate_keys = list(allcandidateDict.keys()) # list of genes 
    allcandidate_values = list(allcandidateDict.values()) # list of screens with which gene is interacting 
    for i in range(Ho_tot):
        sucess = 0
        for k in range(len(allcandidate_keys)): # going through each gene 
            n = len(list(allcandidate_values[k]))
            for j in range(n): # going through each screen with which this gene is associated. 
                picker = random.randint(0, len(screenPool)-1)
                pickedScreen = screenPool[picker]
                filePath = os.path.join(directory_hits_recon, pickedScreen)
                with open(filePath, 'r') as file : 
                    realHits = []
                    for lines in file : 
                        realHits.append(lines.split('\t')[0])
                if allcandidate_keys[k] in realHits : 
                    sucess += 1

        if sucess > threshold : 
            Ha += 1                


    print('p-value: ', Ha/Ho_tot)
    p_value = Ha/Ho_tot
    return p_value

if __name__ == '__main__': 
    # example of usage
    directory_interactions = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'
    directory_hits_recon = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    screen_to_gene_p_value(directory_interactions, directory_hits_recon, 10)

