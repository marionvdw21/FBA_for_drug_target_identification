
def atToall():

    '''
    U S A G E : 
    This function create a dictionary in order to go from Bigg IDs written in this manner: 
         '['36_AT1', '36_AT2', '36_AT3']' 
    to BiGG IDs written in this manner : 
          '36_AT1'

    
    N O T E  : 
    The model is assumed to be Recon3, but that can be changed by changing the lines in the file 'atToallDict.txt'. 

    '''

    atToallDict = {}
    with open('atToallDict.txt', 'r') as file : 
         for lines in file : 
              atToallDict[lines.split('\t')[0]] = lines.split('\t')[1]
         
    return atToallDict

def allToat(): 

    allToatDict = {}
    atToallDict = atToall()
    for key in atToallDict:
         newKey = atToallDict[key]
         allToatDict[str(newKey)] = key
    return allToatDict

if __name__ == '__main__': 
     atToallDict = atToall()
     allToatDict = allToat()
