
def atToall():
    import cobra.io as co
    '''
    U S A G E : 
    This function create a dictionary in order to go from Bigg IDs written in this manner: 
         '['36_AT1', '36_AT2', '36_AT3']' 
    to BiGG IDs written in this manner : 
          '36_AT1'

    The model is assumed to be Recon3, but that can be changed at line 15. 
    '''

    model = co.read_sbml_model('Recon3.xml')

    geneList = [] # geneList : [26_AT1, 103_AT1, 103_AT2, ...]
    for i in range(len(model.genes)): 
        geneList.append(str(model.genes[i]))

    atToallDict = {}
    for gene in geneList: 
        if gene[-3:] == 'AT1': 
                allATs = []
                allATs.append(gene)
                geneName = gene[:-3]
                AT = 1
                for i in range(20): 
                    AT += 1
                    gene2 = '%(geneName)sAT%(AT)s' %{'geneName': geneName, 'AT': AT}
                    if gene2 in geneList : 
                        allATs.append(gene2)
                atToallDict[gene] = allATs
    
    return atToallDict

def allToat(): 
    allToatDict = {}
    atToallDict = atToall()
    for key in atToallDict:
         newKey = atToallDict[key]
         allToatDict[str(newKey)] = key
    return allToatDict

if __name__ == '__main__': 
     atToall()
     allToat()