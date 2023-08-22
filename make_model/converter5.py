import cobra.io as co

def atToall():
    
    human = r"C:\Users\mario\OneDrive\Documents\SURE\Recon3.xml"
    model = co.read_sbml_model(human)

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