import cobra 
import cobra.flux_analysis 
import time 
import cobra.manipulation 
import cobra.flux_analysis.reaction 
import cobra.io as co
import copy 
from GetEssRxn2 import GetEssRxn
from convert2NALM import *

def KOALL3():
    import cobra
    import cobra.flux_analysis 
    import time 
    import cobra.manipulation
    import cobra.flux_analysis.reaction
    import cobra.io as co
    import copy
    model = convertToNalm()[0]
    geneList = convertToNalm()[1] #list of type  [[26_AT1], [314_AT1, 314_AT2], ...] (list of lists)
    totTotTime = 0
    string = ''
    string += 'WT' + '\t' + str(model.optimize()) #line for original model
    #newLines.append(string)
    header = 'second genes knocked out' + '\t' + 'BIOMASS_maintenance flux using 10% technique' + '\t' + 'time' + 'Growth ratio doubleDel:singleDel' #header
    solNoDel = 518.140 #CHANGE WHEN USING NEW NALM6
    newModel = model.copy()

    #first looping over ''essential genes''
    essGenes = GetEssRxn('allSingleGeneDelNALM6.txt')[0] #[26_AT1, 304_AT1, 5069_AT1, ...]
    essGenesAllAT = [] # #--> creating nice list of essential genes of type  [[26_AT1], [314_AT1, 314_AT2], ...] (list of lists)
    for gene in essGenes:
        AT=1
        allATs = []
        allATs.append(gene)
        geneName = gene[:-3]
        for i in range(20): 
            AT +=1
            gene2 = '%(geneName)sAT%(AT)s' %{'geneName':geneName, 'AT':AT}
            if gene2 in model.genes: 
                allATs.append(gene2)
        essGenesAllAT.append(allATs)
    # looping over essential genes : first 'KO'
    for i in range(len(essGenesAllAT)):#change it to (essGenesAllAT)
        print(i, '--------------')
        newLines = []
        newModel = model.copy() 
        KOgenes = essGenesAllAT[i]
        for KOgene in KOgenes : 
            generxn = newModel.genes.get_by_id(KOgene).reactions
            for reaction in generxn: 
                reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*0.1
        newModel.optimize
        solSingleKO = newModel.reactions.BIOMASS_maintenance.flux
        string = 'SingleKO: ' + str(KOgenes) + '\t' + str(solSingleKO) + str(solSingleKO/solNoDel)
        newLines.append(string)

        #looping over all genes > looping over essential genes : second 'KO': 
        for j in range(2): #change to len(geneList) 
            newModel2 = newModel.copy()
            start_time = time.time()
            print(i, j)
            string = ''
            KOgenes = geneList[j]
            for KOgene in KOgenes : 
                print(KOgene)
                generx = newModel2.genes.get_by_id(KOgene).reactions
                for reaction in generxn: 
                    reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*0.1
            newModel2.optimize()
            sol = newModel2.reactions.BIOMASS_maintenance.flux
            string += str(KOgenes) + '\t' + str(sol) + '\t' + str(time.time() - start_time) + str(sol/solSingleKO)
            newLines.append(string)

        OutputFileName = 'DoubleKO_%(EssGene)s_VS_all.txt' %{'EssGene' : str(essGenesAllAT[i][0])}
        with open(OutputFileName, 'w') as file: 
            for lines in range(len(newLines)): 
                file.write(str(newLines[lines]))
                file.write('\n')

            
    return
KOALL3()