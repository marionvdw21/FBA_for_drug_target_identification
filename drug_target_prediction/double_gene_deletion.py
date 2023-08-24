def double_gene_deletion(fileName, myModel):

    '''
    U S A G E : 
    When given a list of gene (crispr hits) and a model, this function will perform double gene deletion of every gene in the list of 
    gene VS all the other genes in the model. Every modified model (that underwent double gene deletion) is then optimized for maximization 
    of biomass production. 

    I N P U T : 
    --> fileName : the name of the file containing all the CRISPR hits that should undergo double gene deletion. It could be outputed from 
        the 'get_hits.py' script, as 'list_of_hits.txt'. It should be a tab-deliited txt file with the gene's Bigg id in first column
        An example is provided as 'list_of_hits.txt'. 
    --> myModel : the metabolic model, of class cobra.core.Model. 

    O U T P U T : 
    --> 'GENEXXX_VS_all.txt'
        For each gene in the 'fileName' file, a file with the name 'GENEXXX_VS_all.txt' wil be outputed. 
        This file will contain the result (biomass production after optimization) of the single gene deletion of that gene, as well as the result
        of the double gene deletion of that gene with every other gene in the model. 
        It will be of format : 
        _______________________________________________________________________________
        |   Wild Type (no gene deletion)        biomass production flux                | 
        |   -----------                                                                |
        |   [gene id of crispr hit deleted (first KO)]     biomass production flux     |
        |   -----------                                                                |
        |   [gene id of gene1 (second KO)]  biomass production flux                    |
        |   [gene id of gene2 (second KO)]  biomass production flux                    |
        |   [gene id of gene3 (second KO)]  biomass production flux                    |
        |   [gene id of gene4 (second KO)]  biomass production flux                    | 
        |                                                                              |
        |   ...                                                                        |
        |                                                                              |
        |   [gene is of genex (secondKO)]   biomass production flux                    |
        |______________________________________________________________________________|   
    '''



    import cobra.flux_analysis 
    import time 
    import cobra.manipulation
    import cobra.flux_analysis.reaction
    import convertMyModel
    
    model = myModel
    model.optimize()
    solNoDel = model.reactions.BIOMASS_maintenance.flux
    string1 = ''
    string1 += 'WT' + '\t' + str(solNoDel) #line for original model
    header = 'second genes knocked out' + '\t' + 'sol_double_KO' + '\t' + 'sol_double_KO/sol_single_KO'


    # making a gene list of type  [[26_AT1], [314_AT1, 314_AT2], ...] (list of lists) 
    # that will be used for the first knock out --> contains only hit genes  
    listHitGenes = []
    with open(fileName, 'r') as file: 
        for lines in file: 
            line_split = lines.split('\t')
            hit_name = line_split[0]
            listHitGenes.append(hit_name)
    listHitGenesAllAT = []
    for gene in listHitGenes: 
        AT = 1
        allATs=[]
        allATs.append(gene)
        geneName = gene[:-3]
        for i in range(20): 
            AT +=1 
            gene2 = '%(geneName)sAT%(AT)s' %{'geneName' : geneName, 'AT': AT}
            if gene2 in model.genes: 
                allATs.append(gene2)
        listHitGenesAllAT.append(allATs)
    
    # making a gene list of type  [[26_AT1], [314_AT1, 314_AT2], ...] (list of lists) 
    # that will be used for the second knock out --> contains all genes in model
    geneListraw = []
    for a in range(len(model.genes)): 
        geneListraw.append(str(model.genes[a]))
    geneList = []
    for gene in geneList : 
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
            geneList.append(allATs)


    #---------------------------------------#
    # ---------- f i r s t    k o  ---------#
    #---------------------------------------#
    #looping over hit genes for 1st KO
    for j in range(len(listHitGenesAllAT)): #change it to for j in range(2), or another integer if it is desired to do this long task in segments
        newLines = []
        newLines.append(string1)
        newLines.append('-------')
        newModel = model.copy()
        KOgenes = listHitGenesAllAT[j]
        #using technique #1 : TOTAL KO
        for KOgene in KOgenes : 
            genesrxn = newModel.genes.get_by_id(KOgene).reactions
            for reaction in genesrxn : 
                reaction.upper_bound = 0
                reaction.lower_bound = 0
        newModel.optimize()
        solKO1 = newModel.reactions.BIOMASS_maintenance.flux
        string = ''
        string += str(KOgenes) + '\t' + str(solKO1) + '\t' + str(solKO1/solNoDel)
        newLines.append(string)
        newLines.append('---------')
        newLines.append(header)
        #---------------------------------------#
        # ---------- s e c o n d   k o  --------#
        #---------------------------------------#
        
        for k in range(len(geneList)): #change it for range(len(geneList))
            newModel2 = newModel.copy()
            KOgenes = geneList[k]
            #technique 1 : total KO
            for KOgene in KOgenes: 
                genesrxn = newModel2.genes.get_by_id(KOgene).reactions
                for reaction in genesrxn :
                    reaction.upper_bound = 0
                    reaction.lower_bound = 0 
            newModel2.optimize()
            solKO2 = newModel2.reactions.BIOMASS_maintenance.flux
            string = ''
            string += str(KOgenes) + '\t' + str(solKO2) + '\t'+ str(solKO2/solKO1)
            newLines.append(string)

        OutputFileName = '%(HitGene)s_VS_all.txt' %{'HitGene': listHitGenes[j]}

        with open(OutputFileName, 'w') as file: 
            for lines in range(len(newLines)): 
                file.write(newLines[lines])
                file.write('\n')
           
    return

if __name__ == '__main__':
    # example of usage 
    import cobra.io as co
    import convertMyModel
    import convertMyMedia
    myModel = convertMyModel.convertMyModel('expression_Recon_ordered2.txt', co.read_sbml_model('Recon3.xml'))[0]
    myModel = convertMyMedia(myModel, 'minimal_media.txt')
    double_gene_deletion('list_of_hits.txt', myModel)