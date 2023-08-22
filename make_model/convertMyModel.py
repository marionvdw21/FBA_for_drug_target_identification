def convertMyModel(fileName, myModel): 

    '''
    U S A G E : 
    The usage is to modify a generic metabolic model, such as Recon3, so that is it resembles a specific cell line, 
    of which the expression file is provided. The way the generic model is modified is by changing the upper and lower 
    bounds of the reaction fluxes, accoding to the exponential function ub = (0.1*exp(0.6908*t))/100, in which 
    ub is the percentage of the original bound that will be applied to the new model, and t is the expression of the 
    associated gene. 
    An example of usage is provided at the end. 

    I N P U T S : 
    --> fileName : the expression file of the cell line. An example is provided as 'expression.txt'. As some genes share 
    the same reactions, the order in which the reaction bounds is modified matters. 
        choices of input file (fileName): 
        - expression_Recon.txt 
        - expression_Recon_ordered.txt 
        - expression_Recon_deordered.txt 

    --> myModel : generic model, which will be modified, of class cobra.core.Model.

    O U T P U T S : 
    --> adapted_model : class cobra.core.Model. It is the input model with the reaction flux boundaries modified according to the 
        expression file. 
    --> geneList : a list of genes of the adapted model, of form [[26_AT1], [314_AT1, 314_AT2], [6768_AT1], [332_AT1, 332_AT2, 332_AT3], ...]
    --> sol : the solution of the modified model once it is optimized for maximum biomass production. 
        The solution corresponds to the flux passing through the reaction 'BIOMASS_maintenance'. 
    '''

    def is_float(string): 
        try: 
            float(string)
            return True
        except ValueError: 
            return False 
    
    import cobra 
    import cobra.flux_analysis 
    import time 
    import cobra.manipulation 
    import cobra.flux_analysis.reaction 
    import cobra.io as co
    import copy 
    from numpy import exp
    from matplotlib import pyplot as plt

    model = myModel
    model.solver = 'glpk'
    model.optimize()

    with open(fileName, 'r') as table : 
        gene_names = []
        expression = []
        for lines in table : 
            lines_split = lines.split('\t')
            gene_names.append(lines_split[0])
            expression.append(lines_split[1])
    
    geneList = [] #geneList is a list of form  [[26_AT1], [314_AT1, 314_AT2], ...]
    for gene in gene_names:
        AT=1
        allATs = []
        allATs.append(gene)
        geneName = gene[:-3]
        for i in range(20): 
            AT +=1 
            gene2 = '%(geneName)sAT%(AT)s' %{'geneName' : geneName, 'AT':AT}
            if gene2 in model.genes:
                allATs.append(gene2)
        geneList.append(allATs)

    adapted_model = model.copy()
    expression_list = []
    name_list = []
    flux_lmt_list = []
    for i in range(len(geneList)): 
        if is_float(expression[i]) == True: 
            t = float(expression[i])
            # the upper and lower bounds are capped according to this exponential function, which can be changed and rearranged.
            # the new upper and lower bounds will be capped at x% of the original bounds (most of them being -1000, 1000). 'x' being the value of 'ub' 
            ub = (0.1*exp(0.6908*t))/100 # x%.
            name_list.append(geneList[i])
            expression_list.append(float(expression[i])) #for the plot
            flux_lmt_list.append(ub) # for the plot
            KOgenes = geneList[i]
            for KOgene in KOgenes: 
                generxn = adapted_model.genes.get_by_id(KOgene).reactions
                for reaction in generxn: 
                    if reaction.lower_bound > model.reactions.get_by_id(reaction.id).upper_bound * ub: 
                        reaction.lower_bound = model.reactions.get_by_id(reaction.id).upper_bound * ub
                    reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*ub 
                    reaction.lower_bound = model.reactions.get_by_id(reaction.id).lower_bound*ub

    adapted_model.slim_optimize()
    sol = str(adapted_model.reactions.BIOMASS_maintenance.flux)

    # uncomment the following section if a plot of the new model's fluxes and bounds is desired. 
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~ ~ ~ ~ p l o t t i n g ~ f l u x e s  ~ ~ ~ ~ #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    x =  expression_list
    y = []     
    
    for i in range(len(name_list)): 
        name = name_list[i]
        myList = []
        for name in name_list[i]: 
            generxn = adapted_model.genes.get_by_id(name).reactions
            for reaction in generxn: 
                myList.append(adapted_model.reactions.get_by_id(reaction.id).flux)
        y.append(myList)

    for xe, ye in zip(x, y):
        plt.scatter([xe]*len(ye), ye, s=1, c='g')
    plt.show()
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~ ~ ~ ~ ~ p l o t t i n g ~ l i m i t s ~ ~ ~ ~ ~ # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    plt.plot(expression_list, flux_lmt_list, '.')
    plt.text(0, 100, sol[:7], bbox = dict(facecolor='yellow', alpha=0.5))
    plt.xlabel('gene expression in NALM6 cells (log2 read counts)') # change 'NALM6' for whichever cell line the expression file corresponds to
    plt.ylabel('flux limit in NALM6 model (percent of Recon3 bound)')  # change 'NALM6' for whichever cell line the expression file corresponds to
    plt.title('NALM6 cobra model based on gene expression')   # change 'NALM6' for whichever cell line the expression file corresponds to
    plt.grid(which='both')
    plt.show()
    '''
    
    return adapted_model, geneList, sol

if __name__ == '__main__':
    # example of usage: 
    import cobra.io as co
    convertMyModel('expression_Recon_ordered2.txt', co.read_sbml_model('Recon3.xml'))