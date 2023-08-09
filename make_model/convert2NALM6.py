import cobra 
import cobra.flux_analysis 
import time 
import cobra.manipulation 
import cobra.flux_analysis.reaction 
import cobra.io as co
import copy 
from numpy import exp
from matplotlib import pyplot as plt

'''
iteration de CONVERT2NALM ou j'utilise une fonction exponentielle comme limites de flux
la meilleure a date


    choices of input file : 
    --> expression_Recon.txt --> allSingleGeneDel_NALM10
    --> expression_Recon_ordered2.txt --> allSingleGeneDel_NALM11
    --> expression_Recon_deordered2.txt --> allSingleGeneDel_NALM12 (la meilleure IMO)
    

'''

human = r"C:\Users\mario\OneDrive\Documents\SURE\Recon3.xml"
model = co.read_sbml_model(human)
model.solver = 'glpk'
model.optimize()

def is_float(string): 
    try: 
        float(string)
        return True
    except ValueError: 
        return False 
    
def convertToNalm6(fileName): 
    
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

    NALM6 = model.copy()
    true = 0
    expression_list = []
    name_list = []
    flux_lmt_list = []
    for i in range(len(geneList)): 
        if is_float(expression[i]) == True: 
            t = float(expression[i])
            ub = (0.1*exp(0.6908*t))/100 # x%
            name_list.append(geneList[i])
            expression_list.append(float(expression[i])) #for the plot
            flux_lmt_list.append(ub) # for the plot
            KOgenes = geneList[i]
            for KOgene in KOgenes: 
                generxn = NALM6.genes.get_by_id(KOgene).reactions
                for reaction in generxn: 
                    if reaction.lower_bound > model.reactions.get_by_id(reaction.id).upper_bound * ub: 
                        reaction.lower_bound = model.reactions.get_by_id(reaction.id).upper_bound * ub
                    reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*ub 
                    reaction.lower_bound = model.reactions.get_by_id(reaction.id).lower_bound*ub

    NALM6.slim_optimize()
    sol = str(NALM6.reactions.BIOMASS_maintenance.flux)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~ ~ ~ ~ p l o t t i n g ~ f l u x e s  ~ ~ ~ ~ #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #x =  expression_list
    #y = []     
    
    #for i in range(len(name_list)): 
    #    name = name_list[i]
    #    myList = []
    #    for name in name_list[i]: 
    #        generxn = NALM6.genes.get_by_id(name).reactions
    #        for reaction in generxn: 
    #            myList.append(NALM6.reactions.get_by_id(reaction.id).flux)
    #    y.append(myList)

    #for xe, ye in zip(x, y):
    #    plt.scatter([xe]*len(ye), ye, s=1, c='g')
    #plt.show()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~ ~ ~ ~ ~ p l o t t i n g ~ l i m i t s ~ ~ ~ ~ ~ # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #plt.plot(expression_list, flux_lmt_list, '.')
    #plt.text(0, 100, sol[:7], bbox = dict(facecolor='yellow', alpha=0.5))
    #plt.xlabel('gene expression in NALM6 cells (log2 read counts)')
    #plt.ylabel('flux limit in NALM6 model (percent of Recon3 bound)')
    #plt.title('NALM6 cobra model based on gene expression')  
    #plt.grid(which='both')
    #plt.show()
    
    return NALM6, geneList, sol

# only this works (aka no forced flux). But this is weird, maybe issue with solver capacity

'''
unordered (expression_Recon.txt) 
3989 (other method : it's 4136)--> above 9
1529 (other method : it's 1574)--> above 10
BIOMASS_maintenance : 389.1870204422617

ordered croissant (expression_Recon_ordered2.txt):
1545 --> above 10
BIOMASS_maintenance : 607.3963558610566

ordered decroissant (expression_Recon_deordered2.txt)
1510 --> above 10
BIOMASS_maintenance : 123.88704956434468
'''