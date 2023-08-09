import numpy as np
import scipy.stats

def fisherModel(model1, model2): 
    with open('fisherModels.txt') as file : 
        for lines in file : 
            if lines.split('\t')[0] == model1 : 
                EE1 = lines.split('\t')[1]
                EN1 = lines.split('\t')[2]
                NE1 = lines.split('\t')[3]
                NN1 = lines.split('\t')[4]
            if lines.split('\t')[0] == model2: 
                EE2 = lines.split('\t')[1]
                EN2 = lines.split('\t')[2]
                NE2 = lines.split('\t')[3]
                NN2 = lines.split('\t')[4]
    #table = np.array([[len(EEn), len(ENn)], [len(NEn), len(NNn)]])
    table1 = np.array([[EE1, EN1], [EE2, EN2]])
    table2 = np.array([[EE1, NE1], [EE2, NE2]])
    res1 = scipy.stats.fisher_exact(table=table1, alternative='two-sided')[1]
    res2 = scipy.stats.fisher_exact(table2, alternative = 'two-sided')[1]
    print('p-value 1 : ', res1)
    print('p-value 2 : ', res2)
    return 

#fisherModel('EssPred18.txt', 'EssPred-FBA-NALM.txt')
fisherModel('EssPred-FBA-NALM.txt', 'EssPred33.txt')

'''
results : 

**p-value 1 comes from table 1 : table1 = np.array([[EE1, EN1], [EE2, EN2]])
**p-value 2 comes from table 2 : table2 = np.array([[EE1, NE1], [EE2, NE2]])

fisherModel('EssPred-FBA-WT.txt', 'EssPred-FBA-NALM.txt')
p-value 1 :  0.011146051376535314
p-value 2 :  0.0055851088499088855

fisherModel('EssPred-FBA-NALM.txt', 'EssPred18.txt')
p-value 1 :  0.4730717099534844
p-value 2 :  0.4678637293532134

fisherModel('EssPred-FBA-WT.txt', 'EssPred_NALM12.txt')
p-value 1 :  0.06860082987782104
p-value 2 :  0.040528846088505786

fisherModel('EssPred-FBA-WT.txt', 'EssPred18.txt')
p-value 1 :  0.08592120251337107
p-value 2 :  0.05189725128541929

fisherModel('EssPred-FBA-NALM.txt', 'EssPred15.txt')
p-value 1 :  0.07729650548050473
p-value 2 :  0.07420962233015409

fisherModel('EssPred22.txt', 'EssPred23.txt')
p-value 1 :  0.0841452784972615
p-value 2 :  0.08158772205834547

fisherModel('EssPred23.txt', 'EssPred20.txt')
p-value 1 :  0.7384153786067779
p-value 2 :  0.736617903498418

fisherModel('EssPred_FBA_NALM.txt', 'EssPred11.txt') (categories vs gradient, prioritizing high expression)
p-value 1 :  0.021148220918024006
p-value 2 :  0.019887833427901676

fisherModel('EssPred-FBA-NALM.txt', 'EssPred33.txt') (categories vs gradient, prioritizinv high expression, with minimal_media)
p-value 1 :  0.01113609335499436
p-value 2 :  0.14315515767884113
'''