from converter5 import *
import scipy.stats 
import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def fisherTest(inputFileList, theCutoff, expCutoff, alternative): 

    allToatDict = allToat()
    essNames = []
    essRank = []

    

    with open('essentiality_Recon.txt', 'r') as file : 
        for lines in file : 
            essNames.append(lines.split('\t')[0])
            essRank.append(lines.split('\t')[1])
    i = 0
    for inputFile in inputFileList:
        cluster = []
        i += 1
        print(inputFile)
        countList, essRankList = [], []
        EE, EEn, EN, ENn, NE, NEn, NN, NNn = [], [], [], [], [], [], [], []
        with open(inputFile, 'r') as file : 
            count = -1
            for lines in file: 
                count += 1
                name = lines.split('\t')[0]
                name = allToatDict[name]
                if name in essNames : 
                    temp = []
                    countList.append(count)
                    index = essNames.index(name)
                    essRankList.append(float(essRank[index]))
                    if count < theCutoff: # Essential for theoretical 
                        if float(essRank[index]) < expCutoff: # Essential for theoretical AND experimental
                            temp.append(count)
                            temp.append(essRank[index])
                            EE.append(temp)
                            EEn.append(name)
                        else : #Essential for theoretical, non essential for experimental
                            temp.append(count)
                            temp.append(essRank[index])
                            EN.append(temp)
                            ENn.append(name)
                    else : # Non essential for therotical
                        if float(essRank[index]) < expCutoff: #non essential for theoretical and essential for experimental
                            temp.append(count)
                            temp.append(essRank[index])
                            NE.append(temp)
                            NEn.append(name)
                        else : # non essential for theoretical and non essential for experimental
                            temp.append(count)
                            temp.append(essRank[index])
                            NN.append(temp)
                            NNn.append(name)
                    #if inputFile == 'EssPred-MOMA-NALM.txt' and count < 1700 and count > 1000 and float(essRank[index]) < -3: 
                    #        cluster.append(name)
                    
                    
        table = np.array([[len(EEn), len(ENn)], [len(NEn), len(NNn)]])
        print(table)
        res = scipy.stats.fisher_exact(table=table, alternative=alternative)
        print(res[1])
        print(inputFile, res)
        print('len NN : ', len(NNn))
        print('len NE : ', len(NEn))
        print('len EN : ', len(ENn))
        print('len EE : ', len(EEn))
        print('---------')
        
        #plt.subplot(2, 2, i)
        ax = plt.gca()
        plt.scatter(countList, essRankList, s = 1)
        plt.axhline(y=expCutoff, color ='r')
        plt.axvline(x=theCutoff, color='r')
        #ax.set_yscale('symlog')
        plt.title(inputFile[8:-4])
        pvalue = res[1]
        #pvalue = round(pvalue, 35)
        text = 'pvalue: ' + str(pvalue)
        plt.annotate(text, (300, 2))
        #if i == 2 or i == 1 : 
        #    ax.get_xaxis().set_visible(False)
        if i == 2 or i == 4 : 
            ax.get_yaxis().set_visible(False)
        title = 'cutoffs: ' + str(theCutoff) + ' and ' + str(expCutoff) + ', ' + str(alternative)
        plt.suptitle(title)
    plt.show()

    #with open('cluster.txt', 'w') as file : 
    #    for i in range(len(cluster)): 
    #        file.write(str(cluster[i]))
    #        file.write('\t')
    #        file.write('\n')

    return 

#fileList = ['EssPred-FBA-WT.txt', 'EssPred-FBA-NALM.txt', 'EssPred-MOMA-WT.txt', 'EssPred-MOMA-NALM.txt']
fileList = ['EssPred33.txt']
fisherTest(fileList, 200, -4, 'two-sided')
# alternative : 'two-sided' or 'greater' or 'less'


'''
fileList = ['EssPred-FBA-NALM.txt'] --> categrical nalm
fisherTest(fileList, 75, -4, 'two-sided')
nalm_fisher.png

fileList = ['EssPred-FBA-WT.txt']
fisherTest(fileList, 56, -4, 'two-sided')
wt_fisher.png

fileList = ['EssPred_NALM11.txt']
fisherTest(fileList, 100, -4, 'two-sided')
exp_nalm_fisher.png
'''