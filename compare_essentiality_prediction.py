def fisherTest(inputFileList, experimental_essentiality, theCutoff, expCutoff, alternative): 

    import converter5
    import scipy.stats 
    import numpy as np
    from matplotlib import pyplot as plt
    import matplotlib
    allToatDict = converter5.allToat()
    essNames = []
    essRank = []

    with open(experimental_essentiality, 'r') as file : 
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
        res = scipy.stats.fisher_exact(table=table, alternative=alternative) # alternative : 'two-sided' or 'greater' or 'less'
        print(res[1])
        print(inputFile, res)
        print('len NN : ', len(NNn))
        print('len NE : ', len(NEn))
        print('len EN : ', len(ENn))
        print('len EE : ', len(EEn))
        print('---------')
        

        ax = plt.gca()
        plt.scatter(countList, essRankList, s = 1)
        plt.axhline(y=expCutoff, color ='r')
        plt.axvline(x=theCutoff, color='r')
        pvalue = res[1]
        text = 'pvalue: ' + str(pvalue)
        plt.annotate(text, (300, 2))
        title = 'cutoffs: ' + str(theCutoff) + ' and ' + str(expCutoff) + ', ' + str(alternative)
        plt.title(title)
    plt.show()

    return 

if __name__ == '__main__': 
    fileList = ['EssPred2.txt']
    fisherTest(fileList, 'essentiality_Recon.txt', 200, -4, 'two-sided')


