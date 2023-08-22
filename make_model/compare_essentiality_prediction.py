def fisherTest(inputFileList, experimental_essentiality, theCutoff, expCutoff, alternative): 

    '''
    U S A G E : 
    This function compare the gene essentiality predicted from single gene deletion to the experimental gene essentiality
    An example of usage is provided at the end. 

    I N P U T : 
    --> inputFileList : a list containing the name(s) of the input file. The input file should be a file outputed from the 'Essentiality_prediction' script, 
        with the most essential gene first and least essential last. 
        An example is provided as 'EssPred.txt'
    --> experimental_essentiality : a tab delimited text file with each line containing the name of the gene and the rank score for essentiality. 
        If the original experimental essentiality file has the gene names as HGNC symbols, it can be passed through the 
        'HGNC_to_Recon' script. 
        An example is provided as 'essentiality_Recon.txt'. 
    --> theCutoff : the cuttof delimiting essential from non-essential genes, regarding the predicted essentiality using single gene deletion (files in 'inputFileList'). 
        If theCutoff = 200, then the genes positioned from 0-199 (in the files from inputFileList) will be considered essential, 
        and the positionned below or equal to 200 will be considered non-essential. 
        It is recommanded to use a cutoff above which the deletion of the associated genes impact minimally the biomass production. 
    --> expCutoff : the cutoff delimiting essential from non-essential genes, regading the experimental essentiality (experimental essentiality file). 
        Expcutoff is a rank score value. For example if the cutoff is -4, then genes having a rank score below -4 will be considered 
        essential, and the genes having a rank score above -4 will be considered non-essential. 
    --> alternative : 'below', 'greater', or 'two-sided'. The nature of the p-value that will be calculated. 

    O U T P U T : 
    --> p-value : the p-value resulting from a fisher test on the 2x2 contigency table : 
        [[predicted essential and experimentally essential,     predicted essential but experimentally non-essential]
        [predicted non-essential but experimentally essential,  predicted non-essential and experimentally non-essential]]
    --> A figure showing the distribution of each gene, according to its predicted and experimental essentiality. 

    '''
    import converter5
    import scipy.stats 
    import numpy as np
    from matplotlib import pyplot as plt
    import matplotlib
    allToatDict = converter5.allToat() # this lines makes the gene name go from [36_AT1, 36_AT2, 36_AT3] to '36_AT1' so that it is cmparable to the experimental essentiality. 
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
    # example of usage : 
    fileList = ['EssPred.txt']
    fisherTest(fileList, 'essentiality_Recon.txt', 200, -4, 'two-sided')


