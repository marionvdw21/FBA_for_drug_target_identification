def interactions(cutoff): 
    import os 
    directory = r'C:\Users\mario\OneDrive\Documents\SURE\CRANKS_HITS_RECON'
    outDirectory = r'C:\Users\mario\OneDrive\Documents\SURE\INTERACTIONS400'

    # create dictionary to get KO1_2
    singleGeneFile = 'allSingleGeneDel_NALM33.txt'
    singleGeneDict = {}
    with open(singleGeneFile, 'r') as file : 
        for lines in file : 
            singleGeneDict[lines.split('\t')[0]] = lines.split('\t')[1]

    # Loop 1 : screens 
    countFile = 0
    for fileName in os.listdir(directory): 
        countFile += 1
        print(countFile)
        filePath = os.path.join(directory, fileName)
        newLines = []
        geneList = []
        with open(filePath, 'r') as file : 
            for lines in file : 
                geneList.append(lines.split('\t')[0])
        
        # Loop 2 : genes in screens 
        for geneName in geneList :
            geneFile = 'HITSdoubleKO6_%(geneName)s_VS_all.txt' %{'geneName' : geneName}
            header = ''
            header += 'Gene_name_1' + '\t' + 'Gene_name_2' + '\t' + 'biomass_flux' + '\t' + 'KO2:KO1_1' + '\t' + 'KO2:KO1_2' + '\t' 
            with open(geneFile, 'r') as file : 
                count = -1
                for lines in file : 
                    count += 1
                    myList_L = []
                    myList_R = []
                    if count == 0 : 
                        sol_no_del = float(lines.split('\t')[1])
                    if count == 2 : 
                        KO1_1 = float(lines.split('\t')[1])
                        newLines.append(header)
                    if count > 4 : 
                        KO2 = float(lines.split('\t')[1])
                        KO1_2 = float(singleGeneDict[lines.split('\t')[0]])

                        # ~ ~ ~ ~ ~ l e t h a l ~ ~ ~ i n t e r a c t i o n s ~ ~ ~ ~ 
                        if abs(KO1_1) > 0.01 : # protection against hit genes that interact with everyone 
                            if abs(KO1_2) == 0.0: # protection against 0 division 
                                KO1_2 = 0.00001
                            minimum = min(KO1_1, KO1_2)
                            if KO2/minimum < cutoff : #THAT WE COULD DO ANOTHER WAY 
                                myList_L.append(geneName) # gene name
                                myList_L.append(lines.split('\t')[0]) # gene name n2
                                myList_L.append(lines.split('\t')[1]) # flux thru biomass rxn
                                myList_L.append(str(KO2/KO1_1)) # KO2 / KO1_1
                                myList_L.append(str(KO2/KO1_2)) # KO2 / KO1_2
                                myList_L.append('lethal')

                        # ~ ~ ~ ~ ~ r e s c u e ~ ~ ~ i n t e r a c t i o n s ~ ~ ~ ~ 
                        if KO2 > 0.1:
                            if KO2 < 0.001 : 
                                KO2 = 0.001
                            if KO1_1 < 0.001 : 
                                KO2 = 0.001
                            if KO1_2 < 0.001 : 
                                KO1_2 = 0.001
                            if KO2/KO1_1 > 1.002 or KO2/KO1_2 > 1.002:
                                myList_R.append(geneName) #gene name no1
                                myList_R.append(lines.split('\t')[0]) #gene name no2
                                myList_R.append(lines.split('\t')[1]) #flux through biomass rxn
                                myList_R.append(str(KO2/KO1_1)) # KO2 / KO1_1
                                myList_R.append(str(KO2/KO1_2)) # KO2 / KO1_2
                                myList_R.append('rescue')

                    if len(myList_L) != 0 : 
                        newLines.append(myList_L)
                    if len(myList_R) != 0 : 
                        newLines.append(myList_R)

        outFileName = fileName[:-15] + '_6_interactions.txt'
        outFilePath = os.path.join(outDirectory, outFileName)
        with open(outFilePath, 'w') as file : 
            for i in range(len(newLines)): 
                for j in range(len(newLines[i])): 
                    file.write(str(newLines[i][j]))
                    file.write('\t')
                file.write('\n')
    return 

interactions(0.97)
