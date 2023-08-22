def PredictEss(FileName, OutFileName):

    '''
    U S A G E : 
    Function used to order a cell's gene expression file, as a function of the gene's expression.
    An example of usage if provided at the end. 
    I N P U T : 
    --> fileName : a tab delimited text file outputed from the 'single_gene_deletion' script. 
        An example is provided as 'allSingleGeneDel.txt'. 
    --> outFileName : the name given to the output file.

    O U T P U T : 
    --> 'outFileName.txt' : a tab delimited text file ordered from most essential to least essential gene, 
        based on the biomass production resulting from single gene deletion. 
        
    '''
    import time 
    from collections import namedtuple
    from operator import attrgetter
    from bisect import bisect, insort

    def is_float(string): 
        try : 
            float(string)
            return True
        except ValueError: 
            return False

    reducedExp = []
    with open('expression_Recon.txt', 'r') as file: 
                gene_names = []
                expression = []
                for lines in file: 
                     lines_split = lines.split('\t')
                     gene_names.append(lines_split[0])
                     expression.append(lines_split[1])
                for i in range(len(gene_names)): 
                     if is_float(expression[i]) == True : 
                          if float(expression[i]) < 7: 
                               reducedExp.append(str(gene_names))

    ordered = []
    Genes = namedtuple('Gene', ('name', 'flux'))
    with open(FileName) as file: 
        count = 0
        for lines in file: 
            count +=1
            myList = []
            line_split = lines.split('\t') 
            if count > 2 : # to skip the header, CHANGE IT DEPENDING ON SINGLE GENE DELETION FILE 
                if line_split[1] not in reducedExp: 
                    print(line_split[0])
                    print(line_split[1])
                    myList.append(line_split[0]) # CHANGE IT DEPENDING ON SINGLE GENE DELETION FILE 
                    myList.append(line_split[3]) # CHANGE IT DEPENDING ON SINGLE GENE DELETION FILE. Here, total KO is used, but myList.append(line_split[1]) would take the 10% value 
                    myTuple = tuple(myList)
                    myTuple = Genes(*myList)
                    ordered.append(myTuple)


                            
    by_flux = attrgetter('flux')
    ordered.sort(key=by_flux)

    with open(OutFileName, 'w') as outputFile: 
        for i in range(len(ordered)): 
            outputFile.write(ordered[i].name)
            outputFile.write('\t')
            outputFile.write(ordered[i].flux)
            outputFile.write('\n') # CHANGE IT DEPENDING ON FILE

    return 

if __name__ == '__main__': 
    # example of usage : 
    PredictEss('allSingleGeneDel.txt', 'EssPred.txt')



