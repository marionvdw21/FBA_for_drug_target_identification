

def PredictEss(FileName, OutFileName):
    import time 
    from collections import namedtuple
    from operator import attrgetter
    from bisect import bisect, insort

    '''
    ~~~ purpose ~~~
        Ordering a file containing fluxes as a result of single gene deletion, in order to predict gene essentiality 
        based on flux thru 'BIOMASS_maintenance' reaction. 

    ~~~ objects ~~~
    myList 
        myList - reducedExp. 
        list of genes that are going to be ranked for essentiality (i.e. NALM6 genes w\ those that got
        turned down).
    reducedExp 
        list of gene names whose expression got turned down/knocked out for simulating NALM6 model. 

    '''

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
    PredictEss('allSingleGeneDel.txt', 'EssPred.txt')



