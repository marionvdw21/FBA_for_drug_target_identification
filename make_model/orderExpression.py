
def orderExpression(expressionFile): 

    '''
    U S A G E : 
    This script is used to order the expression of a cell line from most expressed genes to least expressed and vice versa. 
    This is usefull because when performing modification of flux boundaries in a given model, the order 
    in which is is done matters as some genes share the same reactions. 
    The priority in terms of modifying flux boundaries can therefore be given to least expressed or more expressed genes, 
    depending on the expression file used. 
    Example of usage provided at the end. 

    I N P U T S : 
    --> expressionFile : a tab-delimited text file containing the gene names (in Bigg ID fomat) as the first colums, 
    and expression as the second colomn (in log2 read counts). This file can be ouputed from the 'HGNC_to_Recon' script. 

    O U T P U T S : 
    --> expressionFile_ordered.txt : a txt file which is a copy of the input file exept it is ordered from least expressed gene to most expressed gene
    --> expressionFile_deordered.txt : a txt file which is a copy of the input file exept it is ordered from most expressed gene to least expressed gene
    '''

    from collections import namedtuple
    from operator import attrgetter

    fileName = expressionFile[:-4] + '_Recon.txt'

    Genes = namedtuple('Gene', ('name', 'expression'))
    ordered = []
    

    with open(fileName, 'r') as file : 
        for lines in file : 
            myList = []
            myList.append(lines.split('\t')[0])
            myList.append(lines.split('\t')[1])
            myTuple = Genes(*myList)
            ordered.append(myTuple)

    by_expression = attrgetter('expression')
    ordered.sort(key=by_expression, reverse=True)


    with open('expression_Recon_deordered.txt', 'w') as file : 
        for i in range(len(ordered)): 
            for j in range(len(ordered[i])): 
                file.write(ordered[i][j])
                file.write('\t')
    ordered.sort(key=by_expression, reverse=False)

    with open('expression_Recon_ordered.txt', 'w') as file : 
        for i in range(len(ordered)): 
            for j in range(len(ordered[i])): 
                file.write(ordered[i][j])
                file.write('\t')

    return 


if __name__ == '__main__': 
    # example of usage : 
    orderExpression(ExpressionFile = 'expression_Recon.txt')
