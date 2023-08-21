
def orderExpression(expressionFile): 

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
    orderExpression()
