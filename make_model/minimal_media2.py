import cobra
from convert2NALM6 import convertToNalm6

def minimal_media_fixed(myModel): 
    model = myModel
    # -------------------- limiting the exchange reactions ---------------------------#
    for reaction in model.exchanges: 
        reaction.upper_bound = 1000
        reaction.lower_bound = 0
    for reaction in model.demands: 
        reaction.upper_bound = 0
        reaction.lower_bound = 0
    with open('minimal_media2.txt', 'r') as file : 
        for lines in file : 
            model.reactions.get_by_id(lines.split(':')[0]).upper_bound = 1000
            model.reactions.get_by_id(lines.split(':')[0]).lower_bound = -1000
    print('minimal media model: ', model.optimize())
    sol1 = model.optimize()
    

    # -------------- KO odd ways to make o2 for both models -----------------------------# 
    model.reactions.RE2112C.upper_bound = 0
    model.reactions.RE2112C.lower_bound = 0
    model.reactions.RE2112R.upper_bound = 0
    model.reactions.RE2112R.lower_bound = 0
    model.reactions.EX_h2o2_e.upper_bound = 0
    model.reactions.EX_h2o2_e.lower_bound = 0

    model.optimize()
    with open('minimal_media2.txt', 'r') as file : 
        for lines in file : 
            bound = abs(model.reactions.get_by_id(lines.split(':')[0]).flux)
            if bound == 0.0: 
                model.reactions.get_by_id(lines.split(':')[0]).bounds = (-10 , 10)
            else : 
                model.reactions.get_by_id(lines.split(':')[0]).bounds = (-bound, bound)
    model.optimize()
    print(model.optimize())

    return model

print(minimal_media_fixed(convertToNalm6('expression_Recon_ordered2.txt')[0]))