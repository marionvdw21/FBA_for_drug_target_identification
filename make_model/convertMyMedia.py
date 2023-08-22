
def convertMyMedia(myModel, mediaFile): 

    '''
    U S A G E : 
    This function transforms the metabolic model so that the media is more restrictive, hence more realistic. 
    This is done by modifying the bounds of the exchange and import reactions. 
    An example of usage is provided at the end. 

    I N P U T S : 
    --> myModel : metabolic model (of class cobra.core.Model) on which the media modification will be performed. 
    --> mediaFile : a text file which contains the boundary reactions (exchanges, demands and sinks) that will be allowed in the metabolic model. 
        Apart from those reactions, no other boundary reaction will be allowed.
        Each line of the text file has the following format : 'reaction id : reaction description'. An example is provided as 'minimal_media.txt'.
        The important part if the reaction id, so the text file's lines could also be of format 'reaction id : '
    
    O U T P U T : 
    --> model : a cobra.core.Model object, similar to the inputed model exept the boundary reactions have more restrictive upper and lower bound, 
        which mimics a more realistic media.  

    '''

    model = myModel
    # -------------------- limiting the exchange reactions ---------------------------#
    for reaction in model.exchanges: # all the exchange reaction are rendered export reactions only
        reaction.upper_bound = 1000
        reaction.lower_bound = 0
    for reaction in model.demands: # all the import reactions are knocked out 
        reaction.upper_bound = 0
        reaction.lower_bound = 0
    for reaction in model.sinks : 
        reaction.upper_bound = 0
        reaction.lower_bound = 0
    with open('minimal_media2.txt', 'r') as file : # then, only reactions in the mediaFile are allowed to carry flux
        for lines in file : 
            model.reactions.get_by_id(lines.split(':')[0]).upper_bound = 1000
            model.reactions.get_by_id(lines.split(':')[0]).lower_bound = -1000
    #print('minimal media model: ', model.optimize())

    

    # ------------------ knock out odd ways to make o2  -----------------------------# 
    # if the following reactions were allowed to happen, the cell would find a way to make oxygen 
    # itself, which is not realistic and needs to be fixed by knocking out those reactions. 
    model.reactions.RE2112C.upper_bound = 0
    model.reactions.RE2112C.lower_bound = 0
    model.reactions.RE2112R.upper_bound = 0
    model.reactions.RE2112R.lower_bound = 0
    model.reactions.EX_h2o2_e.upper_bound = 0
    model.reactions.EX_h2o2_e.lower_bound = 0

    # ----------------------fixing bounds ------------------------------------------------# 
    # once the model is optimized, the import and exchange reactions that were allowed by the 
    # 'minimal_media' file are fixed. 
    model.optimize()
    with open(mediaFile, 'r') as file : 
        for lines in file : 
            bound = abs(model.reactions.get_by_id(lines.split(':')[0]).flux)
            if bound == 0.0: 
                model.reactions.get_by_id(lines.split(':')[0]).bounds = (-10 , 10)
            else : 
                model.reactions.get_by_id(lines.split(':')[0]).bounds = (-bound, bound)
    model.optimize()

    return model

if __name__ == '__main__': 
    # example of usage : 
    import convertMyModel
    (convertMyMedia(convertMyModel('expression_Recon_ordered2.txt')[0]), 'minimal_media2.txt')