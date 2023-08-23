def singleGeneDeletion(myModel):

	'''
	U S A G E : 
	the singleGeneDeletion() function take a given metabolic model and performs a single gene deletion simulation 
	for every gene in the model, and then optimizes this new model to maximize biomass production. 

	I N P U T : 
	--> myModel : the cobra.core.Model for which the single gene deletion needs to be executed. 

	O U T P U T : 
	--> 'allSingleGeneDel.txt' : a text file containing, for each gene in the model, the name of the gene and 
		the resulting biomass production upon deletion of this gene. For each gene, two deletion simulations are made : 
		1. the reactions associated to this gene are knocked out to 10% of their original value 
		2. the reactions associated to this are completely knocked out (upper and lower bounds of (0, 0)). 
		The outputed text file is of format : 
		name of knocked out gene | biomass production when 10% knocked out | time taken to simulate 10% knock out (seconds) | biomass production when total knock out | time taken to simulate total knockout (seconds)

	Note : Performing single gene deletion simulation for all the genes in the model can be very long, especially when 
	dealing with large metabolic models such as Recon model. 
	'''


	#should take about 6h to run for all Recon3 genes. 
	import cobra
	import cobra.flux_analysis 
	import time 
	import cobra.manipulation
	import cobra.flux_analysis.reaction
	import cobra.io as co
	import copy

	model = myModel
	#Making a list (geneList2) of genes of form : [[26_AT1], [314_AT1, 314_AT2], ...]
	geneList = []
	for a in range(len(model.genes)):
		geneList.append(str(model.genes[a]))
	geneList2 = []
	for gene in geneList: 
		if gene[-3:] == 'AT1': 
			allATs = []
			allATs.append(gene)
			geneName = gene[:-3]
			AT = 1
			for i in range(20): 
				AT += 1
				gene2 = '%(geneName)sAT%(AT)s' %{'geneName': geneName, 'AT': AT}
				if gene2 in geneList : 
					allATs.append(gene2)
			geneList2.append(allATs)

	newLines = []
	string = ''
	string += 'WT' + '\t' + str(model.optimize()) #line for original model
	newLines.append(string)
	header = 'genes knocked out' + '\t' + 'BIOMASS_maintenance flux using 10% technique' + '\t' + 'time technique 1'#header
	header += '\t' + 'BIOMASS_maintenance flux using total KO' + '\t' + 'time technique 2'
	newLines.append(header)
	#looping over genes
	for i in range(len(geneList2)): 
		string = ''
		KOgenes = geneList2[i] #get gene ID
	
		#TECHNIQUE 1 : setting fluxe of gene-associated rxn to 10% of WT flux
		start_time1 = time.time()
		newModel = model.copy()    
	
		#looping over reactions > looping over genes 
		for KOgene in KOgenes: 
			genesrxn = newModel.genes.get_by_id(KOgene).reactions
			for reaction in genesrxn : 
				reaction.lower_bound = model.reactions.get_by_id(reaction.id).lower_bound*0.1 #setting the lower bound of the rxn to 10% of that of the WT. 
				reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*0.1 #setting the upper bound of the rxn to 10% of that of the WT. 
		newModel.optimize()
		sol = newModel.reactions.BIOMASS_maintenance.flux
		string += str(KOgenes) + '\t' + str(sol) + '\t' + str(time.time()-start_time1)

		#TECHNIQUE 2 : regular (total) KO
		start_time2 = time.time()
		newModel2 = model.copy()
		for KOgene in KOgenes: 
			genesrxn = newModel2.genes.get_by_id(KOgene).reactions
			for reaction in genesrxn: 
				reaction.lower_bound = 0
				reaction.upper_bound = 0
		newModel2.optimize()
		sol2 = newModel2.reactions.BIOMASS_maintenance.flux
		string += '\t' + str(sol2) + '\t' + str(time.time() - start_time2) 
		newLines.append(string)
	

	#Writing in Output File 
	OutputFileName = 'allSingleGeneDel.txt'  # can change ouputFileName
	with open(OutputFileName, 'w') as file : 
		for lines in range(len(newLines)):
			file.write(str(newLines[lines])) 
			file.write('\n')

	return open(OutputFileName, 'r')

if __name__ == '__main__': 
	# example of usage : 
	import cobra.io as co
	singleGeneDeletion(co.read_sbml_model('Recon3.txt'))
