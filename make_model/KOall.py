def KOALL():
	#should take about 6h to run for all Recon3 genes. 
	import cobra
	import cobra.flux_analysis 
	import time 
	import cobra.manipulation
	import cobra.flux_analysis.reaction
	import cobra.io as co
	import copy

	#human = r"C:\Users\mario\OneDrive\Documents\SURE\Recon2.xml"
	human = 'Recon3.xml'
	cobra.io.validate_sbml_model(human)
	print('model validated')
	model = co.read_sbml_model(human)
	
	print('readingModelOk')


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


	totTotTime = 0
	newLines = []
	string = ''
	string += 'WT' + '\t' + str(model.optimize()) #line for original model
	newLines.append(string)
	header = 'genes knocked out' + '\t' + 'BIOMASS_maintenance flux using 10% technique' + '\t' + 'time technique 1'#header
	header += '\t' + 'BIOMASS_maintenance flux using total KO' + '\t' + 'time technique 2'
	newLines.append(header)
	#looping over genes
	for i in range(10): #replace w/ len(geneList2)
		print(i)
		string = ''
		KOgenes = geneList2[i] #get gene ID
	
		#TECHNIQUE 1 : setting fluxe of gene-associated rxn to 10% of WT flux
		start_time1 = time.time()
		newModel = model.copy()    
	
		#looping over reactions > looping over genes 
		for KOgene in KOgenes: 
			#genesrxn = newModel.genes.get_by_id(KOgene.id).reactions 
			genesrxn = newModel.genes.get_by_id(KOgene).reactions
			for reaction in genesrxn : 
				reaction.lower_bound = model.reactions.get_by_id(reaction.id).lower_bound*0.1
				reaction.upper_bound = model.reactions.get_by_id(reaction.id).upper_bound*0.1 #setting the upper bound of the rxn to 10% of that of the WT. 
		newModel.optimize()
		sol = newModel.reactions.BIOMASS_maintenance.flux
		string += str(KOgenes) + '\t' + str(sol) + '\t' + str(time.time()-start_time1)

		#TECHNIQUE 2 : regular KO
		start_time2 = time.time()
		newModel2 = model.copy()
		for KOgene in KOgenes: 
			genesrxn = newModel2.genes.get_by_id(KOgene).reactions
			for reaction in genesrxn: 
				reaction.lower_bound = 0
				reaction.upper_bound = 0
			#cobra.manipulation.knock_out_model_genes(newModel2, KOgene)  #i'm not sure if the model is iterated after each KOgene (does it keep the modifications?)
		newModel2.optimize()
		sol2 = newModel2.reactions.BIOMASS_maintenance.flux
		string += '\t' + str(sol2) + '\t' + str(time.time() - start_time2) 
		newLines.append(string)
	
		#totTime = time.time() - start_time1
		#totTotTime += totTime
		
	print(totTotTime)
	#Writing in Output File 
	OutputFileName = r'C:\Users\mario\OneDrive\Documents\SURE\allSingleGeneDelRecon2.txt'   
	with open(OutputFileName, 'w') as file : 
		for lines in range(len(newLines)):
			file.write(str(newLines[lines])) 
			file.write('\n')

	return open(OutputFileName, 'r')

KOALL()