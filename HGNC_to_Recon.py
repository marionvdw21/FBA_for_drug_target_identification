def HGNC_to_Recon_converter(ExpressionFile, myModel): 
	import numpy as np
	import re
	import cobra.io as co
	
				
	model = myModel
	
	with open('table.txt', mode = 'r') as table: 
		HGNC_names_table = []
		EntrezID_list = []
		previous_symbols = []
		alias_symbols = []
		for lines in table : 
			lines_split = lines.split('\t')
			HGNC_names_table.append(lines_split[1])
			alias_symbols.append(lines_split[5])
			previous_symbols.append(lines_split[4])
			EntrezID_list.append(lines_split[6])
	
	with open(ExpressionFile, mode = 'r') as CRISPR_data_new: 
		new_lines = []
		k = 0
		found2 = 0
		unmatched = []
		new_lines = []
		
		
		for lines in CRISPR_data_new : 	
			k += 1
			vache = 0
			EntrezID = ''
			lines_split = lines.split('\t')
			HGNC_name = (lines_split[0])	#change depending on file
			if HGNC_name in HGNC_names_table: 
				index = HGNC_names_table.index(HGNC_name)
				EntrezID = EntrezID_list[index]
				#new_lines.append(lines.replace(HGNC_name, EntrezID))
				
			elif HGNC_name in previous_symbols: 
				index = previous_symbols.index(HGNC_name)
				EntrezID = EntrezID_list[index]
				#new_lines.append(lines.replace(HGNC_name, EntrezID))
				
			elif HGNC_name in alias_symbols: 
				found += 1
				index = alias_symbols.index(HGNC_name)
				EntrezID = EntrezID_list[index]
				#new_lines.append(lines.replace(HGNC_name, EntrezID))
				
			
			else : 
				found = 0
				for i in range(len(previous_symbols)): 
					previous = previous_symbols[i].split()
					for j in range(len(previous)): 
						previous[j] = re.sub(',', '', previous[j])
					if HGNC_name in previous : 
						EntrezID = EntrezID_list[i]
						#new_lines.append(lines.replace(HGNC_name, EntrezID))
						found = 1
						
						break
				if found == 0 : 
					for i in range(len(alias_symbols)): 
						alias = alias_symbols[i].split()
						for j in range(len(alias)): 
							alias[j] = re.sub(',', '', alias[j])
						if HGNC_name in alias: 
							EntrezID = EntrezID_list[i]
							#new_lines.append(lines.replace(HGNC_name, EntrezID))
							found = 1

							break
				if found == 0 and HGNC_name[0:3] == 'LOC' : #all those genes that start with 'LOC'
					EntrezID = HGNC_name[3:]
					vache = 1
					found = 1
				if found == 0 : 
					if HGNC_name == 'FLJ44635': 
						EntrezID = str(340527)
						found = 1
					if HGNC_name == 'HDGFRP2': 
						EntrezID = str(84717)
						found = 1
					if HGNC_name == 'SGK223': 
						EntrezID = str(157285)
						found = 1
					if HGNC_name == 'FLJ45513': 
						EntrezID == str(729220)
						found = 1
					if HGNC_name == 'ZHX1-C8ORF76': 
						EntrezID = str(100533106)
						found = 1
					if HGNC_name == 'SF3B14': 
						EntrezID = str(51639)
						found = 1
					if HGNC_name == 'HDGFRP3': 
						EntrezID = str(50810)
						found = 1
					if HGNC_name == 'HGC6.3': 
						EntrezID = str(100128124)
						found = 1
					
					
				if found == 0 : 
					EntrezID = HGNC_name
					unmatched.append(HGNC_name)
					#new_lines.append(lines)
				
					
	
			BiGGID = EntrezID + '_AT1'
			if BiGGID in model.genes: 
				found2 += 1
				new_lines.append(lines.replace(HGNC_name, BiGGID))
				if vache == 1: 
					print(HGNC_name)

	
				
	#choosing a relevant name for the output file 
	assayName =  ExpressionFile[0:-4]
	OutputFileName = '%(assayName)s_Recon.txt' %{'assayName': assayName} 			
			
	with open(OutputFileName, 'w') as file:
		for genes in range(len(new_lines)):
			file.write(new_lines[genes])  
	
	return open(OutputFileName, 'r')
	

if __name__ == '__main__': 
	HGNC_to_Recon_converter(ExpressionFile, myModel)