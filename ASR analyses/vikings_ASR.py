#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~USAGE~~~~~~~~~~~~~~#

		#  #Python environment
		# 	#install vistual environment 
		# 	python3 -m pip install virtualenv

		# 	#create env
		# 	cd ~
		# 	python3 -m venv vikings
		# 	source ./vikings/bin/activate
		# 	cd PATH/TO/VIKING/PROJECT
		# 	pip install --upgrade pip 
		# 	pip install -r requirement_python.txt 


		# #R environment
		# 	# sudo apt install r-base-core
		# 	Rscript requirement_R.txt

		# #If the python environment and the R module are installed run: 
		# 	source ./vikings/bin/activate

		# # Then run the script
		# 	./vikings_ASR.py -t 13 -f /home/vincent/Desktop/vikings/csv/



#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
import time
from datetime import datetime
import os
import shutil
from glob import glob
import pandas as pd
import subprocess
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from joblib import Parallel, delayed
import argparse
import sys




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def make_folder(foldername):
	"""Erase and create folders"""
	if os.path.exists(foldername):
		shutil.rmtree(foldername)
		time.sleep(5)
	os.makedirs(foldername)

def make_nexus_file(file_initial_csv,header_columns_characters, folder_nexus_files, df_transformed_to_analyse_with_R, c):
		""" create nexus file for each character whtin the csv file
		file_initial_csv = single initial csv file
		header_columns_characters = from list_headers, single character
		folder_nexus_files = folder to store the new nexus file
		df_transformed_to_analyse_with_R = df with all variable characters
		c = iterate through the df for character (select column)
		"""
		species_name = file_initial_csv.split(".")[0]
		tested_character = header_columns_characters.replace(" ", "_")
		tested_character = tested_character.replace(",", "_")
		tested_character = tested_character.replace("-", "_")
		tested_character = tested_character.replace("(", "_")
		tested_character = tested_character.replace(")", "_")
		tested_character = tested_character.replace("__", "_")
		tested_character = tested_character.replace("_.", ".")
		file_name_nexus = species_name + "_" + tested_character + ".nex"
		outfile = open(folder_nexus_files + "/" + file_name_nexus, "w")
		outfile.write("#NEXUS\n")
		outfile.write("\tBEGIN TAXA;\n")
		outfile.write("\tTITLE Countries;\n")
		outfile.write("\tDIMENSIONS NTAX=6;\n")
		outfile.write("\tTAXLABELS\n")
		outfile.write("\t\tDanish Faroese Icelandic_ST NorwayN NorwayS Swedish_List")
		outfile.write("\nEND;\n\n")
		nex_df = df_transformed_to_analyse_with_R.columns[c]
		nex_df = nex_df.replace(" ", "_")
		nex_df = nex_df.replace("-", "_")
		nex_df = nex_df.replace(",", "_")
		nex_df = nex_df.replace("__", "_")
		outfile.write("BEGIN CHARACTERS;\n\tTITLE  Character_Matrix;\n\tDIMENSIONS  NCHAR=1;\n\tFORMAT DATATYPE = STANDARD GAP = - MISSING = ? SYMBOLS =  \"0 1\";\n\tCHARSTATELABELS \n\t\t1 " + str(tested_character) + " ; \n\tMATRIX\n")
		df_matrix = df_transformed_to_analyse_with_R.iloc[:, c]
		df_matrix = df_matrix.sort_index()
		df_matrix.rename(index={"DEN" : "Danish", "ICE" : "Icelandic_ST", "FOR" : "Faroese", "NORN":"NorwayN", "NORS":"NorwayS", "SWE":"Swedish_List"}, inplace=True)
		#print(df_matrix)
		outfile.write(df_matrix.to_string())
		outfile.write("\n;\nEND;\n\n")
		outfile.write("BEGIN TREES;\n\tTitle 'Language tree';\n\tLINK Taxa = countries;\n\tTRANSLATE\n")
		outfile.write("\t\t1 Danish,\n")
		outfile.write("\t\t2 Faroese,\n")
		outfile.write("\t\t3 Icelandic_ST,\n")
		outfile.write("\t\t4 NorwayN,\n")
		outfile.write("\t\t5 NorwayS,\n")
		outfile.write("\t\t6 Swedish_List;\n")
		outfile.write("\tTREE pruned_tree_nordic_languages = ((2:354.782664,3:354.782664):581.405169,(((4:168,5:168):168.858007,1:336.858007):385.501458,6:722.359464):967.645246);;\n")
		outfile.write("\nEND;\n\n")
		outfile.close()
		return species_name, tested_character, file_name_nexus

def make_r_sript(species_name,tested_character,folder_R, folder_with_initial_csv, folder_nexus_files, file_name_nexus):
	""" create a file with the R script with aceER and aceARD analyses
	"""
	#log
	r_script_name_log = species_name + "_" + tested_character + "_log.txt"
	r_script_name_log = r_script_name_log.replace("__", "_")
	#r script
	file_name_r_script = species_name + "_" + tested_character + ".r"
	#jpg-aceER
	r_script_name_jpg_aceER = species_name + "_" + tested_character + "_aceER.jpg"
	#jpg-aceAR
	r_script_name_jpg_aceARD = species_name + "_" + tested_character + "_aceARD.jpg"
	outfile_r = open(folder_R + "/" + file_name_r_script, "w")
	#outputs log file
	outfile_r.write("sink(\"" + folder_with_initial_csv + folder_R + "/" + r_script_name_log + "\"" + ", append=FALSE, split=TRUE)\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ Import packages ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("library(\"ape\")\n")
	#outfile_r.write("library(\"phytools\")\n")
	outfile_r.write("library(\"abind\")\n")
	outfile_r.write("library(\"geiger\")\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("character<-\"" + tested_character + "\"\n")
	outfile_r.write("cat(character)\n")
	outfile_r.write("folder_nexus_files=\"" + folder_with_initial_csv + folder_nexus_files + "\"\n")
	outfile_r.write("setwd(folder_nexus_files)\n")
	outfile_r.write("tree = read.nexus(\"" + file_name_nexus + "\")\n")
	outfile_r.write("tree= multi2di(tree)\n")
	outfile_r.write("tree$edge.length <- tree$edge.length / 1000 # rescale tree to keep ML surface happy. Branches are now in 1000’s of years not years.\n")
	outfile_r.write(tested_character + " = read.nexus.data(\"" + file_name_nexus + "\")\n")
	outfile_r.write("states <- as.numeric(" + tested_character + ") # ensure states are in tree-wise order\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ ML RECONSTRUCTION ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ aceER analyse data ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("analysis<-\"aceER_analysis\"\n")
	outfile_r.write("aceER = ace(states, tree, type=\"discrete\", model=\"ER\", use.expm=TRUE, use.eigen=TRUE)\n")
	outfile_r.write("aceER # view results\n")
	outfile_r.write("cat(analysis)\n")
	outfile_r.write("matrix_node <-aceER$lik.anc[,1]\n")
	outfile_r.write("matrix_node\n")
	outfile_r.write("\n#outputs tree - aceER\n")
	outfile_r.write("png(\"" + folder_with_initial_csv + folder_R + "/" + r_script_name_jpg_aceER + "\",height=870,width=574)\n")
	outfile_r.write("plot(tree)\n")
	outfile_r.write("tiplabels(pch=22, bg=as.numeric(" + tested_character + ")+1, cex=2, adj=0.2) \n")
	outfile_r.write("nodelabels(pie=aceER$lik.anc, piecol=c(\"black\", \"red\"))\n")
	outfile_r.write("dev.off ();\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ acARD analyse data ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("analysis<-\"aceARD_analysis\"\n")
	outfile_r.write("aceARD = ace(states, tree, type=\"discrete\", model=\"ARD\", use.expm=TRUE, use.eigen=TRUE)\n")
	outfile_r.write("aceARD\n")
	outfile_r.write("cat(analysis)\n")
	outfile_r.write("aceARD$lik.anc\n")
	outfile_r.write("\n#outputs tree - aceARD\n")
	outfile_r.write("png(\"" + folder_with_initial_csv + folder_R + "/" + r_script_name_jpg_aceARD + "\",height=870,width=574)\n")
	outfile_r.write("plot(tree)\n")
	outfile_r.write("tiplabels(pch=22, bg=as.numeric(" + tested_character + ")+1, cex=2, adj=0.2) \n")
	outfile_r.write("nodelabels(pie=aceARD$lik.anc, piecol=c(\"black\", \"red\"))\n")
	outfile_r.write("dev.off ();\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ Marginal vs Join ancestral states ~~~~~~~~~~~~~~~~#\n")

	outfile_r.write("#plot(aceER$lik.anc,aceER$lik.anc,xlab=\"ace(...,marginal=TRUE)\",ylab=\"ace(...,marginal=FALSE)\")\n")
	outfile_r.write("#plot(aceARD$lik.anc,aceARD$lik.anc,xlab=\"ace(...,marginal=TRUE)\",ylab=\"ace(...,marginal=FALSE)\")\n")

	outfile_r.write("\n#~~~~~~~~~~~~~~~~ Test different evolutionary models (AIC) ~~~~~~~~~~~~~~~~#\n")
	outfile_r.write(tested_character + "_transformed<-abind(lapply(" + tested_character + ", unlist), along=1) #reformat the character matrix for AIC test.\n")

	outfile_r.write("# Create a table to store de results\n")
	outfile_r.write("results.anc <- data.frame(model=c(\"ER\",\"ARD\"), lnL=numeric(2),AICc=numeric(2),params=numeric(2))\n")
	outfile_r.write("# Fit ER\n")
	outfile_r.write("ER_fit <- fitDiscrete(tree, " + tested_character + "_transformed,model=\"ER\")\n")
	outfile_r.write("# Store the results\n")
	outfile_r.write("results.anc[1,-1]<- c(lnL=ER_fit$opt$lnL,AICc=ER_fit$opt$aicc,ER_fit$opt$k)\n")
	outfile_r.write("# Fit ARD\n")
	outfile_r.write("ARD_fit <- fitDiscrete(tree," + tested_character + "_transformed,model=\"ARD\")\n")
	outfile_r.write("# Store the results\n")
	outfile_r.write("results.anc[2,-1]<- c(lnL=ARD_fit$opt$lnL,AICc=ARD_fit$opt$aicc,ARD_fit$opt$k)\n")
	outfile_r.write("# Order the results by AIC\n")
	outfile_r.write("results.anc <- results.anc[order(results.anc$AICc),]\n")
	outfile_r.write("results.anc\n")
	outfile_r.close()
	return file_name_r_script

def run_rscript(folder_R, file_name_r_script):
	"""run the R script as a process"""
	out = subprocess.Popen(["Rscript", folder_R + "/" + file_name_r_script],  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	stdout, stderr = out.communicate()

def make_clustermap(df, folder_with_initial_csv, folder_final_results, ancestral_reconstruction_analyses_type):
	""" make cluster map """

	sns.set(font_scale=0.5)
	plt.subplots(figsize=(50, 10))
	g = sns.heatmap(df, xticklabels=True, yticklabels=True, cmap="OrRd", cbar=True, linewidth=0.5, annot=True)
	g.set_ylim(len(df.index), 0)
	#plt.show()
	g.figure.savefig(folder_with_initial_csv + folder_final_results + "/all_species_heatmap_" + ancestral_reconstruction_analyses_type + ".svg", format='svg', dpi=600)
	print("Heatmap for " + ancestral_reconstruction_analyses_type + " saved")

def return_result_lines_from_the_R_log(ancestral_reconstruction_analyses_type, file_name_results_csv,files_log_r_analyses_list, constant_characters_list, species_name,folder_with_initial_csv, folder_results_per_species, ancestral_state_values_dictionnary, df_constant_character):
		""" select the results lines within the R log script"""
		file_name_results_csv = file_name_results_csv + "_results_" + ancestral_reconstruction_analyses_type + ".csv"
		for logfile in files_log_r_analyses_list:# get the line before the match for each log file
			#print(logfile)
			f = open(logfile)
			columns_header = f.readline()
			with open(logfile) as log:
				for input_line in log:
					if not input_line.startswith(ancestral_reconstruction_analyses_type):
						previous_line = input_line
					if input_line.startswith(ancestral_reconstruction_analyses_type):
						result_line = str(previous_line.split()[1])
						columns_header = columns_header.rstrip()
						dictionnary_new_value = {columns_header: [result_line]}
						ancestral_state_values_dictionnary.update(dictionnary_new_value)
					else:
						continue
		if len(constant_characters_list) != 0: #make the df with or without the constant characters
			df_compiled_results = pd.DataFrame.from_dict(ancestral_state_values_dictionnary)
			df_compiled_results = df_compiled_results.rename(index={0: species_name})
			print_number_analyses = "\t\t\t\t\t" + str(len(df_compiled_results.columns)) + " characters were analysed."
			df_constant_character = df_constant_character.rename(index={"DEN": species_name})
			df_constant_character = df_constant_character.head(1)
			df_compiled_results = pd.concat([df_constant_character, df_compiled_results], axis=1, sort=False)
		else:
			df_compiled_results = pd.DataFrame.from_dict(ancestral_state_values_dictionnary)
			df_compiled_results = df_compiled_results.rename(index={0: species_name})
			print_number_analyses = "\t\t\t\t\t" + str(len(df_compiled_results.columns)) + " characters were analysed."
		df_compiled_results.to_csv(folder_with_initial_csv + folder_results_per_species + "/" + file_name_results_csv)
		return df_compiled_results, print_number_analyses

def extract_nodes_values(files_log_r_analyses_list,file_name_results_csv,folder_with_initial_csv, folder_node_values):
	global df_nodes_values
	df_empty = pd.DataFrame(index=["Root","Icelandic-Faroe", "Swedish", "Danish", "Norway"])
	file_name_results_csv = file_name_results_csv + "_results_node_values_aceER_analysis.csv"
	for logfile in files_log_r_analyses_list:
		list_index = ["Root","Icelandic-Faroe", "Swedish", "Danish", "Norway"]
		f = open(logfile)
		columns_header = f.readline()
		columns_header = columns_header.rstrip()
		with open(logfile) as log:
			for input_line in log:
				if input_line.startswith("aceER_analysis"):
					input_line = input_line.rstrip()
					input_line = input_line.split(" ")
					input_line.pop(0)
					df_nodes_values = pd.DataFrame({columns_header: input_line}, index=list_index)
					df_empty = pd.concat([df_empty, df_nodes_values.reindex(df_empty.index)], axis=1)
				else:
					continue
		df_nodes_values = df_empty
	columns_headers = list(df_nodes_values.columns.values)
	columns_headers.sort()
	df_nodes_values = df_nodes_values[columns_headers]
	df_nodes_values.to_csv(folder_with_initial_csv + folder_node_values + "/" + file_name_results_csv)
	return df_nodes_values

def process_csv(folder_with_initial_csv, folder_with_transformed_csv, folder_node_values, df_clustering_map_aceER, df_clustering_map_aceARD, file_initial_csv):
		""" This function is going through all the steps in order make the Ancestral reconstruction:
			1 transform the orgininal csv file into a df with presence/absence coded with 0 and 1.
			2 create nexus file for each character within the presence/abscence matrix.
			3 create nexus file for each character of the presence/abscence matrix.
			4 run the r script as process.
			5 return df for the aceER and acARD analyses
		"""
		# 1 transform the orgininal csv file into a df with presence/absence coded with 0 and 1.
		ancestral_state_values_dictionnary = {}
		new_df = pd.DataFrame()
		df_compiled_results = pd.DataFrame()
		df = pd.read_csv(file_initial_csv)
		df = df[["AREA", "GENUSE", "SPUSE"]]
		area = ["DEN", "FOR", "ICE", "NORN", "NORS", "SWE"]
		area.sort()
		genuse = df['GENUSE'].unique()
		spuse = df['SPUSE'].unique()
		uses_list = [*genuse, *spuse]
		df_transformed = []
		for countries in area:
			df_countries = df.loc[df['AREA'] == countries]
			df_countries_genuse = df_countries[['AREA', 'GENUSE']]
			df_countries_genuse = df_countries_genuse.rename({"AREA": "AREA", 'GENUSE': "USES"}, axis='columns')
			df_countries_spuse = df_countries[['AREA', 'SPUSE']]
			df_countries_spuse = df_countries_spuse.rename({"AREA": "AREA", 'SPUSE': "USES"}, axis='columns')
			df_Concatenated = df_countries_genuse.append(df_countries_spuse, sort=False)
			columns_headers = df_Concatenated["USES"].values
			columns_headers = list(set(columns_headers))
			columns_number = len(columns_headers)
			dfObj = pd.DataFrame([["1"] * columns_number], columns=columns_headers, index=[countries])
			df_transformed.append(dfObj)
		df_transformed = pd.concat(df_transformed, sort=False)
		df_transformed = df_transformed.fillna(0)
		df_transformed = df_transformed.sort_index(axis=1)
		print_current_analyses = "\nCurently analysing: " + file_initial_csv.split(".")[0]
		#remove constant characters from the AR analyse and put them in the result
		#u = df with the character all the same
		u = pd.DataFrame(df_transformed.nunique())
		u = u[u == 1].dropna()
		constant_characters_list = list(u.index.values)
		if len(constant_characters_list) != 0:
			df_transformed_to_analyse_with_R = df_transformed.drop(columns=constant_characters_list)
			df_constant_character = df_transformed[constant_characters_list]
			print_constant_characters = "\t\t\t\t\t" + str(df_constant_character.shape[1]) + " constant characters found, they will be reported in the final results."
			print_variable_characters = "\t\t\t\t\t" + str(df_transformed_to_analyse_with_R.shape[1]) + " characters will be analysed."
		else:
			df_constant_character = pd.DataFrame()
			df_transformed_to_analyse_with_R = df_transformed.drop(columns=constant_characters_list)
			print_constant_characters = ("\t\t\t\t\tNo constant characters found")
			print_variable_characters = ("\t\t\t\t\t" + str(df_transformed.shape[1]) + " characters will be analysed.")
		#save the presence absence matrix into a csv file
		name_csv = file_initial_csv.split(".")[0]
		df_transformed.rename(index={"DEN":"Danish", "ICE":"Icelandic_ST", "FOR":"Faroese", "NORN":"NorwayN", "NORS":"NorwayS", "SWE":"Swedish_List"}, inplace=True)
		df_transformed.to_csv(folder_with_initial_csv + "/" + folder_with_transformed_csv + "/" + name_csv + ".csv")
		# 2 create nexus file for each character within the presence/abscence matrix.
		c = 0
		list_headers = list(df_transformed_to_analyse_with_R)
		for header_columns_characters in list_headers:
			species_name, tested_character, file_name_nexus = make_nexus_file(file_initial_csv,header_columns_characters,folder_nexus_files, df_transformed_to_analyse_with_R, c)
			c = c + 1
		# 3 create nexus file for each character of the presence/abscence matrix.
			#1 make the rscript
			file_name_r_script = make_r_sript(species_name,tested_character,folder_R, folder_with_initial_csv, folder_nexus_files, file_name_nexus)
		# 4 run the r script as process.
			run_rscript(folder_R, file_name_r_script)
		# 5 return df for the aceER and acARD analyses
		file_name_results_csv = file_initial_csv.split(".")[0]
		file_name_results_csv = file_name_results_csv.replace(" ", "_")
		file_name_results_csv = file_initial_csv.split(".")[0]
		files_log_r_analyses_list = glob(folder_with_initial_csv + folder_R + "/" + file_name_results_csv + "*_log.txt")
		ancestral_reconstruction_analyses_type = "aceER_analysis"
		df_compiled_results_aceER, print_number_analyses = return_result_lines_from_the_R_log(ancestral_reconstruction_analyses_type, file_name_results_csv,files_log_r_analyses_list, constant_characters_list, species_name,folder_with_initial_csv, folder_results_per_species, ancestral_state_values_dictionnary, df_constant_character)
		df_clustering_map_aceER = pd.concat([df_compiled_results_aceER, df_clustering_map_aceER], axis=0, sort=False)
		ancestral_reconstruction_analyses_type = "aceARD_analysis"
		df_compiled_results_aceARD, print_number_analyses = return_result_lines_from_the_R_log(ancestral_reconstruction_analyses_type, file_name_results_csv,files_log_r_analyses_list, constant_characters_list, species_name,folder_with_initial_csv, folder_results_per_species, ancestral_state_values_dictionnary, df_constant_character)
		df_clustering_map_aceARD = pd.concat([df_compiled_results_aceARD, df_clustering_map_aceARD], axis=0, sort=False)
		extract_nodes_values(files_log_r_analyses_list, file_name_results_csv, folder_with_initial_csv, folder_node_values)
		print(print_current_analyses)
		print(print_constant_characters)
		print(print_variable_characters)
		print(print_number_analyses)
		return df_clustering_map_aceER, df_clustering_map_aceARD

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Script starting here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
folder_with_initial_csv = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--folder", type=str, help="Folder with the csv files", default=folder_with_initial_csv, required=True)
parser.add_argument("-t", "--thread", type=int, help="Input Number of threads", default=2)
args = parser.parse_args()
folder_with_initial_csv = args.folder
core_number = int(args.thread)
if not os.path.isdir(folder_with_initial_csv):
    print('The path specified does not exist.')
    sys.exit()
else:
	os.chdir(folder_with_initial_csv)

#~~~~~~~~~~~~~ print the start time and parameters ~~~~~~~~~~~~~#
now = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
print("Starting analyse: " + now + '\n')
print("\nThe path to the CSV files is the following: " + folder_with_initial_csv)
print(str(core_number) + " threads will be used. \n")

#~~~~~~~~~~~~~ create the folders to store the different results ~~~~~~~~~~~~~#
folder_with_transformed_csv = "1_character_file_transformed_csv"
make_folder(folder_with_transformed_csv)
# make the folder for nexus files
folder_nexus_files = "2_nexus_files"
make_folder(folder_nexus_files)
# make the folder for R analyses
folder_R = "3_R_analyses"
make_folder(folder_R)
# make the folder for compiled results
folder_results_per_species = "4_compiled_results_per_species"
make_folder(folder_results_per_species)
# make the folder for compiled results
folder_final_results = "5_final_results"
make_folder(folder_final_results)
# make the folder for brachlength values
folder_node_values = "6_node_values"
make_folder(folder_node_values)

#~~~~~~~~~~~~~ gather the list of csv files + empty variables ~~~~~~~~~~~~~#
file_initial_csv_list = glob('*.csv')
file_initial_csv_list.sort()
print(str(len(file_initial_csv_list)) + " csv files will be analyzed.")
for elem in file_initial_csv_list:
		print("\t" + elem)
df_clustering_map_aceER = pd.DataFrame()# define empty df for function return_result_lines_from_the_R_log()
df_clustering_map_aceARD = pd.DataFrame()
df_aceER = pd.DataFrame()
df_aceARD = pd.DataFrame()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ run the script in parallel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
res = Parallel(n_jobs=core_number, prefer="threads")(delayed(process_csv)(folder_with_initial_csv, folder_with_transformed_csv, folder_node_values, df_clustering_map_aceER, df_clustering_map_aceARD, file_initial_csv) for file_initial_csv in file_initial_csv_list)
df_clustering_map_aceER = [item[0] for item in res]
df_clustering_map_aceARD = [item[1] for item in res]
df_aceER = df_aceER.append(df_clustering_map_aceER, ignore_index=False)
df_aceARD = df_aceARD.append(df_clustering_map_aceARD, ignore_index=False)

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAKE THE GRAPHIC OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ancestral_reconstruction_analyses_type = "aceER_analysis"
df_clustering_map = df_aceER.fillna(0)
df_clustering_map = df_clustering_map.sort_index(axis=1)
df_clustering_map = df_clustering_map.sort_index(axis=0)
df_clustering_map.to_csv(folder_with_initial_csv + folder_final_results + "/all_species_aceER.csv")
df = pd.read_csv(folder_with_initial_csv + folder_final_results + "/all_species_aceER.csv", index_col=0)
columns_headers = list(df.columns.values)
columns_headers_general_uses = ['M', 'AF', 'IC', 'F', 'V', 'SSR', 'Fu', 'A', 'Co']
columns_headers_other_uses = list(set(columns_headers) - set(columns_headers_general_uses))
columns_headers_other_uses.sort()
column_order = columns_headers_general_uses + columns_headers_other_uses
df = df[column_order]
df.to_csv(folder_with_initial_csv + folder_final_results + "/all_species_ordered_aceER.csv")
#make_clustermap(df, folder_with_initial_csv, folder_final_results, ancestral_reconstruction_analyses_type)

ancestral_reconstruction_analyses_type = "aceARD_analysis"
df_clustering_map = df_aceARD.fillna(0)
df_clustering_map = df_clustering_map.sort_index(axis=1)
df_clustering_map = df_clustering_map.sort_index(axis=0)
df_clustering_map.to_csv(folder_with_initial_csv + folder_final_results + "/all_species_aceARD.csv")
df = pd.read_csv(folder_with_initial_csv + folder_final_results + "/all_species_aceARD.csv", index_col=0)
columns_headers = list(df.columns.values)
columns_headers_general_uses = ['M', 'AF', 'IC', 'F', 'V', 'SSR', 'Fu', 'A', 'Co']
columns_headers_other_uses = list(set(columns_headers) - set(columns_headers_general_uses))
columns_headers_other_uses.sort()
column_order = columns_headers_general_uses + columns_headers_other_uses
df = df[column_order]
df.to_csv(folder_with_initial_csv + folder_final_results + "/all_species_ordered_aceARD.csv")
make_clustermap(df, folder_with_initial_csv, folder_final_results, ancestral_reconstruction_analyses_type)

#Make heat map for node values

file_node_list = glob(folder_with_initial_csv + folder_node_values + "/*.csv")
for file_csv in file_node_list:
	file_csv_df=file_csv
	file_csv = file_csv.split("/")[-1]
	file_csv = file_csv.split(".")[0]
	file_name_results_svg = file_csv + ".svg"
	df_nodes_values = pd.read_csv(file_csv_df, index_col=0)
	make_clustermap(df_nodes_values, folder_with_initial_csv, folder_node_values, file_name_results_svg)

print("\nanalyses done!")
print(datetime.now())
