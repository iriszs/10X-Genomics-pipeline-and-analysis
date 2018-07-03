#!/usr/bin/python3

import os
import subprocess


#inputdirectory = ""

#if run_id is empty, use the name of the folder where the data is in
#run_id = ""

#will be auto-generated when knowing the input directory. (is output of CellRanger)

#output = ""
#generate by mm10 or GRCh38
#species = ""


#Annotations can be empty.
annotations = ""


inputdirectory = "/media/mdubbelaar/BD_1T/Iris/Data/BE_10X_pilot_II/"
run_id = "10X_pilot_II"
output = "/media/mdubbelaar/BD_1T/Iris/Human_Sample/outs/filtered_gene_bc_matrices/GRCh38/"
species = "Human"

#run_id = "Tung data"
#output = "/media/mdubbelaar/BD_1T/Iris/Tung_data/"
#species = "Human"
#annotations = "annotations.txt"




#user kiezen of ze willen filteren of niet


def perform_Analysis(output, run_id):
	#subprocess.call(["Rscript", "--vanilla", "/media/imgorter/BD_1T/Iris/Scripts/Pipeline/basic_analysis.R", output, run_id, species, annotations])
	subprocess.call(["Rscript", "--vanilla", "/media/mdubbelaar/BD_1T/Iris/Scripts/Pipeline/basic_analysis.R", output, run_id, species, annotations])
	
	


def main():
	perform_Analysis(output, run_id)


if __name__ == "__main__":
	main()	

