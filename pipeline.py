#!/usr/bin/python3

import os
import subprocess


inputdirectory = "/media/imgorter/BD_1T/Iris/Data/BE_10X_pilot_II/"
#if run_id is empty, use the name of the folder where the data is in
#run_id = "10X_pilot_2"
#output = "/media/imgorter/BD_1T/Iris/Human_Sample/outs/filtered_gene_bc_matrices/GRCh38/"
#species = "Human"

run_id = "10X_pilot_1"
output = "/media/imgorter/BD_1T/Iris/Mouse_Sample/outs/filtered_gene_bc_matrices/mm10/"
species = "Mouse"



def perform_Analysis(output, run_id):
	subprocess.call(["Rscript", "--vanilla", "/media/imgorter/BD_1T/Iris/Scripts/Pipeline/basic_analysis.R", output, run_id, species])
	
	


def main():
	perform_Analysis(output, run_id)


if __name__ == "__main__":
	main()	

