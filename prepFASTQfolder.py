import os
from datetime import date
import re
import urllib
from urllib import request
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

today = date.today().strftime("%Y%m%d")
cwd = os.getcwd()
def checkcwd(cwd):
    print("Is this the directory you want to work on? "+cwd)
    answercwd = input("Type [y/n]: ")
    if answercwd == "y":
        return cwd
    else:
        newpath = input("New path: ")
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        checkcwd(cwd)

cwd = checkcwd(cwd)

### ENTER FOLDER NAME and location of ILLUMINA or FASTQ source files
name = input("Enter folder name for project: ") #"Juno_Project"
folder = cwd+"/"+name+"/" #set folder path
illumina = folder+"fastq/"

### FOR ILLUMINA
def prepilluminafiles(fastqfolder):
  samples = []
  other_files = []
  counter = 0
  for file in os.listdir(fastqfolder):
		try:
			if file.endswith("R1_001.fastq.gz"):
				file_rn=str(file)
				file_rn = re.sub('_S.+?R','_R', file_rn)
				file_rn = re.sub('_001','',file_rn)
				os.rename(fastqfolder+file,fastqfolder+file_rn)
				samples.append(re.sub('_R1.fastq.gz','',file_rn))
				counter += 1
			elif file.endswith("R2_001.fastq.gz"):
				file_rn=str(file)
				file_rn = re.sub('_S.+?R','_R', file_rn)
				file_rn = re.sub('_001','',file_rn)
				os.rename(fastqfolder+file,fastqfolder+file_rn)
				counter += 1
			elif file.endswith("_R1.fastq.gz"):
				file_rn=str(file)
				samples.append(re.sub('_R1.fastq.gz','',file_rn))
				counter += 1
			elif file.endswith("_R2.fastq.gz"):
				counter += 1
			else:
				other_files.append(str(file))
				counter += 1
		except Exception as e:
			raise e
			print("No files found...")
	print("Number of files renamed = ",counter)
	print("================ FINISHED RENAMING FASTQ TO PROPER INPUT NAMES ================")

prepilluminafiles(illumina)

### FOR DOWNLOAD OF FILES FOR JUNO_PROJECT
listoflinks = pd.read_csv('JunoMetaData.csv').Download_Link
# Function to automatically get fastqfiles from a DL link
def getfastq(link):
	url = str(link)
	soup = BeautifulSoup(request.urlopen(url),features='html.parser')
	templist = soup.find_all(href=re.compile('fastq'))
	isolate = url.split(sep='/')[-1]
	count = 1
	for each in templist:
		download = str(each).split(sep='"')[1]
		templink = url+'/'+download
		name = re.sub('JN_US_','','isolate)
		filename = final_directory+name+"_R"+str(count)+".fastq.gz"
		request.urlretrieve(templink, filename)
		count += 1
for link in linklist:
	getfastq(link)

