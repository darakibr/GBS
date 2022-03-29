# DOWNLOAD the JUNO data set

import os
from datetime import date
import re
import pandas as pd
import numpy as np
import urllib
from urllib import request
from bs4 import BeautifulSoup

cwd = os.getcwd()+'/'
fastqloc = cwd+"fastq/"
if not os.path.exists(fastqloc):
	os.makedirs(fastqloc)

meta = pd.read_csv('JunoMetaData.csv')
linklist = meta.Download_Link

# Downloads fastq.gz files using links
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

# CODE to actually Download the files. Do not run unless original DL files lost. #
for link in linklist:
	getfastq(link)

meta['Study_Name'] = "Juno Project"
