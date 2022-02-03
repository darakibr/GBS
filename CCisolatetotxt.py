# Extracting the isolate numbers by clonal complex from the MLST and Clonal Complex Excel sheet

import os
import pandas as pd
import re
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

os.chdir("/Users/daraki/Documents/py_codes/")
cwd = os.getcwd()
#print(cwd)
from datetime import date
today = date.today().strftime("%Y%m%d")
random.seed(1)

def randcolorlist(n,type=1):
    if n < 0 | type < 0:
        raise ValueError('First argument n must be an integer above 0')
    elif n>3 and n<11:
        if type == 1:
            clist = mpl.cm.tab10(np.linspace(0,1,n))
        elif type%2==0:
            clist = mpl.cm.Set1(np.linspace(0,1,n))
        else:
            clist = mpl.cm.tab20b(np.linspace(0,1,n))
    else:
        if type==1:
            clist = mpl.cm.nipy_spectral(np.linspace(0, 1, n))
        else:
            clist = mpl.cm.gist_rainbow(np.linspace(0, 1, n))
    return clist


ccdf = pd.read_excel('MLSTandCC.xlsx')
ccdf = ccdf.rename(columns={"Clonal Complex (6 of 7 matches)":"CC"})
def clonalgroup(cc):
    unique = ccdf['CC'].unique()
    if cc not in unique:
        raise ValueError('There is no Clonal Complex ' + cc + ' in the dataset.')
    for c in unique:
        group = ccdf[ccdf['CC'] == c]
        group = group[['Isolate']].transpose().to_string(header=False,index=False)
        if c == cc:
            with open(str(c)+'clonalgroup.txt', 'w+') as file:
                file.write(group)
        else:
            print("Clonal Complex "+cc+" not found.")
#clonalgroup(459)
#clonalgroup(452)
#clonalgroup(196)

amr = pd.read_excel('GBSIVoutput.xlsx',sheet_name='AMRall')
amr = amr.loc[:,('gene','tetM','tet(W/N/W)','tetO','Erm(A)','ErmB','ErmC','ErmT')]
amr.replace(0,-1,inplace=True)
amr = amr.rename(columns={"gene":"Isolate"})

header = "DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,ABCD\nCOLOR,#ff0000\n\
SHOW_STRIP_LABELS,1\nSTRIP_LABEL_SIZE_FACTOR,0.5\nSTRIP_LABEL_COLOR,#ffffff\n\
COLOR_BRANCHES,TRUEFALSE\nSTRIP_LABEL_SHIFT,-11\nSHOW_LABELS,1\SIZE_FACTOR,1\nLEGEND_DATA\n"
amrheader = "DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,AMR\nCOLOR,#ff0000\n\
LEGEND_TITLE,AMR gene types\nLEGEND_POSITION_X,100\nLEGEND_POSITION_Y,250\nLEGEND_SHAPES,3,3,3,2,2,2,2\n\
LEGEND_COLORS,#FF0000,#0000FF,#FFEA00,#FF0000,#0000FF,#4CBB17,#FFEA00\n\
LEGEND_LABELS,tetM,tet(W/N/W),tetO,Erm(A),ErmB,ErmC,ErmT\nFIELD_SHAPES,3,3,3,2,2,2,2\n\
FIELD_COLORS,#FF0000,#0000FF,#FFEA00,#FF0000,#0000FF,#4CBB17,#FFEA00\n\
FIELD_LABELS,tetM,tet(W/N/W),tetO,Erm(A),ErmB,ErmC,ErmT\nDATA\n"
leafheader = "DATASET_SYMBOL\nSEPARATOR COMMA\DATASET_LABEL,Group\nCOLOR,#ff0000\nMAXIMUM_SIZE,10\n\
LEGEND_TITLE,Group by hierbaps\nLEGEND_POSITION_X,100\nLEGEND_POSITION_Y,250\nLEGEND_SHAPES,3,3,3,2,2,2,2\n\
LEGEND_COLORS,_COLOURS_\n\
LEGEND_LABELS,_NUMBER_\nFIELD_SHAPES,3,3,3,2,2,2,2\n\
FIELD_COLORS,_COLOURS_\n\
FIELD_LABELS,_NUMBER_\nDATA\n"

allcdc = pd.read_csv(cwd+'/CDC/allcdc.csv')
temp1 = allcdc[allcdc['Serotype']=="IV"]
temp1 = temp1.loc[:,('Run','Year')]
temp1 = temp1.rename(columns={"Run":"Isolate"})
bcjb = pd.read_excel('CJBGBSIVmeta.xlsx')
temp2 = bcjb.astype(str)
temp2 = temp2.loc[:,('IsolateName','Year')]
temp2 = temp2.rename(columns={"IsolateName":"Isolate"})
combined = temp1.append(temp2,ignore_index=True, sort=False)
Toronto = pd.read_csv('Toronto.csv')
combined = combined.append(Toronto,ignore_index=True, sort=False)
found = pd.read_csv('Yearsfound.csv')
combined = combined.append(found,ignore_index=True, sort=False)
meta = ccdf.loc[:,('Isolate','ST','CC','Source')]
meta = pd.merge(meta,combined,how='outer')
metayeargroup = meta.fillna(0)
metayeargroup['Year'] = metayeargroup.loc[:,'Year'].astype(int)
cat = np.array(['Pre 2000','2000-2005','2006-2010','2011-2015','2016+'])
metayeargroup.loc[metayeargroup['Year'] >1, 'Yeargroup'] = cat[0]
metayeargroup.loc[metayeargroup['Year'] >2000, 'Yeargroup'] = cat[1]
metayeargroup.loc[metayeargroup['Year'] >2010, 'Yeargroup'] = cat[2]
metayeargroup.loc[metayeargroup['Year'] >2010, 'Yeargroup'] = cat[3]
metayeargroup.loc[metayeargroup['Year'] >2010, 'Yeargroup'] = cat[4]

meta5year = meta.fillna(0)
meta5year['Year'] = meta5year.loc[:,'Year'].astype(int)
cat1 = np.array(['2000-2005','2006-2010','2011-2015','2016+'])
meta5year.loc[meta5year['Year'] >1999, 'Yeargroup'] = cat1[0]
meta5year.loc[meta5year['Year'] >2005, 'Yeargroup'] = cat1[1]
meta5year.loc[meta5year['Year'] >2010, 'Yeargroup'] = cat1[2]
meta5year.loc[meta5year['Year'] >2015, 'Yeargroup'] = cat1[3]

metayeargroup = pd.merge(metayeargroup,amr,how='outer')
meta5year = pd.merge(meta5year,amr,how='outer')

### Adding data from hierbaps, only created for specific 452 and 459 groups ###
baps452 = pd.read_csv(cwd+'/hierbaps_GBS-IV-452.csv')
baps459 = pd.read_csv(cwd+'/hierbaps_GBS-IV-459.csv')
baps452 = baps452.rename(columns={"level 1":"Group"}).drop('level 2', axis=1)
baps459 = baps459.rename(columns={"level 1":"Group"}).drop('level 2', axis=1)
baps459196 = pd.read_csv(cwd+'/hierbaps_GBS-IV-459196.csv')
baps459196 = baps459196.rename(columns={"level 1":"Group"}).drop('level 2', axis=1)
baps196 = pd.read_csv(cwd+'/hierbaps_GBS-IV-196.csv')
baps196 = baps196.rename(columns={"level 1":"Group"}).drop('level 2', axis=1)

### Create meta data file for phylogenetic tree ###
# currently creates Year, Sequence type, hierbaps group, year category and AMR resistance
### HIERBAPS - Way to cluster strings into groups by statistics of shared SNPs.
#       Hierarichal structure of tree (R scrip already written in StrepLab-UBUNTU).

#### Why have the 452 and 459 emerged as the predominant strains
#### Mobile genetic element vs Resistance
#### Identify representative strains that are from 459, 452
#           two clusters, do these differer in ability to cause disease.
#           A few that have same Mobile Genetic Element.
#           Wild type in known regulators (Do not have mutations on CovRS).
#           Filter VCF file to find mutations at specific regions.

# This is used in newdataset()
def newdata(df,column,l=False,c=1):
    try:
        list = df[column].unique()
        colors = randcolorlist(len(list),type=c)
        count = 0
        shapes = ""
        labels = ','.join([str(item) for item in list])
        colornames = ""
        for item in list:
            c = mpl.colors.to_hex(colors[count],keep_alpha=False)
            df.loc[:,'Style'] = np.where(df.loc[:,(column)].str.contains(item),c,df.loc[:,'Style'])
            shapes = shapes+"1,"
            colornames = colornames+str(c)+","
            count +=1
        temp = df['Isolate']+','+df['Style']+','+df[column]
        if l ==True:
            shapes = shapes.rstrip(",")
            legendtext = 'LEGEND_TITLE,'+df[column].name+'\n\
            LEGEND_SHAPES,'+shapes+'\
            \nLEGEND_COLORS,'+colornames[:-1]+'\
            \nLEGEND_LABELS,'+labels+'\n'
            with open('TempLEGEND.txt','w') as file2:
                file2.write(legendtext)
        return temp
    except:
        print("Check that you assigned valid arguments: df=dataframe, column='col_name', l=True/False, c=int()")

# This is used to create .txt file for ITOL tree visualization
def newdataset(df,column,namevar:str,colorbranch=False,legend=False,colorscheme=1):
    try:
        data = newdata(df,column,l=legend,c=colorscheme)
        data.to_csv('temp.csv',index=False, header=False)
        with open('temp.csv', 'r') as file:
            lines = file.read()
            lines = re.sub('"','',lines)
            tempheader = re.sub('ABCD',str(namevar),header)
        if colorbranch == True:
            tempheader = re.sub('TRUEFALSE',str(1),tempheader)
        else:
            tempheader = re.sub('TRUEFALSE',str(0),tempheader)
        if legend == True:
            with open('TempLEGEND.txt', 'r') as file:
                temptext = file.read()
                tempheader = re.sub('LEGEND_',temptext,tempheader)
        else:
            tempheader = re.sub('LEGEND_','',tempheader)
        with open(str(namevar)+'.txt','w') as file2:
            file2.write(tempheader)
            file2.write(lines)
    except:
        print("Check that you assigned valid arguments: df=dataframe, column='col_name', l=True/False, c=int()")


# This creates the Antimicrobial layer .txt file for ITOL tree visualization
def amrset(df, namevar:str):
    try:
        temp = df[['Isolate','tetM','tet(W/N/W)','tetO','Erm(A)','ErmB','ErmC','ErmT']]
        temp.to_csv('temp.csv',index=False, header=False)
        with open('temp.csv', 'r') as file:
            lines = file.read()
            lines = re.sub('"','',lines)
            with open(namevar+'.txt','w') as file2:
                file2.write(amrheader)
                file2.write(lines)
    except:
        print("Check if df = dataset with AMR data, and namevar = str(output_file_name)")

def GroupST(df,n=10):
    try:
        #print('Original Sequence Types:',df.ST.unique())
        df.ST = df.ST.fillna(0)
        for item in df.ST.unique():
            if df.loc[df['ST']==item,'ST'].count() > n:
                df.loc[df['ST']==item,'ST'] = item
            else:
                df.loc[df['ST']==item,'ST'] = 0
        #print("Grouped all ST with less than",n,"occurances into 'other':",df.ST.unique())
        return df
    except:
        print("Check df=dataframe with meta data including sequence type column named 'ST' and n=int() to only select ST with more than n occurances.")

# FOR CC452 tree
temp452 = meta5year[meta5year['CC']==452]
temp452 = pd.merge(temp452,baps452,how='outer')
csvdata1 = pd.DataFrame(data=temp452,columns=['Isolate','ST','CC','Source','Year','Style','Yeargroup','Group'], dtype=str)

temp452.Name = 'CC 452'
newdataset(csvdata1,'Yeargroup',"452YearGroup",colorbranch=False,legend=True,colorscheme=2)
newdataset(csvdata1,'Group',"452Group",colorbranch=True,legend=True)
newdataset(csvdata1,'ST',"452SeqType")
amrset(temp452,"452amr")

temp459 = meta5year[meta5year['CC']==459]
temp459196 = metayeargroup[metayeargroup['ST']==196]
temp459196 = temp459196.append(metayeargroup[metayeargroup['CC']==459])

# FOR CC459 tree
temp459 = pd.merge(temp459,baps459,how='outer')
temp459.Name = 'CC 459'
csvdata2 = pd.DataFrame(data=temp459,columns=['Isolate','ST','CC','Source','Year','Style','Yeargroup','Group'], dtype=str)
newdataset(csvdata2,'Yeargroup',"459YearGroup",colorbranch=False,legend=True)
newdataset(csvdata2,'Group',"459Group",colorbranch=True,legend=True)
newdataset(csvdata2,'ST',"459SeqType",colorscheme=2)
amrset(temp459,"459amr")
CC459groupedST = GroupST(csvdata2)
newdataset(CC459groupedST,'ST','459GroupST',colorscheme=2,legend=True)

temp196 = metayeargroup[metayeargroup['CC']==196]
temp196.Name = 'CC 196'

# FOR CC196 tree
temp196 = pd.merge(temp196,baps196,how='outer')
csvdata4 = pd.DataFrame(data=temp196,columns=['Isolate','ST','CC','Source','Year','Style','Yeargroup','Group'], dtype=str)
newdataset(csvdata4,'Yeargroup',"196Yeargroup",colorbranch=False,legend=True)
newdataset(csvdata4,'Group',"196Group",colorbranch=True,legend=True)
newdataset(csvdata4,'ST',"196SeqType",colorscheme=2)
amrset(temp196,"196amr")

# FOR combined CC459+ST196 tree
temp459196 = pd.merge(temp459196,baps459196,how='outer')
temp459196.Name = 'CC 459 + ST-196'
csvdata5 = pd.DataFrame(data=temp459196,columns=['Isolate','ST','CC','Source','Style','Yeargroup','Group'], dtype=str)
newdataset(csvdata5,'Yeargroup',"459196YearGroup",colorbranch=False,legend=True)
newdataset(csvdata5,'Group',"459196Group",colorbranch=True,legend=True)
newdataset(csvdata5,'ST',"459196SeqType",colorscheme=1,legend=True)
amrset(temp459196,"459196amr")
CC459196groupedST = GroupST(csvdata5)
newdataset(CC459196groupedST,'ST','459196GroupST',colorscheme=3,legend=True)

temp0 = metayeargroup.loc[metayeargroup['Year']==0, ['Isolate','Year']]
#print(metayeargroup.loc[metayeargroup['Year']==0, ['Isolate','ST']])
temp0.to_csv('NoYear'+today+'.csv')

width = 0.4

def getbarheights(df,yrlist,col='Year'):
    try:
        npercat = yrlist.copy()
        counter = 0
        for i in yrlist:
            npercat[counter] = (df.loc[df['Yeargroup']==i,col].count())
            counter +=1
        barheight = list(npercat.astype(int))
        return barheight
    except:
        print("df must be a valid dataframe including column 'Yeargroup' with same number of categories as yrlist argumen, and col='col_name' in df.")

def barbycol(df,column,yrl):
    try:
        uniqueslist = df[column].unique()
        uniqueslist.sort()
        listofheights = list()
        listoftotals = getbarheights(df,yrlist=yrl,col=column)
        for uitem in uniqueslist:
            tempdf = df[df[column]==uitem]
            tempy = getbarheights(tempdf,yrlist=yrl)
            percent = list()
            for b,m in zip(tempy, listoftotals):
                if m>0:
                    percent.append(b/m)
                else:
                    percent.append(0)
            listofheights.append(percent)
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.85,0.85])
        colorlist = randcolorlist(len(listofheights),2)
        for x in range(len(uniqueslist)):
            if x==0:
                ax.bar(yrl, listofheights[x], width, color=colorlist[x])
                h=listofheights[x]
            else:
                ax.bar(yrl, listofheights[x], width,bottom=h, color=colorlist[x])
                h = [sum(x) for x in zip(h, listofheights[x])]
        ax.set_xlabel('Year range')
        ax.set_ylabel('Percent')
        ax.set_title(str(df.Name)+' by '+column+' by year')
        ax.legend(labels=uniqueslist)
        plt.savefig(str(df.Name)+'by'+column+'.png', format='png')
        #plt.show()
    except:
        print("df must be a valid dataframe including column='col_name' with same number of categories in yrl list argument as 'Yeargroup' col in df")

meta5year.Name = 'GBS-IV'

barbycol(temp452,'Group',cat1)
barbycol(temp459,'Group',cat1)
barbycol(meta5year,'CC',cat1)
barbycol(temp459196,'Group',cat)
barbycol(temp196,'Group',cat)

os.remove('temp.csv')
os.remove('TempLEGEND.txt')