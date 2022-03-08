### EXAMPLE
### Codes used to clean the data produced by snpEff
### Actual files used were modified to fit certain aspects (index and header positions)

### FILE.vcf needs to be cleaned of extra data added by snpEff and transposed to a .csv file. This new.csv file will be called as the starting df.
df = pd.read_csv( '<path_to_file>' , sep='\t', index_col=0, header=17) # Read .csv data made from .vcf file using row and col names as indexes.
df = df.replace('.', np.nan) # Programs use '.' placeholder for missing values

### IF there is additional META data in another file, load that as well, and use pd.merge to join the data together
  #Example
  meta = pd.read_csv( '<path_to_meta_file>' )
  dfmerged =  pd.concat([meta, df], join='outer', axis=1, verify_integrity=True)

count = dfmerged.append(dfmerged.count().rename('Total')) # Adds a row with the total number of isolates per category (snp)
n = len(count)-2 # Total number of isolates, subtract 2 because first 2 rows contain meta data associated with columns
percent = count.append([count.iloc[n,:]/n*100].rename('Percent'))
percent.to_csv(' <path_to_METAwTotals.csv> ' )
  
### EXAMPLE 2, working without transposing first
### pd.join/merge/concat will not work with dataframes with duplicated entries (descriptive header rows)
### Keep data like this for working with single dataframe

df = pd.read_csv( '<path_to_file>' , sep='\t', index_col=0, header=0)

def total(data, add_percent = True, descriptiverows = 23): # There are 23 descriptive rows before individual isolates are listed, other vcf files may vary.
  temp = data.transpose().replace(0,np.nan) # Ensure all cell not to be counted are empty
  data['Total'] = temp.count() # counts not empty values so will count 0s
  if add_percent == True:
    n = len(data)-descriptiverows
    data['Percent'] = data['Total']/n*100
  return data.transpose()

def hotencoder(col,data): # hot encode categorical data for machine learning.
  from sklearn.preprocessing import OneHotEncoder
  encoder = OneHotEncoder()
  enc = encoder.fit_transform(data[[col]])
  catlist = encoder.get_feature_names_out()
  print("Categories encoded:", catlist)
  df= pd.DataFrame(enc.toarray(), columns=catlist)
  return df

def rfcaller(data, features=list(), outcome=str(), n_estimators= 500, perc_test = 0.8, accuracy= True, show_features=10):
  from sklearn.ensemble import RandomForestClassifier
  from sklearn.model_selection import train_test_split
  x = data[features]
  y = data[[outcome]]
  xtrain, xtest, ytrain, ytest = train_test_split(x,y,test_size= perc_test)
  rfc = RandomForestClassifier(n_estimators = n_estimators)
  rfc.fit(xtrain,ytrain)
  if accuracy == True:
    from sklearn.metrics import accuracy_score
    pred = rfc.predict(xtest)
    print("Accuracy:", accuracy_score(ytest,pred))
  if show_features >0:
    features = pd.Series(rfc.feature_importances_, index=features).sort_values(ascending=False)
    print(features[0:show_features])
  print("Random Forest Model done as rfc.")

def addgenecountcol(df, countsdf, ct):
  n = list(df.index).index('GeneName') # get index location of GeneName column in original dataframe
  for colname in countsdf.column.unique():
    templist = list(countsdf[colname])
    templist.insert(n,colname)
    df[ct+colname] = templist
  return df

def mutationsbyisolate(df, counttype=('all','stop')):
  if counttype = 'all':
    snpcol = [col for col in df.columns if 'p.' in col] # list of all columns that contain if a mutation is present in isolate (1= present, np.nan= not found).
    snpcoldf = df[snpcoldf]
    tdf = snpcoldf.transpose().replace(0,np.nan).groupby('GeneName').count()
    allcounts = tdf.transpose()
    addgenecountcol(df, allcounts, ct = counttype)
  elif counttype = 'stop':
    # To get which columns correspond to stop mutations as labelled by snpEff.
    stopcol = [val for val in snpcol if 'fs' in val].extend([val for val in snpcol if '*' in val]) # makes a list subset from snpcol that disreguards mutations that do not cause early stop to protein.
    stopcoldf = df[stopcol] # USE this to get counts for mutations that disrupt the protein formation by gene in question.
    tdf = stopcoldf.transpose().replace(0,np.nan).groupby('GeneName').count()
    stopcounts = tdf.transpose()
    addgenecountcol(df, stopcounts, ct = counttype)
  else:
    print("Gene list type must be either 'all' or 'stop'."

def yearsubsetquery(df, colname = 'serotype', colvalue='IV', ncols=-8): # default search for type 'IV' 'serotype' and returns last 8 columns.
  querydf = df.query(colname+" == '"+colvalue+"'")
  yeardf = querydf.groupby('year').count
  colnameloc = list(df.columns).index(colname)
  finaldf = pd.concat([t4year.iloc[:,colnameloc],t4year.iloc[:,ncols:]], axis=1)
  return finaldf

