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

def search_column(data,searchcol, val, new_col = True, newname=0):
  if new_col == True:
    if newname ==0:
      newname = input('New Column Name: ')
    list = data[searchcol]
    data[newname] = [1 if val in row else np.nan for row in list] # iterates over rows adding a 1 if val found, otherwise adding NaN
    return data
  else:
    df = [1 if val in row else np.nan for row in list]
    return df

def total(data, add_percent = True):
  temp = data.transpose()
  data['Total'] = temp.count()
  if add_percent == True:
    n = len(data)-23 ### There are 23 descriptive rows before individual isolates are listed, other vcf files may vary. ###
    data['Percent'] = data['Total']/n*100
  return data.transpose()

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

def hotencoder(col,data):
  from sklearn.preprocessing import OneHotEncoder
  encoder = OneHotEncoder()
  enc = encoder.fit_transform(data[[col]])
  catlist = encoder.get_feature_names_out()
  print("Categories encoded:", catlist)
  df= pd.DataFrame(enc.toarray(), columns=catlist)
  return df
