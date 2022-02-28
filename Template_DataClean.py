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
  
