import gwaslab as gl
import pandas as pd
import numpy as np
import os
import sys
import gzip
import re

def check_empty_columns(df, threshold=0.75, sample_size=1000000):
  total_rows = len(df)
  sample_size = min(sample_size, total_rows)  # Ensure sample size is not larger than the total number of rows
    
  random_indices = np.random.choice(total_rows, size=sample_size, replace=False)
  sampled_df = df.iloc[random_indices]
    
  for column in sampled_df.columns:
    missing_percentage = sampled_df[column].isnull().mean()
    if missing_percentage >= threshold:
      return False  
    
  return True

def checkrsid(df):
  required_match_percentage = 0.8
  subset = df.SNP.head(1000).astype(str)
  matches = subset.str.contains('rs\d+', regex=True, na=False).sum()

  ## check rsid columns
  if (matches / len(subset)) >= required_match_percentage:
    print("Found rsID SNPs in SNP column")
    return True
  else:
    print("SNP column doesn't seem to contain rsID. Exiting")
    return False

def getrsIDCol(ss):
  required_match_percentage = 0.75  # 75% of the rows must match the pattern
  for column in ss.data.columns:
    subset = ss.data[column].head(1000).astype(str)
        
    if subset.isnull().all() or len(subset) == 0:
      continue  # Skip this column if it has no valid data
        
    matches = subset.str.contains('rs\d+', regex=True, na=False).sum()
    if matches > 0 and (matches / len(subset)) >= required_match_percentage:
      return column
    return "not found"

def guess_separator(filepath):
  with gzip.open(filepath, 'rt') as f: 
    first_line = f.readline()
  separators = [',', '\t', ';', ' ']

  separator_counts = {sep: first_line.count(sep) for sep in separators}
  guessed_separator = max(separator_counts, key=separator_counts.get)

  return guessed_separator if separator_counts[guessed_separator] > 0 else None

def getbuild(ss, given_build):
  if given_build != "NA":
    if given_build == "37":
      return "19"
    if given_build == "38":
      return "38"
  else:
    print("No build found. Automatically detecting using HapMap SNPS")
    ss.infer_build()
    return ss.meta['gwaslab']['genome_build']


def check_and_reorder_dataframe(df):
  # Predefined sets of columns
  cols_sets = [
    ['SNP', 'A1', 'A2', 'BETA', 'SE'],  
    ['SNP', 'A1', 'A2', 'OR', 'SE'],    
    ['SNP', 'A1', 'A2', 'BETA', 'P'],   
    ['SNP', 'A1', 'A2', 'OR', 'P']      
    ]

  colnames_set = set(df.columns)

    # Find the first matching set and reorder DataFrame
  for cols_set in cols_sets:
    if set(cols_set).issubset(colnames_set):
      return df.loc[:, cols_set]

    # If no matching set is found, raise an error
  raise ValueError("None of the predefined column sets are fully present in the DataFrame columns.")


def writefile(ss, outpath):
  ss.to_csv(outpath, sep="\t", index=False, compression="gzip")


def main(args):

  sumstats_path = args[1]
  build = args[2]

  output_path = sumstats_path.split(".")[0]
  output_path = output_path.replace("sumstats/", "sumstats/ldsc/")
  output_path = output_path + "_ldsc.txt.gz"

  #### run some checks first we don't want to repeat if the file already exists

  if os.path.exists(output_path):
    print("Output file exists, checking contents")
    ldscfile = pd.read_csv(output_path, sep="\t", compression='gzip')
    
    ecols = check_empty_columns(ldscfile)
    checkrsidcols = checkrsid(ldscfile)

    if not ecols or not checkrsidcols:
      print("Something wrong with the output file.. Trying again")
    else:
      print("Input files seem OK. Exiting program")
      exit()
  
  ref_rsid_vcf_19 = "/well/ckb/users/aey472/projects/pgs_subtype/data/GCF_000001405.25.gz"
  ref_rsid_vcf_38 = "/well/ckb/users/aey472/projects/pgs_subtype/data/GCF_000001405.40.gz"

  sep = guess_separator(sumstats_path)

  # Use the first command line argument as the input to gl.Sumstats
  ss = gl.Sumstats(sumstats_path, fmt="auto", sep=sep)
  ss.basic_check()
  ss.harmonize()

  ## something here we need to fix the NEA and EA alleles so that we can use them for the adding rsID

  ## check whether there is a column which seems to contain the rsid
  
  print("Checking whether RSID col is present")
  rsIDcol = getrsIDCol(ss)
  print("Finished checking RSID cols");


  ## if not, then try and add one
  if rsIDcol == "not found":
    ## first we need to figure out the build of the dataset
    print("Didn't find any RSID col, first checking build\n")
    build = getbuild(ss, build)
    print(f"Finished checking build, identified as {build}\n")

    ## doing so depends on what the build of the data is
    if build == '19':
        ss.assign_rsid(ref_rsid_vcf = ref_rsid_vcf_19, chr_dict = gl.get_number_to_NC(build="19"), n_cores=4)
        ss.data = ss.data.rename(columns={'rsID': 'SNP'})
    elif build == "38":
        ss.assign_rsid(ref_rsid_vcf = ref_rsid_vcf_38, chr_dict = gl.get_number_to_NC(build="38"), n_cores=4)
        ss.data = ss.data.rename(columns={'rsID': 'SNP'})
    else:
        sys.exit("Build detected as neither 19 nor 38")
  ## if we did find it, then make sure it's named 'SNP'
  else:
    print(f"RSID cols identified as {rsIDcol}, renaming to SNP")
    ss.data = ss.data.rename(columns={rsIDcol: 'SNP'})

  ## rename the EA and NEA alleles to A1 and A2 
  print("Renaming allele columns")
  ss.data = ss.data.rename(columns={'EA': 'A1'})
  ss.data = ss.data.rename(columns={'NEA': 'A2'})

  ## finally drop any columns which aren't needed ##
  
  print("Finally checking the right columns are present")
  ss = pd.DataFrame(ss.data)
  ss = check_and_reorder_dataframe(ss)
  
  print("Writing output files")
  writefile(ss, output_path)
  print("Completed!")

if __name__ == "__main__":
  main(sys.argv)
