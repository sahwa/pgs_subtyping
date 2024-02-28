import gwaslab as gl
import pandas as pd
import sys
import gzip
import re

sumstats_path = sys.argv[1]
build = sys.argv[2]

def getrsIDCol(ss):
  for column in ss.data.columns:
    if ss.data.SNPID.astype(str).str.contains('rs\d+', regex=True, na=False).any():
      return column
  else:
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


def parsefiles(sumstats_path):
  output_path = sumstats_path.split(".")[0]
  stem = re.search(pattern, sumstats_path).group(0)
  output_path = output_path.replace("sumstats/", "sumstats/ldsc/")
  output_path = output_path + + "_ldsc.txt.gz"

  ref_rsid_vcf_19 = "/well/ckb/users/aey472/projects/pgs_subtype/data/GCF_000001405.25.gz"
  ref_rsid_vcf_38 = "/well/ckb/users/aey472/projects/pgs_subtype/data/GCF_000001405.40.gz"

  sep = guess_separator(sumstats_path)

  # Use the first command line argument as the input to gl.Sumstats
  ss = gl.Sumstats(sumstats_path, fmt="auto", sep=sep)
  #ss.basic_check()
  #ss.harmonize()

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
        ss.assign_rsid(ref_rsid_vcf = ref_rsid_vcf_19, chr_dict = gl.get_number_to_NC(build="19"), ncores=4)
        ss.data = ss.data.rename(columns={'rsID': 'SNP'})
    elif build == "38":
        ss.assign_rsid(ref_rsid_vcf = ref_rsid_vcf_38, chr_dict = gl.get_number_to_NC(build="38"), ncores=4)
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
    
parsefiles(sumstats_path)
