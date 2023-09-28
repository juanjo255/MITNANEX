from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
import pandas as pd
from pysradb.sraweb import SRAweb

# THIS CODE IS TO RETRIEVE METADATA FROM MITOCHONDRIAL GENOME BY BIOPROJECT ACCESSION
# THE GOAL IS TO CHECK IF THESE MITOCHONDRIAL GENOMES HAVE SRA

# GET METADATA BY BIOPROJECT ACCESSION USIN GENOMEAPI
def get_metadata_genome_api (bioprojects) -> list:
    result = list()
    #CONNECT WITH GENOME API
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        assemblies = genome_api.assembly_descriptors_by_bioproject(bioprojects)
        #print (assemblies)
        assemblies_dict = assemblies.to_dict()
        
        # ITER THROUGH ASSEMBLIES
        for assembly in assemblies_dict["assemblies"]:
            # CHECK KEYS WHICH STORE SRA EXIST
            if not ("biosample" in assembly["assembly"].keys()):
                continue
            elif not("sample_ids" in assembly["assembly"]["biosample"].keys()):
                continue
            for record in assembly["assembly"]["biosample"]["sample_ids"]:
                if "SRA" in record.values():
                    #print(assembly["assembly"]["assembly_accession"])
                    #print (assembly["assembly"]["bioproject_lineages"])
                    #print (assembly["assembly"]["biosample"]["sample_ids"])
                    #print (assembly["assembly"]["biosample"]["description"]["organism"]["organism_name"])
                    result.append((assembly["assembly"]["assembly_accession"],
                                assembly["assembly"]["bioproject_lineages"],
                                assembly["assembly"]["biosample"]["sample_ids"],
                                assembly["assembly"]["biosample"]["description"]["organism"]["organism_name"]
                                ))
        return result

# PARSE CSV FILE DOWNLOADED FROM NCBI
def parse_data (file) -> list:
    data = pd.read_csv(file, delimiter=",")
    bioprojects=data["BioProject"].to_list()
    return bioprojects

# ASSOCIATE GENOME ACCESSION TO A RECORD IN THE INITIAL DATASET
def get_more_info_filtered_mito (created_file, original_file):
    my_data= pd.read_excel(created_file, index_col=0)
    data_original= pd.read_csv(original_file, delimiter=",")
    data_original= data_original[data_original["Assembly"].notna()]
    #print (data_original.head())
    new_data= my_data.set_index(0).join(data_original.set_index('Assembly'))
    new_data.rename(columns = {1:'bioprojects', 2:'SRA_accession', 3:'organism_name_bioproject'}, inplace = True)
    print (new_data.head())
    new_data.to_excel("mito_ncbi_final.xlsx")

# GO THROUGH EVERY RECORD AND FILTER MITOCONDRIAL GENOMES BY SRA
# TRIGGER FUNCTION
def filt_mito():
    result = list()
    data = parse_data("organelles.csv")
    
    # ITERATE EVERY 500 FILES SINCE API DOES NOT SUPPORT MORE THAN THAT
    for span in range (0, len (data), 500):
        bioprojects=data[span:span+500]
    #bioprojects = ["PRJNA48091"]
        result = result + get_metadata_genome_api(bioprojects)
    
    dataframe= pd.DataFrame(result)
    dataframe.to_excel("mito_ncbi.xlsx")
    print ("MITOCHONDRIAS FILTERED BY SRA")
    
def get_sra_info_filtered_mito():
    db = SRAweb()
    sra= pd.read_excel("mito_ncbi_final.xlsx", index_col=0)
    #print (sra.sort_values(by=3).iloc[:,:3])
    sra_list = list()
    for sras in sra["SRA_accession"]:
        sras= sras.strip('][').split(', ')
        #print (sras)
        for index in range(0,len(sras)):
            if "SRA" in sras[index]:
                #print (sras[index+1].split("'")[3])
                sra_list.append(sras[index+1].split("'")[3])
                break
    #print(sra_list)
    df = db.sra_metadata(sra_list)
    df= df[['organism_name',"instrument",'instrument_model',
            'total_size','run_total_bases']]
    print(df.sort_values(by='organism_name').set_index('organism_name').head())
    #print(df.head())
    

if __name__ == "__main__":
    #filt_mito()
    #get_more_info_filtered_mito("mito_ncbi.xlsx", 'organelles.csv')
    get_sra_info_filtered_mito()