import pandas as pd
from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from pysradb.sraweb import SRAweb
from tqdm import trange
import matplotlib.pyplot as plt

# PARSE CSV FILE DOWNLOADED FROM NCBI
def parse_data(file) -> list:
    data = pd.read_csv(file, delimiter=",")
    # bioprojects = list()
    # for i in range(len(data.index)):
    #     for replicon in data.loc[i,'Replicons'].split('; '):
    #         if 'mitochondrion' in replicon:
    #             bioprojects.append(data.loc[i,'BioProject'])
    bioprojects = data["bioproject_s"].to_list()
    return bioprojects

# GO THROUGH EVERY RECORD AND FILTER MITOCONDRIAL GENOMES BY SRA
def filt_mito(csv: str, output_file: str):
    print("FILTERING MITOCHONDRIAS BY SRA")
    result = list()
    data = parse_data(csv)

    # ITERATE EVERY 500 FILES SINCE API DOES NOT SUPPORT MORE THAN THAT
    for span in trange(0, len(data), 500):
        bioprojects = [i for i in data[span : span + 500] if isinstance(i, str)]
        # bioprojects = ["PRJNA48091"]
        result = result + [
            (*v, bioprojects[k])
            for k, v in enumerate(get_metadata_genome_api(bioprojects))
        ]
        # print(result)
    dataframe = pd.DataFrame(result)
    dataframe.rename(
        columns={
            0: "assembly_accesion",
            1: "bioprojects",
            2: "SRA_accession",
            3: "organism_name",
            4: "bioproject",
        },
        inplace=True,
    )
    dataframe.to_excel(output_file)
    return dataframe.head()

# GET DATA FROM EVERY BIOPROJECT
def get_metadata_genome_api(bioprojects) -> list:
    result = list()
    # CONNECT WITH GENOME API
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        assemblies = genome_api.assembly_descriptors_by_bioproject(bioprojects)
        assemblies_dict = assemblies.to_dict()
        if "assemblies" in assemblies_dict.keys():
            # ITER THROUGH ASSEMBLIES
            for assembly in assemblies_dict["assemblies"]:
                # CHECK THAT KEYS STORING SRA EXIST
                if not ("biosample" in assembly["assembly"].keys()):
                    continue
                elif not ("sample_ids" in assembly["assembly"]["biosample"].keys()):
                    continue
                for record in assembly["assembly"]["biosample"]["sample_ids"]:
                    if "SRA" in record.values():
                        # print (assembly["assembly"]["bioproject_lineages"])
                        # print (assembly["assembly"]["biosample"]["sample_ids"])
                        # print (assembly["assembly"]["biosample"]["description"]["organism"]["organism_name"])
                        result.append(
                            (
                                assembly["assembly"]["assembly_accession"],
                                assembly["assembly"]["bioproject_lineages"],
                                assembly["assembly"]["biosample"]["sample_ids"],
                                assembly["assembly"]["biosample"]["description"][
                                    "organism"
                                ]["organism_name"],
                            )
                        )
    return result

# ASSOCIATE GENOME ACCESSION TO A RECORD IN THE INITIAL DATASET
def extend_info_filt_data(created_file, original_file, output):
    my_data = pd.read_excel(created_file, index_col=0)
    data_original = pd.read_csv(original_file, delimiter=",")
    new_data = my_data.merge(
        data_original[["bioproject_s", "create_date_dt"]],
        right_on="bioproject_s",
        left_on="bioproject",
    )
    new_data.dropna(subset="SRA_accession", inplace=True)
    new_data.drop_duplicates(subset="bioproject", inplace=True)
    new_data.reset_index(drop=True, inplace=True)
    new_data.to_excel(output)
    return new_data.head()

# GET SRA METADATA OF FILTERED MITOCHONDRIAL GENOMES
def get_sra_info_filtered_mito(input_excel_file: str, output_file: str, original_file:str):
    db = SRAweb()
    sra = pd.read_excel(input_excel_file, index_col=0)
    data_original = pd.read_csv(original_file, delimiter=",")
    sra_list = list()
    for i in range(0, len(sra)):
        sras = (sra.loc[i, "SRA_accession"]).strip("][").split(", ")
        ## get the sra value
        for index in range(0, len(sras)):
            if "SRA" in sras[index]:
                sra_list.append(sras[index + 1].split("'")[3])
                break
    df = db.sra_metadata(sra_list)
    df = df[
        [
            "organism_name",
            "instrument",
            "instrument_model",
            "total_size",
            "run_accession",
            "bioproject",
        ]
    ]
    df = df.merge(
        data_original[["bioproject_s", "create_date_dt"]],
        right_on="bioproject_s",
        left_on="bioproject"
    )
    # df = df[
    #     (df["instrument"].isin(["GridION", "MinION", "PromethION"]))
    #     & (df["total_size"].astype(int) > 0)
    # ]
    # df['date'] = date_sra_list
    df = df[df["total_size"].astype(int) > 0]
    #df.drop_duplicates(subset="bioproject", inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.to_excel(output_file)
    return df


# DESCRIPTIVE STATISTICS
def analyse_dataset(file):
    df = pd.read_excel(file, index_col=0)
    df = (
        df.groupby(["instrument"])
        .agg({"instrument": "count"})
        .assign(percentage=lambda x: x.instrument / x.instrument.sum() * 100)
    )
    # df.to_excel("sra_metadata_analysis_results.xlsx")
    print(df)
    return df.head()


if __name__ == "__main__":
    filt_mito("datasets_metadata/wgs_selector_animal.csv", "datasets_metadata/sra_per_bioproject.xlsx")
    # extend_info_filt_data("datasets_metadata/sra_per_bioproject.xlsx", "datasets_metadata/wgs_selector.csv", "datasets_metadata/sra_per_bioproject.xlsx")
    get_sra_info_filtered_mito(
        "datasets_metadata/sra_per_bioproject.xlsx",
        "datasets_metadata/sra_metadata.xlsx",
        "datasets_metadata/wgs_selector_animal.csv"
    )
    analyse_dataset('datasets_metadata/sra_metadata.xlsx')
