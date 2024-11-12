# Differential Genea Expression Analysis of Glioblastoma Multiforme (GBM) Tissue

### Rationale
In this analysis, I aim to identify genes that are differentially expressed in Glioblastoma Multiforme (GBM) tumors compared to normal tissue. Understanding differential gene expression in GBM may reveal insights into tumor biology, potential biomarkers, or therapeutic targets. I will perform the following series of steps:
 that include data normalization, application of statistical tests to identify significant genes, and visualization of results.


### Processing TCGA DATA for Bulk Rna-seq Analysis
* The Cancer Genome Atlas (TCGA) is a comprehensive and coordinated effort to accelerate our understanding of the molecular basis of cancer through the application of genome analysis technologies, including large-scale genome sequencing. 
* Data was collected from the GDC (Genomic Data Commons Data Portal) 
    * Glioblastoma Multiforme Project (Project ID: TCGA-GBM)
* Glioblastoma Multiforme dataset
    * Glioblastoma Multiforme (GBM) is a fast-growing type of malignant brain tumor that is the most common brain tumor in adults.
    * Patients with GBM have a poor prognosis and usually survive less than 15 months following diagnosis. Currently there are no effective long-term treatments for this disease.
### Objective 
* In this project I will be performing Bulk RNA-Seq Analysis of 157 brain tumor samples and 5 normal tissue samples. The goal of the analysis his ellucidate signaficantlly differential expressed genes. 


#### Load, Inspect, and Prepare Data
* Process data by assay type before running EDA and RNA-seq analaysis. Each tissue sample has a corresponding file_id, file_name, and project_id. However the files are nested and need to be placed into a single data frame for analysis.


```python
file_manifest = pd.read_csv(
    "data/gdc_sample_sheet.2024-08-24.tsv",
    sep="\t",
    header=0,
    names=[
        "file_id",
        "file_name",
        "data",
        "data_type",
        "project_id",
        "case_id",
        "sample_id",
        "sample_type",
    ],
    skiprows=[97],  # remove duplicate sample
)
allowed_sample_types = ["Primary Tumor", "Solid Tissue Normal"]
file_manifest = file_manifest[file_manifest["sample_type"].isin(allowed_sample_types)]
file_manifest.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>file_id</th>
      <th>file_name</th>
      <th>data</th>
      <th>data_type</th>
      <th>project_id</th>
      <th>case_id</th>
      <th>sample_id</th>
      <th>sample_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>69531d59-6508-4074-a63a-27b047121681</td>
      <td>03ddea6f-d377-4e2a-a59d-91d5c5b6e4ed.rna_seq.a...</td>
      <td>Transcriptome Profiling</td>
      <td>Gene Expression Quantification</td>
      <td>TCGA-GBM</td>
      <td>TCGA-27-2521</td>
      <td>TCGA-27-2521-01A</td>
      <td>Primary Tumor</td>
    </tr>
    <tr>
      <th>1</th>
      <td>e676c661-cae3-4755-a90f-e632a3d82850</td>
      <td>113df549-f143-4c50-8e18-7c2159bda4eb.rna_seq.a...</td>
      <td>Transcriptome Profiling</td>
      <td>Gene Expression Quantification</td>
      <td>TCGA-GBM</td>
      <td>TCGA-19-1390</td>
      <td>TCGA-19-1390-01A</td>
      <td>Primary Tumor</td>
    </tr>
    <tr>
      <th>2</th>
      <td>56c0c169-5480-4d14-b0d7-bd4214d0e7cf</td>
      <td>d088027b-e691-4962-bf1c-85abd24ac76c.rna_seq.a...</td>
      <td>Transcriptome Profiling</td>
      <td>Gene Expression Quantification</td>
      <td>TCGA-GBM</td>
      <td>TCGA-27-1830</td>
      <td>TCGA-27-1830-01A</td>
      <td>Primary Tumor</td>
    </tr>
    <tr>
      <th>3</th>
      <td>69f7dc91-32d1-4359-850f-860aab36e433</td>
      <td>0ce1eb02-a90c-484a-b861-c359562d3be1.rna_seq.a...</td>
      <td>Transcriptome Profiling</td>
      <td>Gene Expression Quantification</td>
      <td>TCGA-GBM</td>
      <td>TCGA-32-1970</td>
      <td>TCGA-32-1970-01A</td>
      <td>Primary Tumor</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ba9afa8e-405c-49b7-bc66-8ab1af1f942e</td>
      <td>dad61e18-e3f1-4beb-b3c3-ae434e35af2d.rna_seq.a...</td>
      <td>Transcriptome Profiling</td>
      <td>Gene Expression Quantification</td>
      <td>TCGA-GBM</td>
      <td>TCGA-06-0190</td>
      <td>TCGA-06-0190-01A</td>
      <td>Primary Tumor</td>
    </tr>
  </tbody>
</table>
</div>



#### Distrubution of Sample Types


```python
sample_type = file_manifest["sample_type"].value_counts().to_dict()
plt.bar(*zip(*sample_type.items()), align="center")
# Set Title
plt.title("Distrubution of Sample Type")
# Set x-axis label
plt.xlabel("Sample Type")
# Set y-axis label
plt.ylabel("Sample Type Count")

plt.show()
```


    
![png](/output_files/output_8_0.png)
    


Create framework to properly match file ids, file names, case ids and sample ids into the same data frame to get read counts for each sample.


```python
file_id = file_manifest["file_id"].values
file_names = file_manifest["file_name"].values
case_id = file_manifest["case_id"].values
sample_id = file_manifest["sample_id"].values
sample_type = file_manifest["sample_type"].values
condition = {i: k for i, k in zip(sample_id, sample_type)}

mainfest_elements = zip(file_id, file_names, case_id, sample_id)


manifest_map = {
    id: {
        "file_name": fn,
        "case_id": c_id,
        "sample_id": s_id,
    }
    for id, fn, c_id, s_id in mainfest_elements
}


assay_type = [
    "unstranded",
    "stranded_first",
    "stranded_second",
    "tpm_unstranded",
    "fpkm_unstranded",
    "fpkm_uq_unstranded",
]


assays = {a_type: {s_id: [] for s_id in sample_id} for a_type in assay_type}
```

The following creates a pickle file of processed RNA-seq data (by assay type), gene ids, and gene names.


```python
path = "data/processedRnaSeq.pickle"
check_file = os.path.isfile(path)

if check_file:
    print("File exsist already, loading now")
    with open("data/processedRnaSeq.pickle", "rb") as f:
        db = pkl.load(f)
else:
    for id in manifest_map:
        file_name = manifest_map[id]["file_name"]
        sample_id = manifest_map[id]["sample_id"]
        partial_data = pd.read_csv(
            f"data\gdc_final\{id}\{file_name}",
            sep="\t",
            header=1,
            skiprows=[2, 3, 4, 5],
        )
        for a_t in assay_type:
            assays[a_t][sample_id] = partial_data[a_t].values.tolist()

    gene_id = partial_data["gene_id"].values.tolist()
    gene_names = partial_data["gene_name"].values.tolist()

    for a_t in assay_type:
        assays[a_t]["gene_id"] = gene_id

    db = {}
    db["assays"] = assays
    db["gene_id"] = gene_id
    db["gene_names"] = gene_names

    with open("data/processedRnaSeq.pickle", "wb") as f:
        pkl.dump(db, f)

assay_data = db["assays"]  # Dictionary of RNA Seq Data
```

    File exsist already, loading now
    

Create a dataframe with gene ids, gene names, and gene lengths used in the Illumina platforms that completed the RNA sequcencing. 


```python
gene_info = pd.read_csv(
    "data/wgs.ASCAT.gene_level.copy_number_variation.tsv",
    sep="\t",
    header=0,
    usecols=["gene_id", "gene_name", "start", "end"],
)
gene_info["gene_length"] = gene_info.apply(lambda x: x["end"] - x["start"], axis=1)
gene_info.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_id</th>
      <th>gene_name</th>
      <th>start</th>
      <th>end</th>
      <th>gene_length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000223972.5</td>
      <td>DDX11L1</td>
      <td>11869</td>
      <td>14409</td>
      <td>2540</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000227232.5</td>
      <td>WASH7P</td>
      <td>14404</td>
      <td>29570</td>
      <td>15166</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ENSG00000278267.1</td>
      <td>MIR6859-1</td>
      <td>17369</td>
      <td>17436</td>
      <td>67</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ENSG00000243485.5</td>
      <td>MIR1302-2HG</td>
      <td>29554</td>
      <td>31109</td>
      <td>1555</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ENSG00000284332.1</td>
      <td>MIR1302-2</td>
      <td>30366</td>
      <td>30503</td>
      <td>137</td>
    </tr>
  </tbody>
</table>
</div>



#### Assay Choice Rationale
* Stranded_first, refers to data where the first strand of the RNA transcript is kept, preserving the directionality of transcription. This allows you to identify which DNA strand the RNA was transcribed from.
* This is particularly important if you're studying genes that overlap on opposite strands. This however is beyond the scope of this project. 



```python
def get_rna_data(data: dict, assay_type: str) -> pd.DataFrame:
    """
    This function returns desired rna seq count data in a dataframe
    """
    df = pd.DataFrame.from_dict(data[assay_type]).set_index("gene_id")
    return df


count_data = get_rna_data(assay_data, "stranded_first")
count_data.shape
count_data.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>TCGA-27-2521-01A</th>
      <th>TCGA-19-1390-01A</th>
      <th>TCGA-27-1830-01A</th>
      <th>TCGA-32-1970-01A</th>
      <th>TCGA-06-0190-01A</th>
      <th>TCGA-02-2486-01A</th>
      <th>TCGA-14-0790-01B</th>
      <th>TCGA-14-0789-01A</th>
      <th>TCGA-06-0645-01A</th>
      <th>TCGA-76-4926-01B</th>
      <th>...</th>
      <th>TCGA-27-2519-01A</th>
      <th>TCGA-26-5134-01A</th>
      <th>TCGA-06-0675-11A</th>
      <th>TCGA-06-1804-01A</th>
      <th>TCGA-14-1829-01A</th>
      <th>TCGA-06-0184-01A</th>
      <th>TCGA-32-2638-01A</th>
      <th>TCGA-12-3652-01A</th>
      <th>TCGA-28-5209-01A</th>
      <th>TCGA-76-4925-01A</th>
    </tr>
    <tr>
      <th>gene_id</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000000003.15</th>
      <td>174</td>
      <td>704</td>
      <td>2065</td>
      <td>1815</td>
      <td>1897</td>
      <td>3069</td>
      <td>1841</td>
      <td>1309</td>
      <td>1401</td>
      <td>3808</td>
      <td>...</td>
      <td>2326</td>
      <td>387</td>
      <td>306</td>
      <td>1692</td>
      <td>7138</td>
      <td>3416</td>
      <td>3838</td>
      <td>3262</td>
      <td>2334</td>
      <td>2684</td>
    </tr>
    <tr>
      <th>ENSG00000000005.6</th>
      <td>17</td>
      <td>1</td>
      <td>4</td>
      <td>13</td>
      <td>2</td>
      <td>7</td>
      <td>5</td>
      <td>1</td>
      <td>5</td>
      <td>33</td>
      <td>...</td>
      <td>13</td>
      <td>3</td>
      <td>3</td>
      <td>3</td>
      <td>66</td>
      <td>7</td>
      <td>16</td>
      <td>1</td>
      <td>5</td>
      <td>4</td>
    </tr>
    <tr>
      <th>ENSG00000000419.13</th>
      <td>921</td>
      <td>517</td>
      <td>544</td>
      <td>382</td>
      <td>580</td>
      <td>718</td>
      <td>402</td>
      <td>296</td>
      <td>652</td>
      <td>593</td>
      <td>...</td>
      <td>764</td>
      <td>98</td>
      <td>514</td>
      <td>604</td>
      <td>1109</td>
      <td>679</td>
      <td>851</td>
      <td>865</td>
      <td>671</td>
      <td>1555</td>
    </tr>
    <tr>
      <th>ENSG00000000457.14</th>
      <td>504</td>
      <td>544</td>
      <td>467</td>
      <td>383</td>
      <td>368</td>
      <td>436</td>
      <td>292</td>
      <td>262</td>
      <td>379</td>
      <td>548</td>
      <td>...</td>
      <td>382</td>
      <td>137</td>
      <td>386</td>
      <td>412</td>
      <td>492</td>
      <td>388</td>
      <td>574</td>
      <td>267</td>
      <td>500</td>
      <td>598</td>
    </tr>
    <tr>
      <th>ENSG00000000460.17</th>
      <td>428</td>
      <td>553</td>
      <td>377</td>
      <td>316</td>
      <td>432</td>
      <td>298</td>
      <td>232</td>
      <td>187</td>
      <td>322</td>
      <td>366</td>
      <td>...</td>
      <td>282</td>
      <td>94</td>
      <td>159</td>
      <td>310</td>
      <td>483</td>
      <td>310</td>
      <td>407</td>
      <td>356</td>
      <td>434</td>
      <td>549</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 162 columns</p>
</div>



#### Sample read count Distrubution


```python
# calcualte total counts per sample and log transform counts
total_counts = count_data.sum(axis=0).to_frame()
total_counts.columns = ["Read Count per Sample"]

log_counts = count_data.apply(lambda x: np.log2(x + 1))

sns.histplot(
    data=total_counts,
    x="Read Count per Sample",
    kde=True,
    # bins=50,
    # bw_adjust=5,
)
plt.show()
```


    
![png](/output_files/output_18_0.png)
    



```python
h_clustering = linkage(log_counts.T, "ward")
plt.figure(figsize=(20, 10))
dendrogram(h_clustering, labels=count_data.columns)
plt.xticks(rotation=90)
plt.ylabel("Distance")
plt.title("Hierarchical Clustering of Samples")
plt.tight_layout()
plt.show()
```


    
![png](/output_files/output_19_0.png)
    


The image above shows the results of hierarchical clustering, which can be visualized via a dendrogram. Samples clustered together are more similar to each other, and the length of the branches (vertical lines) connecting clusters represents the distance or dissimilarity between clusters. The above image shows that the normal tissue samples are more closely related , than the brain tumor samples. 

#### TPM Normalization Method
* TPM normalizes for gene length, but it also ensures that the sum of the TPM values in each sample is the same, making it more comparable across samples. TPM is calculated by normalizing read counts by gene length first, and then scaling by sequencing depth.
* Before applying tpm normalization some genes are missing length information so need to populate using Esmbl API.


```python
gene_map = {gi: gn for gn, gi in zip(db["gene_names"], db["gene_id"])}

gene_length_map = {
    gi: len
    for gi, len in zip(gene_info["gene_id"].values, gene_info["gene_length"].values)
}

gene_to_fix = []
for gene in count_data.index:
    if gene not in gene_length_map:
        gene_to_fix.append(gene)

gene_name = []
for gene in gene_to_fix:
    gene_name.append(gene_map[gene])


def get_gene_name(name):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{name}?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    stable_gene = r.json()[0]["id"]
    return stable_gene


def get_gene_length(ensmbl):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensmbl}?expand=1"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    length = r.json()["Transcript"][0]["length"]
    # print(repr(decoded))
    return length


fix_genes = {}
for idx, name in enumerate(gene_name):
    stable = get_gene_name(name)
    length = get_gene_length(stable)
    fix_genes[gene_to_fix[idx]] = length
gene_length_map.update(fix_genes)
```

TPM Normalized read count data.


```python
count_data = count_data[count_data.sum(axis=1) > 0]  # filter out 0 reads
count_data_tpm = count_data.apply(lambda x: x / gene_length_map[x.name], axis=1)
count_data_tpm = count_data_tpm.apply(lambda x: x / (x.sum() / 1_000_000), axis=0)
print(count_data_tpm.shape)
count_data_tpm.head()
```

    (56543, 162)
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>TCGA-27-2521-01A</th>
      <th>TCGA-19-1390-01A</th>
      <th>TCGA-27-1830-01A</th>
      <th>TCGA-32-1970-01A</th>
      <th>TCGA-06-0190-01A</th>
      <th>TCGA-02-2486-01A</th>
      <th>TCGA-14-0790-01B</th>
      <th>TCGA-14-0789-01A</th>
      <th>TCGA-06-0645-01A</th>
      <th>TCGA-76-4926-01B</th>
      <th>...</th>
      <th>TCGA-27-2519-01A</th>
      <th>TCGA-26-5134-01A</th>
      <th>TCGA-06-0675-11A</th>
      <th>TCGA-06-1804-01A</th>
      <th>TCGA-14-1829-01A</th>
      <th>TCGA-06-0184-01A</th>
      <th>TCGA-32-2638-01A</th>
      <th>TCGA-12-3652-01A</th>
      <th>TCGA-28-5209-01A</th>
      <th>TCGA-76-4925-01A</th>
    </tr>
    <tr>
      <th>gene_id</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ENSG00000000003.15</th>
      <td>3.362298</td>
      <td>11.785398</td>
      <td>40.523517</td>
      <td>37.115810</td>
      <td>43.434346</td>
      <td>77.321314</td>
      <td>56.744073</td>
      <td>34.730010</td>
      <td>26.870381</td>
      <td>71.147335</td>
      <td>...</td>
      <td>38.945130</td>
      <td>19.850860</td>
      <td>4.443098</td>
      <td>27.895187</td>
      <td>143.514382</td>
      <td>51.281836</td>
      <td>66.472749</td>
      <td>68.244192</td>
      <td>54.680790</td>
      <td>46.785993</td>
    </tr>
    <tr>
      <th>ENSG00000000005.6</th>
      <td>0.283101</td>
      <td>0.014427</td>
      <td>0.067648</td>
      <td>0.229103</td>
      <td>0.039464</td>
      <td>0.151987</td>
      <td>0.132813</td>
      <td>0.022865</td>
      <td>0.082644</td>
      <td>0.531350</td>
      <td>...</td>
      <td>0.187582</td>
      <td>0.132616</td>
      <td>0.037540</td>
      <td>0.042624</td>
      <td>1.143583</td>
      <td>0.090563</td>
      <td>0.238816</td>
      <td>0.018030</td>
      <td>0.100951</td>
      <td>0.060089</td>
    </tr>
    <tr>
      <th>ENSG00000000419.13</th>
      <td>9.679106</td>
      <td>4.707071</td>
      <td>5.805967</td>
      <td>4.248487</td>
      <td>7.222417</td>
      <td>9.838194</td>
      <td>6.738782</td>
      <td>4.271157</td>
      <td>6.800986</td>
      <td>6.025665</td>
      <td>...</td>
      <td>6.957055</td>
      <td>2.733903</td>
      <td>4.058974</td>
      <td>5.415698</td>
      <td>12.126600</td>
      <td>5.543756</td>
      <td>8.015984</td>
      <td>9.842070</td>
      <td>8.549586</td>
      <td>14.741862</td>
    </tr>
    <tr>
      <th>ENSG00000000457.14</th>
      <td>2.810925</td>
      <td>2.628465</td>
      <td>2.645061</td>
      <td>2.260543</td>
      <td>2.431899</td>
      <td>3.170447</td>
      <td>2.597654</td>
      <td>2.006311</td>
      <td>2.098006</td>
      <td>2.955112</td>
      <td>...</td>
      <td>1.846029</td>
      <td>2.028246</td>
      <td>1.617646</td>
      <td>1.960459</td>
      <td>2.855063</td>
      <td>1.681160</td>
      <td>2.869342</td>
      <td>1.612222</td>
      <td>3.380926</td>
      <td>3.008613</td>
    </tr>
    <tr>
      <th>ENSG00000000460.17</th>
      <td>0.554730</td>
      <td>0.620937</td>
      <td>0.496226</td>
      <td>0.433431</td>
      <td>0.663438</td>
      <td>0.503581</td>
      <td>0.479629</td>
      <td>0.332780</td>
      <td>0.414231</td>
      <td>0.458663</td>
      <td>...</td>
      <td>0.316696</td>
      <td>0.323405</td>
      <td>0.154850</td>
      <td>0.342800</td>
      <td>0.651353</td>
      <td>0.312146</td>
      <td>0.472807</td>
      <td>0.499554</td>
      <td>0.681984</td>
      <td>0.641884</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 162 columns</p>
</div>



#### Quality Control and Filtering
* Filtering out lowly expressed genes across all samples. For more significant differential gene expression analysis. This will help to reduce noise in the data and hopefully focus on relevant genes. Will also improve computing efficiency.


```python
thresholds = np.arange(0, 5, 0.1)
genes_retained = []
for min_tpm in thresholds:
    mask = (count_data_tpm > min_tpm).sum(axis=1) >= 5
    genes_retained.append(mask.sum())  # Count number of genes retained

plt.figure(figsize=(20, 10))
plt.plot(thresholds, genes_retained, marker="o", color="green")
plt.axvline(
    x=1.0, color="red", linestyle="--", label="TPM = 1"
)  # put line at 1 as rough heuristic
plt.xlabel("Threshold (TPM)")
plt.ylabel("Num Genes Retained")
plt.title("Filtering for Genes to Retain")
plt.legend()
plt.show()
```


    
![png](/output_files/output_26_0.png)
    


Based on chart above i'll filter genes w/ expression threshold of < 1 TPM, 1 TPM threshold is a common heuristic, again the purpose of filtering genes is to increase the statiscal power of the DE analysis. By adding a TPM  filter I'll only consider 30,463 genes out of the original 56,543.


#### Principal Component Analysis of Sample Types


```python
# Separating out the features
x = count_data.loc[:,].T.values

# Separating out the target
metadata = pd.DataFrame(
    zip(count_data.T.index, [condition[x] for x in count_data.T.index]),
    columns=["Sample", "Condition"],
).set_index("Sample")
y = metadata["Condition"].values


x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)

principalComponents = pca.fit_transform(x)

finalDf = pd.DataFrame(
    data=principalComponents, columns=["principal component 1", "principal component 2"]
)
finalDf["Target"] = y
```


```python
fig = plt.figure(figsize=(8, 5))
pc1 = pca.explained_variance_ratio_[0]
pc2 = pca.explained_variance_ratio_[1]
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(f"Principal Component 1 {pc1*100:.1f}%", fontsize=15)
ax.set_ylabel(f"Principal Component 2 {pc2*100:.1f}%", fontsize=15)
ax.set_title("Two Component PCA of Samples", fontsize=20)

targets = ["Primary Tumor", "Solid Tissue Normal"]
colors = ["r", "g"]
for target, color in zip(targets, colors):
    indicesToKeep = finalDf["Target"] == target
    ax.scatter(
        finalDf.loc[indicesToKeep, "principal component 1"],
        finalDf.loc[indicesToKeep, "principal component 2"],
        c=color,
        s=25,
    )
ax.grid()
ax.legend(targets)
```




    <matplotlib.legend.Legend at 0x1bb72d67790>




    
![png](/output_files/output_30_1.png)
    


* The Principal Component Graph (PCA) graph shown above is the result of two-component analysis showing the two sample groups being study. The first principal component captures the most variablity 11.6% and the second component captures 5.9%. Take note I used the raw count data to create pca plot. Using TPM normalized counts would change the amount of variance explain by each component. Filtering the genes would also alter the variance as well. 
* The plot visualizes the relationship between the samples based on their pc scores. You can the seperation between the two sample groups and potential clustering within samples that are more similar to one another. There also looks like their may be some outliers in the two sample groups as there are samples that seem to be away from clusters.

#### Differential Expression Analysis


```python
mask = (count_data_tpm > 1).sum(axis=1) >= 5  # filtering huerstic
filtered_data = count_data_tpm[mask]
print(filtered_data.shape)

normal_sample = metadata.index[metadata["Condition"] == "Solid Tissue Normal"]
tumor_sample = metadata.index[metadata["Condition"] == "Primary Tumor"]
results = []

for gene in filtered_data.index:
    base_mean = filtered_data.loc[gene, :].mean()
    normal = filtered_data.loc[gene, normal_sample]
    tumor = filtered_data.loc[gene, tumor_sample]
    normal_mean = normal.mean()
    tumor_mean = tumor.mean()
    fc = (tumor_mean - normal_mean) / normal_mean  # fold change
    log2fc = np.log2((tumor_mean + 1) / (normal_mean + 1))  # Adding 1 to avoid log of 0
    t_stat, p_val = ttest_ind(normal, tumor)
    results.append(
        {
            "geneid": gene,
            "fc": fc,
            "base_mean": base_mean,
            "log2fc": log2fc,
            "t_stat": t_stat,
            "p_val": p_val,
        }
    )
# filtered_data.head()
```

    (30463, 162)
    

    C:\Users\osuch\AppData\Local\Temp\ipykernel_12912\244121798.py:15: RuntimeWarning: divide by zero encountered in scalar divide
      fc = (tumor_mean - normal_mean) / normal_mean  # fold change
    

Final filtering for signficant genes. Another heurstic I'll be using is eleminating genes with a base mean less than 10 also I will filter for adusted pvalues < 0.05 and absolute log2 fold changes > 1.


```python
results_df = pd.DataFrame(results)
results_df["padj"] = multipletests(results_df["p_val"], method="fdr_bh")[1]
results_df["abs_log2fc"] = results_df["log2fc"].abs()
results_df["symbol"] = results_df.apply(lambda row: gene_map[row["geneid"]], axis=1)

results = results_df[results_df.base_mean >= 10]
sigs = results[(results.padj < 0.05) & (abs(results.log2fc) > 1)]
print(
    f"After running differential gene expression analysis and filtering for significant genes based on the above criteria we are left with {len(sigs)} genes"
)
sigs.head()
```

    After running differential gene expression analysis and filtering for significant genes based on the above criteria we are left with 2092 genes
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>geneid</th>
      <th>fc</th>
      <th>base_mean</th>
      <th>log2fc</th>
      <th>t_stat</th>
      <th>p_val</th>
      <th>padj</th>
      <th>abs_log2fc</th>
      <th>symbol</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000000003.15</td>
      <td>9.739217</td>
      <td>45.924452</td>
      <td>3.159546</td>
      <td>-3.759455</td>
      <td>0.000238</td>
      <td>0.002085</td>
      <td>3.159546</td>
      <td>TSPAN6</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ENSG00000001036.14</td>
      <td>1.784385</td>
      <td>21.818367</td>
      <td>1.370717</td>
      <td>-3.319123</td>
      <td>0.001118</td>
      <td>0.007373</td>
      <td>1.370717</td>
      <td>FUCA2</td>
    </tr>
    <tr>
      <th>19</th>
      <td>ENSG00000002549.12</td>
      <td>2.311172</td>
      <td>34.013409</td>
      <td>1.636995</td>
      <td>-3.286768</td>
      <td>0.001246</td>
      <td>0.008065</td>
      <td>1.636995</td>
      <td>LAP3</td>
    </tr>
    <tr>
      <th>20</th>
      <td>ENSG00000002586.20</td>
      <td>11.054143</td>
      <td>75.951787</td>
      <td>3.402884</td>
      <td>-4.047245</td>
      <td>0.000080</td>
      <td>0.000842</td>
      <td>3.402884</td>
      <td>CD99</td>
    </tr>
    <tr>
      <th>43</th>
      <td>ENSG00000004455.17</td>
      <td>3.132240</td>
      <td>11.604690</td>
      <td>1.732940</td>
      <td>-5.006171</td>
      <td>0.000001</td>
      <td>0.000022</td>
      <td>1.732940</td>
      <td>AK2</td>
    </tr>
  </tbody>
</table>
</div>




```python
# view distribution of scores
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.histplot(sigs["abs_log2fc"], bins=50, kde=True)
plt.title("Distribution of Absolute log2FC from Differential Expression Analysis")
plt.xlabel("Absolute log2Fc")
plt.ylabel("Frequency (# of genes)")
plt.show()
```


    
![png](/output_files/output_36_0.png)
    


The distribution graph displays the absolute log2 fold change (log2FC) values, which measure the magnitude of gene expression changes between tumor and normal samples. Most genes fall within the lower range of absolute log2FC (0.5 to 1.5), indicating relatively small expression changes. The highest frequency of genes has a log2FC value around 1, reflecting a roughly 2-fold expression change. As log2FC increases, the frequency of genes decreases, showing that large changes in expression (log2FC > 2) are less common. The graph follows a right-skewed distribution, with a long tail suggesting that while small fold changes are typical, some genes still exhibit substantial expression changes (log2FC > 4), albeit less frequently.

### Volcano Plots of DE genes


```python
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=results_df,
    x="log2fc",
    y="padj",
    hue="log2fc",
    palette="viridis",
    alpha=0.9,
    edgecolor=None,
)
plt.title("Volcano Plot of All Genes")
plt.axhline(y=0.05, color="red", linestyle="-", linewidth=1)
plt.axvline(x=1, color="blue", linestyle="-", linewidth=1)
plt.axvline(x=-1, color="blue", linestyle="-", linewidth=1)
plt.xlabel("log2 Fold Change")
plt.ylabel("Adjusted P-value")
plt.legend(title="log2 Fold Change", loc="lower left")
plt.show()
```


    
![png](/output_files/output_39_0.png)
    


The volcano plot displays the log2 fold change (x-axis) to show gene expression differences between tumor and normal samples, with positive values indicating upregulation in tumors and negative values indicating downregulation. The y-axis represents the adjusted p-value, where lower values indicate higher statistical significance, and a red line marks the significance cutoff. Genes with larger fold changes are highlighted in brighter colors (yellow/green for upregulated, blue/purple for downregulated). Vertical blue lines denote fold change thresholds, indicating biologically significant changes. Genes on the far right upregulated and genes to the far left are downregulated, while those near the center exhibit little to no differential expression.


```python
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=sigs,
    x="log2fc",
    y="padj",
    hue="log2fc",
    palette="viridis",
    alpha=0.9,
    edgecolor=None,
)
plt.title("Volcano Plot of Significant Genes")
plt.axhline(y=0.05, color="red", linestyle="-", linewidth=1)
plt.axvline(x=1, color="blue", linestyle="-", linewidth=1)
plt.axvline(x=-1, color="blue", linestyle="-", linewidth=1)
plt.xlabel("log2 Fold Change")
plt.ylabel("Adjusted P-value")
plt.legend(title="log2 Fold Change", loc="lower left")
plt.show()
```


    
![png](/output_files/output_41_0.png)
    


This volcano plot focuses on significant genes and provides a clearer view of the genes with both substantial fold changes and strong statistical significance. There appears to be a stronger pattern of upregulation in tumor samples compared to downregulation, with fewer downregulated genes meeting the significance threshold. 


```python
# Extract the list of gene names for DEGs

gene_list = sigs[sigs["abs_log2fc"] >= 2].symbol.tolist()

# Initialize g:Profiler
gpro = GProfiler(return_dataframe=True)

# Perform GO analysis using the significant gene list
go_results = gpro.profile(organism="hsapiens", query=gene_list)

# Display the first few results
go_results.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source</th>
      <th>native</th>
      <th>name</th>
      <th>p_value</th>
      <th>significant</th>
      <th>description</th>
      <th>term_size</th>
      <th>query_size</th>
      <th>intersection_size</th>
      <th>effective_domain_size</th>
      <th>precision</th>
      <th>recall</th>
      <th>query</th>
      <th>parents</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GO:CC</td>
      <td>GO:0030054</td>
      <td>cell junction</td>
      <td>1.248963e-27</td>
      <td>True</td>
      <td>"A cellular component that forms a specialized...</td>
      <td>2230</td>
      <td>385</td>
      <td>119</td>
      <td>22097</td>
      <td>0.309091</td>
      <td>0.053363</td>
      <td>query_1</td>
      <td>[GO:0110165]</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GO:CC</td>
      <td>GO:0031982</td>
      <td>vesicle</td>
      <td>3.023965e-24</td>
      <td>True</td>
      <td>"Any small, fluid-filled, spherical organelle ...</td>
      <td>4004</td>
      <td>385</td>
      <td>159</td>
      <td>22097</td>
      <td>0.412987</td>
      <td>0.039710</td>
      <td>query_1</td>
      <td>[GO:0043227]</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GO:CC</td>
      <td>GO:0045202</td>
      <td>synapse</td>
      <td>2.348335e-23</td>
      <td>True</td>
      <td>"The junction between an axon of one neuron an...</td>
      <td>1468</td>
      <td>385</td>
      <td>89</td>
      <td>22097</td>
      <td>0.231169</td>
      <td>0.060627</td>
      <td>query_1</td>
      <td>[GO:0030054]</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GO:CC</td>
      <td>GO:0005737</td>
      <td>cytoplasm</td>
      <td>6.472974e-18</td>
      <td>True</td>
      <td>"The contents of a cell excluding the plasma m...</td>
      <td>12345</td>
      <td>385</td>
      <td>301</td>
      <td>22097</td>
      <td>0.781818</td>
      <td>0.024382</td>
      <td>query_1</td>
      <td>[GO:0005622, GO:0110165]</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GO:CC</td>
      <td>GO:0022626</td>
      <td>cytosolic ribosome</td>
      <td>1.099537e-16</td>
      <td>True</td>
      <td>"A ribosome located in the cytosol." [GOC:mtg_...</td>
      <td>117</td>
      <td>385</td>
      <td>24</td>
      <td>22097</td>
      <td>0.062338</td>
      <td>0.205128</td>
      <td>query_1</td>
      <td>[GO:0005829, GO:0005840]</td>
    </tr>
  </tbody>
</table>
</div>




```python
# perform pathway analysis
pathway_analysis_results = gpro.profile(
    organism="hsapiens", query=gene_list, sources=["KEGG", "REAC"]
)
# display the results
pathway_analysis_results.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source</th>
      <th>native</th>
      <th>name</th>
      <th>p_value</th>
      <th>significant</th>
      <th>description</th>
      <th>term_size</th>
      <th>query_size</th>
      <th>intersection_size</th>
      <th>effective_domain_size</th>
      <th>precision</th>
      <th>recall</th>
      <th>query</th>
      <th>parents</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>REAC</td>
      <td>REAC:R-HSA-156842</td>
      <td>Eukaryotic Translation Elongation</td>
      <td>1.668325e-16</td>
      <td>True</td>
      <td>Eukaryotic Translation Elongation</td>
      <td>94</td>
      <td>245</td>
      <td>24</td>
      <td>10842</td>
      <td>0.097959</td>
      <td>0.255319</td>
      <td>query_1</td>
      <td>[REAC:R-HSA-72766]</td>
    </tr>
    <tr>
      <th>1</th>
      <td>REAC</td>
      <td>REAC:R-HSA-156902</td>
      <td>Peptide chain elongation</td>
      <td>9.262719e-16</td>
      <td>True</td>
      <td>Peptide chain elongation</td>
      <td>90</td>
      <td>245</td>
      <td>23</td>
      <td>10842</td>
      <td>0.093878</td>
      <td>0.255556</td>
      <td>query_1</td>
      <td>[REAC:R-HSA-156842]</td>
    </tr>
    <tr>
      <th>2</th>
      <td>REAC</td>
      <td>REAC:R-HSA-192823</td>
      <td>Viral mRNA Translation</td>
      <td>1.485265e-14</td>
      <td>True</td>
      <td>Viral mRNA Translation</td>
      <td>90</td>
      <td>245</td>
      <td>22</td>
      <td>10842</td>
      <td>0.089796</td>
      <td>0.244444</td>
      <td>query_1</td>
      <td>[REAC:R-HSA-168273]</td>
    </tr>
    <tr>
      <th>3</th>
      <td>REAC</td>
      <td>REAC:R-HSA-1799339</td>
      <td>SRP-dependent cotranslational protein targetin...</td>
      <td>1.702170e-14</td>
      <td>True</td>
      <td>SRP-dependent cotranslational protein targetin...</td>
      <td>113</td>
      <td>245</td>
      <td>24</td>
      <td>10842</td>
      <td>0.097959</td>
      <td>0.212389</td>
      <td>query_1</td>
      <td>[REAC:R-HSA-72766]</td>
    </tr>
    <tr>
      <th>4</th>
      <td>REAC</td>
      <td>REAC:R-HSA-2408557</td>
      <td>Selenocysteine synthesis</td>
      <td>4.064438e-14</td>
      <td>True</td>
      <td>Selenocysteine synthesis</td>
      <td>94</td>
      <td>245</td>
      <td>22</td>
      <td>10842</td>
      <td>0.089796</td>
      <td>0.234043</td>
      <td>query_1</td>
      <td>[REAC:R-HSA-2408522]</td>
    </tr>
  </tbody>
</table>
</div>




```python
sns.barplot(x="intersection_size", y="name", data=pathway_analysis_results.head(10))
plt.title("Top 10 Pathways")
plt.xlabel("Term size")
plt.ylabel("Pathway")
plt.show()
```


    
![png](/output_files/output_45_0.png)
    


The Gene Ontology (GO) enrichment analysis chart you provided highlights pathways enriched in your dataset, likely from the genes identified as differentially expressed in GBM compared to normal tissue. Here’s an interpretation of the key pathways in relation to Glioblastoma Multiforme (GBM):

1. Eukaryotic Translation Elongation and Peptide Chain Elongation
These pathways are involved in protein synthesis, specifically in elongating the amino acid chain during translation. Overactivation of translation processes is common in cancers, including GBM, where rapid cell growth and division demand increased protein synthesis.

2. Viral mRNA Translation and Host Translation Machinery Modulation
Cancer cells, similar to viral infections, often hijack the host's translational machinery for their own growth. GBM cells may exploit these pathways to support uncontrolled proliferation. Interestingly, the presence of pathways associated with viral translation highlights how tumor cells might mimic viral strategies for survival and growth.

3. SRP-Dependent Cotranslational Protein Targeting to Membrane
This pathway is involved in targeting proteins to the endoplasmic reticulum membrane, a critical process for producing membrane-bound and secretory proteins. In GBM, these processes could facilitate cell signaling changes and tumor growth, as many oncogenic proteins are membrane-associated.

4. Selenocysteine Synthesis
Selenoproteins play a role in redox regulation and oxidative stress. Cancer cells, including GBM cells, often manipulate oxidative stress pathways to promote survival and resistance to treatments.

5. Nonsense-Mediated Decay (NMD)
NMD is a quality control pathway that degrades mRNA with premature stop codons, preventing translation of potentially harmful truncated proteins. In cancer, alterations in NMD can affect gene expression stability, leading to the dysregulation of genes involved in cell growth and survival.

6. Response to Amino Acid Deficiency (EIF2AK4/GCN2)
Cancer cells frequently experience metabolic stress, including nutrient deprivation. GBM cells may upregulate pathways that allow them to survive under amino acid deprivation, supporting their adaptation to the tumor microenvironment.

These enriched pathways suggest that GBM tumors may depend heavily on translational and metabolic regulation to support tumor growth and survival. The findings indicate potential targets for therapy, such as inhibitors of protein synthesis pathways or regulators of oxidative stress. These pathways could also inform biomarkers for GBM, especially in distinguishing GBM cells' adaptations compared to normal brain tissue.


```python
ranking = sigs[["symbol", "t_stat"]].dropna().sort_values("t_stat", ascending=False)
manual_set = {"things": sigs.symbol}


pre_res = gp.prerank(
    rnk=ranking,
    gene_sets=["GO_Biological_Process_2021", manual_set],
    seed=6,
    permutation_num=100,
)
```


```python
out = []
for term in list(pre_res.results):
    out.append(
        [
            term,
            pre_res.results[term]["fdr"],
            pre_res.results[term]["es"],
            pre_res.results[term]["nes"],
        ]
    )
out_df = (
    pd.DataFrame(out, columns=["Term", "fdr", "es", "nes"])
    .sort_values("fdr")
    .reset_index(drop=True)
)
out_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Term</th>
      <th>fdr</th>
      <th>es</th>
      <th>nes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GO_Biological_Process_2021__protein targeting ...</td>
      <td>0.0</td>
      <td>-0.273625</td>
      <td>-2.138840</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GO_Biological_Process_2021__mRNA processing (G...</td>
      <td>0.0</td>
      <td>-0.356528</td>
      <td>-2.219854</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GO_Biological_Process_2021__RNA splicing, via ...</td>
      <td>0.0</td>
      <td>-0.383762</td>
      <td>-2.336474</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GO_Biological_Process_2021__mRNA splicing, via...</td>
      <td>0.0</td>
      <td>-0.361982</td>
      <td>-2.202237</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GO_Biological_Process_2021__protein N-linked g...</td>
      <td>0.0</td>
      <td>-0.493344</td>
      <td>-2.341351</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>268</th>
      <td>GO_Biological_Process_2021__mitochondrion orga...</td>
      <td>1.0</td>
      <td>0.217855</td>
      <td>0.699544</td>
    </tr>
    <tr>
      <th>269</th>
      <td>GO_Biological_Process_2021__regulation of phos...</td>
      <td>1.0</td>
      <td>0.172608</td>
      <td>0.576438</td>
    </tr>
    <tr>
      <th>270</th>
      <td>GO_Biological_Process_2021__receptor-mediated ...</td>
      <td>1.0</td>
      <td>0.254994</td>
      <td>0.764255</td>
    </tr>
    <tr>
      <th>271</th>
      <td>GO_Biological_Process_2021__DNA repair (GO:000...</td>
      <td>1.0</td>
      <td>0.147265</td>
      <td>0.532827</td>
    </tr>
    <tr>
      <th>272</th>
      <td>GO_Biological_Process_2021__negative regulatio...</td>
      <td>1.0</td>
      <td>0.156316</td>
      <td>0.602732</td>
    </tr>
  </tbody>
</table>
<p>273 rows × 4 columns</p>
</div>



### Gene Set Enrichment Analysis Example


```python
gseaplot(
    rank_metric=pre_res.ranking,
    term=out_df.sort_values("nes", ascending=False).iloc[:5].Term,
    **pre_res.results[out_df.sort_values("nes", ascending=False).iloc[5].Term]
)
```




    [<Axes: xlabel='Gene Rank', ylabel='Ranked metric'>,
     <Axes: >,
     <Axes: >,
     <Axes: ylabel='Enrichment Score'>]




    
![png](/output_files/output_50_1.png)
    


### Conclusion 
The analysis identified several genes with significant differential expression in GBM versus normal tissue. The upregulated genes may play roles in tumor growth, invasion, and resistance to therapy, while downregulated genes could be associated with disrupted normal cell function in the tumor environment. These findings provide a foundation for further investigation into the molecular mechanisms of GBM and may guide the development of targeted therapies or biomarkers for early detection and treatment monitoring.

Future Directions This analysis has identified several promising avenues for further research. Investigating specific genes within these enriched pathways could provide insights into their precise roles in GBM pathology. Additionally, experimental studies validating the effects of inhibiting key pathways—such as translation, cellular communication, and metabolic stress response—would help assess their therapeutic potential. Further analysis on patient outcomes associated with these gene signatures could also refine biomarker discovery for prognosis and treatment response in GBM.


```python

```
