<div align="center">
  <a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/Logo.png" width="600"></a>
	<p><b>A Bayesian mixture model approach to zero-inflation of counts for scRNA-seq data</b></p>
</div>

---

<details open="open">
	<summary><b>Table of Contents</b></summary>
	<ul>
		<li><a href="#about">About</a></li>
		<li><a href="#features">Features</a></li>
		<li><a href="#run">Running the Code</a></li>
	</ul>
</details>

---

<a id="about"></a>
<h2>About</h2>

<p>Single cell RNA-seq data exhibit large numbers of zero count values, that as demonstrated in our <a href="https://www.biorxiv.org/content/10.1101/2021.05.19.444841v2">paper</a>, for a subset of transcripts, be better modelled by a zero inflated negative binomial distribution.
We develop a novel Dirichlet process mixture model which employs both a mixture at the cell level to model multiple cell types, and a mixture of single cell RNA-seq counts at the transcript level to model the transcript specific zero-inflation of counts.
It is shown that this approach outperforms previous approaches that applied multinomial distributions to model single cell RNA-seq counts, and also performs better or comparably to existing top performing methods.
By taking a Bayesian approach we are able to build interpretable models of expression within clusters, and to quantify uncertainty in cluster assignments.
Applied to a publicly available data set of single cell RNA-seq counts of multiple cell types from the mouse cortex and hippocampus, we demonstrate how our approach can be used to distinguish sub-populations of cells as clusters in the data, and to identify gene sets that are indicative of membership of a sub-population.
The methodology is implemented as an open source Snakemake pipeline hosted on this repository.</p>

<div>
  <a href="https://www.biorxiv.org/content/10.1101/2021.05.19.444841v2"><img align="left" width="100" src="https://github.com/tt104/scmixture/blob/assets/Images/bioRxiv_Logo.png"></a>
  <p><a href="https://www.biorxiv.org/content/10.1101/2021.05.19.444841v2"><b>Identifying sub-populations of cells in single cell transcriptomic data – a Bayesian mixture model approach to zero-inflation of counts</b></a></p>
</div>

---

<a id="features"></a>
<h2>Features</h2>

After running the pipeline, the following graphs are produced allowing the results of clustering by the Bayesian nonparametric zero-inflated negative binomial model to be interpreted:

<details open="open">
	<summary><b>UMAP</b></summary>
	
<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/UMAP.jpg" width="450" align="right"></a>
	
All cells are projected onto a 2D plane through UMAP embeddings of their gene expression allowing the clusters assigned by the model to be visualised.
	
> This graph is created under the directory results/plots
	
<br clear="right"/>
	
</details>

<details open="open">
	<summary><b>Cluster Means</b></summary>
	
For each cluster, an interactive graph is produced plotting the mean gene expression within that cluster against the mean expression across all clusters. This allows the identification of genes with unusual expression within the discovered sub-populations of cells.
	
> These graphs are created under the directory results/{data}_exp
	
<p align="center"><sub><b>Average Gene Expression of All Clusters Against Cluster 1</b></sub></p>
<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/Cluster_Means.gif" width="100%"></a>
	
</details>

<details open="open">
	<summary><b>g:Profiler</b></summary>
	
For the up and down regulated genes within each cluster, the results of gene set enrichment analysis with gprofiler2 are presented through interactive Manhattan plots.
This enables the identification of significantly enriched biological functions and pathways from well established sources such as Gene Ontology (GO) and KEGG.
	
> This graph is created under the directory results/plots
	
<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/gProfiler.gif" width="100%"></a>
	
</details>

<a id="run"></a>
<h2>Running the Code</h2>

<h3>Dependencies</h3>

<details>
	<summary><b>Singularity</b></summary>
	<p>Singularity is a container platform. Within a container, required programs, libraries, data and scripts can be packaged together and distributed to create a reproducable runtime environment.</p>
	<p><b>Installation</b></p>
	<p>Singularity runs on Linux, though can also be run on Windows and Mac through a virtual machine. For installation instructions, see: https://sylabs.io/guides/3.0/user-guide/installation.html</p>
</details>
<details>
	<summary><b>Snakemake</b></summary>
	<p>Snakemake is a workflow management system for reproducible and scalable data analysis through the creation of pipelines. The workflow is defined in a ‘Snakefile’ consisting of rules that represent how to create output files from input files.
	</p>
	<p><b>Installation</b></p>
	<p>To install see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html</p>
</details>

<h3>Instructions</h3>

<details open="open">
	<summary><b>Set Up</b></summary>
	
**1.** Download the code from the main brach
	
- Either directly from GitHub by pressing Code ➜ Download ZIP</li>

<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/Download_ZIP.png" width="350"></a>
	
- Or through a Linux terminal using the command below
	
```console
(base) user@terminal:~$ wget https://github.com/tt104/scmixture/archive/refs/heads/main.zip
```

**2.** Extract the contents of zip file

- Either locate the file and pressing "Extract"

<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Extract_ZIP.png" width="600"></a>

- Or if using the terminal, use the command below, where `scmixture.zip` is the file path of the zip file and `~/new_file_path` is the location where you want to save the file

```console
(base) user@terminal:~$ unzip scmixture.zip -d ~/new_file_path
```

**3.** Once the project has been unzipped, navigate to the project folder, where `scmixture` is the file path

```console
(base) user@terminal:~$ cd scmixture
```

</details>
	
<details open="open">
	<summary><b>Configuration</b></summary>

**5.** Place one or more scRNA-seq csv dataset(s) into folder `scmixture/data`

- If the folder does not exist, you can create it using the command below
	
```console
(base) user@terminal:~$ mkdir data
```
	
<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Data_Folder.png" width="600"></a>
	
**6.** Optionally edit config/config.yaml to change the algorithm's settings
	
```yaml
# K-means initialisation - no. of clusters
kmeans: 40
# Controls the number of transcripts selected to consider in clustering
filter: 1000
# Burn-in, total number of MCMC chain iterations, thinning out interval and total number of chains
burnin: 500
iter: 1000
thin: 10
chains: 4
# model - nb or mult
model: "nb"
```
	
</details>

<details open="open">
	<summary><b>Running the Code</b></summary>

**7.** Activate the snakemake environment through conda.
```console
(base) user@terminal:~/scmixture$ conda activate snakemake
```

**8.** Run the code through the command below:
* Replace `n` in `--cores n` with the number of cores you wish to use
* Replace `DATASET` with the name of your dataset, ignoring the .csv file extension, e.g. `--config data=GBM`
```console
(snakemake) user@terminal:~/scmixture$ snakemake --use-singularity --cores n --config data=DATASET
```

<a src="https://github.com/tt104/scmixture/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Code_Run.gif" width="800"></a>

</details>
