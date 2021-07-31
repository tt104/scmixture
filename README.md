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
		<li><a href="#dependencies">Dependencies</a></li>
		<li><a href="#instructions">Instructions</a></li>
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
	
<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/GBM.umapClusters-1.jpg" width="450" align="right"></a>
	
All cells are projected onto a 2D plane through UMAP embeddings allowing the clusters assigned by the model to be visualised.
	
> This graph is created under the directory results/plots
	
<br clear="right"/>
	
</details>

<details open="open">
	<summary><b>Cluster Means</b></summary>
	
For each cluster, an interactive graph is produced plotting the mean gene expression within that cluster against the mean expression across all clusters. This allows the identification of genes with unusual expression within the discovered sub-populations of cells.
	
> These graphs are created under the directory results/{data}_exp
	
<p align="center"><sub><b>Average Gene Expression of All Clusters Against Cluster 1</b></sub></p>
<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmixture/blob/assets/Images/Cluster_Means.gif" width="100%"></a>
	
</details>

<details open="open">
	<summary><b>g:Profiler</b></summary>
	
For the up and down regulated genes within each cluster, the results of gene set enrichment analysis are presented through interactive Manhattan plots.
	
> This graph is created under the directory results/plots
	
<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Manhattan.gif" width="100%"></a>
	
</details>

<a id="dependencies"></a>
<h2>Dependencies</h2>

> For more information about a dependency, click the dropdown button to the left.
> 
<details>
	<summary><a href="https://sylabs.io/guides/3.0/user-guide/installation.html"><b>Singularity</b></a></summary>
	<p><b>What is Singularity?</b></p>
	<p>Singularity is a container platform. Within a container, required programs, libraries, data and scripts can be packaged together and distributed to create a reproducable runtime environment.</p>
	<p><b>Installation</b></p>
	<p>Singularity runs on Linux, though can also be run on Windows and Mac through a virtual machine. For installation instructions, see: https://sylabs.io/guides/3.0/user-guide/installation.html</p>
</details>
<details>
	<summary><a href="https://snakemake.readthedocs.io/en/stable/getting_started/installation.html"><b>Snakemake</b></a></summary>
	<p><b>What is Snakemake?</b></p>
	<p>Snakemake is a workflow management system for reproducible and scalable data analysis through the creation of pipelines. The workflow is defined in a ‘Snakefile’ consisting of rules that represent how to create output files from input files.
	</p>
	<p><b>Installation</b></p>
	<p>To install see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html</p>
</details>

<a id="instructions"></a>
<h2>Instructions</h2>

> Before following the instructions, make sure you have installed the <a href="#dependencies">dependencies.</a>

<details open="open">
	<summary><b>Set Up</b></summary>
	
**1.** Download the code from the main brach
	
- Either directly from GitHub by pressing Code ➜ Download ZIP</li>

<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Extract_ZIP.png" width="400"></a>
	
- Or through a Linux terminal using the command below
	
```console
(base) user@terminal:~$ wget https://github.com/tt104/scmix/archive/refs/heads/main.zip
```

**2.** Extract the contents of zip file

- Either locate the file and pressing "Extract"

<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Extract_ZIP.png" width="600"></a>

- Or if using the terminal, use the command below, where `scmix-main.zip` is the file path of the zip file and `~/new_file_path` is the location where you want to save the file

```console
(base) user@terminal:~$ unzip scmix-main.zip -d ~/new_file_path
```

**3.** Once the project has been unzipped, navigate to the project folder, where `scmix-main` is the file path

```console
(base) user@terminal:~$ cd scmix-main
```

**4.** Build the singularity container `scmix.sif` through command:

```console
(base) user@terminal:~/scmix-main$ sudo singularity build containers/scmix.sif containers/scmix.def
```

<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Container_Build.gif" width="700"></a>

</details>
	
<details open="open">
	<summary><b>Dataset Configuration</b></summary>

**5.** Place one or more scRNA-seq csv dataset(s) into `scmix-main/data`
	
<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Data_Folder.png" width="600"></a>
	
**6.** Optionally specify the dataset you wish to cluster when the code is run by adding its name into config.yaml
* Do not include the file extension `.csv` in the name
* This will not work unless the data is located in `scmix-main/data`
	
<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Config_Update.png" width="800"></a>
	
</details>

<details open="open">
	<summary><b>Running the Code</b></summary>

**7.** Activate the snakemake environment through conda.
```console
(base) user@terminal:~/scmix-main$ conda activate snakemake
```

**8.** Run the code through the command below:
* Replace `n` in `--cores n` with the number of cores you wish to use
* Optionally add the command option `--config data=DATASET`, where `DATASET` is the name of your dataset, if you did not specify a dataset in the config.yaml file
```console
(snakemake) user@terminal:~/scmix-main$ snakemake --use-singularity --cores n
```

<a src="https://github.com/tt104/scmix/archive/refs/heads/main.zip"><img src="https://github.com/tt104/scmix/blob/assets/Images/Code_Run.gif" width="800"></a>

</details>
