# CSHMM-TF-for-time-series-scRNA-Seq

This github repository provides the source code for the paper CSHMM-TF for time-series single-cell RNA-Seq data

## Before you begin
* Please make sure that you have **A Docker installation**. Follow the documentation here: https://docs.docker.com/engine/installation/

* Download the Docker file for building docker image: https://raw.githubusercontent.com/jessica1338/CSHMM-TF-for-time-series-scRNA-Seq/master/Dockerfile

## Install the Docker container

Build the docker image, make sure that the Dockerfile is under the current directory
```
docker build -t cshmm_tf_release .
```

Run and connect to the container:
```
 docker run -it cshmm_tf_release /bin/bash
```

## Run the initialization part within the container

Now within the container, run the example script:
```
python scdiff_init.py -d treutlein2014 -tf tfDNA_predicted_100.txt.update 
```

Sometimes the scdiff initialization result will be different than the result in our paper when you run on different environment
In this case you can skip this step use the provided initialization file "init_cluster_treutlein2014_lung.txt"

## Run the CSHMM-TF training

after the model initialization, you need to run the training script "run_TF.sh" using the following commend

```
sh run_TF.sh
```

Here we only run 2 iterations to save time, we provide the .pickle file for 10th iteration so that you can run the analysis step without training the full 10 epochs
you can use the same setting to run on your dataset

## Run the CSHMM-TF analysis 

You can use the following command to generate some of the result/visualization in paper
```
python example_train_and_analysis.py
```

You can also look at the ipython notebook: example_train_and_analysis.ipynb


## Run your data on docker container

To mount your data on a local disk to a location within the docker filesystem, use the ```-v``` option:

```
docker run -it -v [data folder on your comupter]:[data path you want on docker container] cshmm_tf_release /bin/bash
```
For example:
```
docker run -it -v ~/my_data:/my_data_dc cshmm_tf_release /bin/bash
```
Then you can access the files in ```~/my_data``` from ```/my_data_dc``` on docker container

# INPUTS AND PRE-PROCESSING (We use the same input format as SCDIFF, so most of the following description are from their github page: https://github.com/phoenixding/scdiff/blob/master/README.md)

CSHMM-TF takes two required input files (-d for data file and -tf for tf-target information), two optional files (-k/--cluster, -e/--etfListFile) and a few other optional parameters. 

* __-d__  
This specifies the single cell RNA-Seq expression data.  
If the RNA-Seq data is not processed, the instruction about how to calculate expression based on RNA-Seq raw reads can be found in many other studies, e.g (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/).
For example, users can use Tophat + Cufflink to calculate the gene expression in terms of FPKM.  Please refer to corresponding tools for instructions. 
Once we get the RNA-Seq gene expression, the expression data should be transformed to log space for example by log2(x+1) where x could represent the gene expression in terms of RPKM, FPKM or TPM depending
on what tools are used to precoess the RNA-Seq expression data.  
The input file has the following formatting requirements:
	* __Header Row__  
	First 3 columns are "Cells","Time","Label" and the remaining columns are gene names.   
	* __Data Rows__  
		* __1st column__: Cell ID, represents the ID for the cell.
		* __2nd column__: Cell time, Integer, represents the measurement time of the cell. 
		* __3rd column__: Cell label, represents the label of the cell (e.g cell type if known). In most cases, we don't have any prior knowledge of the cell type. In this case, use "NA" instead.
		Or, you can use any name you want to label each cell. We don't use this information in our model and it's only used to mark the cells with 
		the same attributes (any known attributes users are interested in, for example, cell type, time point, WT/Treatment, etc.) in the visualization. 
		Please avoid too many different labels, which will make the visualization very crowded. It's recommended to use <20 different labels. 
		If, however, users want to use more labels, please use the 'NA' as the label for all cells and use the cell IDs to infer the label composition of each node. 
		* __4th- columns__: Gene expression values.  
	
	Example input:     
	[example data file (lung dataset)](/treutlein2014)

* __-tf__  
This specifies the TF-gene interaction data.  In other words, it specifies the TF targets. 
Under the tf_dna directory, we provided a [human TF-gene interaction file](tf_dna/Human_TF_targets.txt) and a [mouse TF-gene interaction file](tf_dna/Mouse_TF_targets.txt) inferred using the strategy in our previous study (https://www.ncbi.nlm.nih.gov/pubmed/20219943). 
Although this TF-gene interactions are collected in human and mouse, they should be also able to apply to other close species.
Besides, in our previous work DREM (http://sb.cs.cmu.edu/drem/), we did collected the TF-gene interactions for common species including human, mouse, fry, E.coli, yeast, Arabidopsis. 
Please refer to  http://sb.cs.cmu.edu/drem/DREMmanual.pdf appendix B for complete details. 
Those TF-gene interaction files can be downloaded from our DREM software (https://github.com/phoenixding/idrem/tree/master/TFInput).
You might need to unzip and re-format the file to satisfy the requirements. The TF-gene interaction file has the following formatting requirements:  
 
	* __Header Row__  
	```
	TF	Gene	Input
	```
	* __Data Rows__  
		* __1st column__: TF ID (gene symbol)
		* __2rd column__: gene ID (gene symbol)
		* __3rd column__: Input, optional, the interaction strength between TF and target gene. If missing, by default it's 1.  
		This column is not used in scdiff. 
		 	
	Example file:   
	[example TF gene interaction file (we also use this file in our work)](/tfDNA_predicted_100.txt.update)

