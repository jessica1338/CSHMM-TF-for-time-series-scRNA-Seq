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


