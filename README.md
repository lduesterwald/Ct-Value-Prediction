# Ct-Value-Prediction
This repository contains scripts to train and score a Random Forest regression model to predict the Ct values of input genome data. It contains the following Python scripts:
* *runKMC.py* - running KMC (k-mer counter) on all genome data files in a directory
* *createDataFrame.py* - creating a feature matrix of k-mer frequencies
* *trainModel.py* - training and evaluating a model to predict the Ct values from the k-mer matrix
* *predictCt.py* - using the model to predict the Ct value of an individual genome.

This repo also includes the *sample* directory containing the data and model files for testing and running the scripts including instructions on how to download relevant data.


## Requirements
In order to run the scripts:
1. KMC must be installed and the kmc.sh script must be downloaded
2. Python and the required libraries must be installed
3. A directory containing genome data as *.fasta* files for the model training must be available

### 1. Installing KMC
This repo requires the KMC (k-mer counter) executable package (kmc, kmc_dump, and kmc_tools) to be installed. KMC can be downloaded [here](https://refresh-bio.github.io/). Furthermore, the kmc.sh script to run kmc and kmc_dump in tandem (found in this [GitHub repo](https://github.com/Tinyman392/GenomicModelCreator/blob/edf2c58a616a2bdb7bbaa19aee87dd7b795afcba/KMC/kmc.sh)) must be downloaded into the same directory as the KMC executive package. This directory path should be added to your PATH environment variable.

### 2. Python Packages
This repo is designed to be run with python3 (version 3.9.12).
This repo requires the following packages to be installed:
* numpy (version 1.21.5)
* pandas (version 1.4.2)
* sklearn (version 1.0.2)

### 3. Genome Data Directory
This repo requires a directory containing genome data files. The name of the directory can be passed into the scripts (option -g). The genome data files should have *.fasta* extension and be named <genome_id>.fasta. 

Two scripts in this repo (*createDataFrame.py* and *predictCt.py*) also require a comma-separated (*.csv*) metadata file with information about each genome. This file should contain at least the genome_id (corresponding to the file names), testing instrument, and Ct value for each genome file. If the file includes additional columns, they will be ignored by the scripts. The path to this file must be passed into the scripts (option -c). 



## Running the Scripts
The scripts in this repo implement a pipeline and should be run in the following order: 
1. *runKMC.py*
2. *createDataFrame.py*
3. *trainModel.py*
4. *predictCt.py*

This repo also contains the *ct_value_prediction.sh* bash script to run the entire pipeline.

### *runKMC.py*
The *runKMC.py* script is used to run KMC on every *.fasta* file in the genome directory and to move the output files produced by KMC containing the frequencies of every unique k-mer in the genome to a specified output directory. 

An example run would be:
~~~
python3 runKMC.py -g <genome directory>
~~~

This script takes in the following options:
* -g  --genomes_dir: Specify the directory containing the *.fasta* files to train the model. The default is the current directory “./”.
* -k --kmc_out_dir: Specify the directory to be created by the script for storing the output files produced by KMC containing k-mer frequencies. The default is “~/<genomes_dir>/kmc_output/”.
* -s --kmr_size: Specify the size of the k-mer with which to run KMC. The default is 10.


### *createDataFrame.py*
The *createDataFrame.py* script is used to create and store the k-mer matrix used to train the model as a Pandas DataFrame. The DataFrame represents the features consisting of the frequencies of every unique k-mer across all genomes and the testing instrument of the genome one-hot encoded. The DataFrame also contains the Ct values as labels. 

An example run would be:
~~~
python3 createDataFrame.py -g <genome directory>  -c ~/<metadata file>
~~~

The script takes in the following options:
* -g --genomes_dir: Specify the directory containing the *.fasta* files to train the model. This should be the same directory used for *runKMC.py*. The default is the current directory “./”.
* -k --kmc_out_dir:  Specify the directory to be created by the script for storing the output files produced by KMC containing k-mer frequencies. This should be the same directory used for *runKMC.py*. The default is “~/<genomes_dir>/kmc_output/”.
* -s --kmr_size: Specify the size of the k-mer with which to run KMC. This should be the same as used in *runKMC.py*. The default is 10.
* -c --csv_path: Specify the path to the comma-separated (*.csv*) metadata file containing the genome_id, testing instrument, and Ct value of each genome in the genome directory. There is no default for this option.
* -o --output_dir: Specify the directory for storing the k-mer DataFrame created by the script. If the output directory does not exist, it will be created by the script. The default is “~/<genomes_dir>/output/”. 
* -d --df_name: Specify the name that the k-mer DataFrame created by the script will be stored as. This must be a *.csv* file. The default is "kmr_df.csv".
* -i --dictionary_name: Specify the name that the dictionary of { k-mer : column number } used to create the DataFrame will be stored as. This must be a *.pkl* file. The default is "kmr_dictionary.pkl". 
An example dictionary used in the createDataFrame.py script is included in the sample folder.


### *trainModel.py*
The *trainModel.py* script trains and stores a Random Forest regression model using the DataFrame created by *createDataFrame.py* to predict the Ct value. The script then evaluates the model’s accuracy and calculates the R2 score, RMSE (root mean squared error), and the model’s accuracy within certain intervals and writes the results to an output file. 
An example model trained by this script is included in the sample folder. 

An example run would be:
~~~
python3 trainModel.py -o ~/<genome_directory>/output -m ct_prediction_model.sav
~~~

The script takes in the following options:
* -o --output_dir: Specify the directory containing the k-mer DataFrame that will be used for storing the model trained in this script. This must be the same as in *createDataFrame.py*. The default is “./output”.
* -d --df_name: Specify the name that the DataFrame created by *createDataFrame.py* was stored as. Must be a *.csv* file. The default is "kmr_df.csv".
* -m --model_name: Specify the name that the Ct value prediction model will be stored as. This must be a .sav file. The default is "ct_model.sav".
* -f --output_file_name: Specify the name of the output file for this script. This file will be created in the <output_dir> directory. The default is "output_file_trainModel".
* -ts --test_size: Specify the size of the test set to be used in the train_test_split during model training and evaluation. The default is 0.2.
* -nt --num_trees: Specify the ‘n_estimators’ (number of trees) parameter in the Random Forest regression model. The default was established through hyperparameter tuning and is 400.
* -td --tree_depth: Specify the ‘max_depth’ (tree depth) parameter in the Random Forest regression model. The default was established through hyperparameter tuning and is None.
* -rs --row_subsampling: Specify the ‘max_samples’ (row subsampling) parameter in the Random Forest regression model. The default was established through hyperparameter tuning and is 0.25.

An example output file would be as follows:
~~~
Model Accuracy:
R2: 0.64141243234
RMSE: 5.5342345432740915

Accuracy Within Intervals:
within 6: 0.7904761904761904
within 5: 0.6761904761904762
within 4: 0.580952380952381
within 3: 0.5047619047619047
within 2: 0.38095238095238093
within 1: 0.23853211009174313
~~~


### *predictCt.py*
The *predictCt.py* script takes in one genome *.fasta* file and predicts its Ct value using the model created by *trainModel.py*. The script runs KMC and uses the frequencies of the genome’s k-mers as features in the same way as the *createDataFrame.py* script.
This script also requires a *.csv* file containing the genome_id of the genome and the testing instrument.

An example run would be:
~~~
python3 predictCt.py -g <genome directory> -c ~/<metadata file> -m ct_prediction_model.sav -n <genome id>.fasta
~~~

The script takes in the following options:
* -g --genomes_dir: Specify the directory containing the genome data as a *.fasta* file to predict the Ct value. The default is the current directory “./”.
* -s --kmr_size: Specify the size of the k-mer with which to run KMC. This must be the same as was used to create the model. The default is 10.
* -c --csv_path: Specify the path to the *.csv* file containing information about the genome whose Ct value to predict. This file must contain the <genome_id> of the genome (matching the name of the file) and the testing instrument of the genome. There is no default for this option.
* -o --output_dir: Specify the directory containing the Ct value prediction model and the dictionary {k-mer : column number} created by *trainModel.py* and *createDataFrame.py*. The default is "~/<genomes_dir>/output”.
* -i --dictionary_name: Specify the name of the dictionary {k-mer : column number} used to create the DataFrame in *createDataFrame.py*. Must be a .pkl file. The default is “kmr_dictionary.pkl".
* -m --model_name: Specify the name of the Ct value prediction model created by the *trainModel.py* script. This should be a .sav file. The default is "ct_model.sav".
* -n --genome_name: Specify the genome *.fasta* file name for which to predict the Ct value. Must be in <genomes_dir>. This file should be named <genome_id>.fasta. There is no default for this option.



### *ct_value_prediction.sh*
The *ct_value_prediction.sh* script runs all 4 scripts in a sequence. This script takes in the union of the arguments of the individual component scripts. Running the script with the -h option will list all optional and required arguments. 

An example run would be:
~~~
bash ct_value_prediction.sh -g <genome directory> -c ~/<metadata file> -m ct_prediction_model.sav -n <genome id>.fasta
~~~
