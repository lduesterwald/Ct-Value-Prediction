import sys
import os
import numpy as np
import pandas as pd
import pickle

# this function parses any parameters passed in through the command line or sets them to a default value
# parameters:
#    args: the list of arguments passed in through the command line
# returns: genomes_dir, genome_name, kmr_size, csv_path, model_name, dictionary_name
def parseParams(args, start_dir):
    # setting default values for parameters:
    genomes_dir = start_dir + "/" # (-g) the directory containing the genome to predict the Ct value for
    kmc_out_dir = genomes_dir + "kmc_output/" # (-k) the directory in which the output files of KMC are stored
    kmr_size = 10 # (-s) the size of the k-mer to run KMC with. Must be identical to the size used to construct the model
    dictionary_name = start_dir + "/kmr_dictionary.pkl" # (-i) the name of the stored dictionary (k-mer : column number)
    model_name = start_dir + "/ct_model.sav" # (-m) the name of the stored Ct value prediction model

    # required_parameter:
    genome_name = "" # (-n) required, the name of the file containing the genome (must be in .fasta format)
    csv_path = "" # (-c) required, the path to the .csv file containg the MCoV-id and instrument for all files in genomes_dir

    # parsing any parameters passed in through the command line
    for i in range(len(args)):
        if(args[i] == "-h" or args[i] == "--help"):
            print(helpOption())
            sys.exit()
        if (i == len(args) - 1):
            break
        elif (args[i] == "-g" or args[i] == "--genomes_dir"):
            genomes_dir = args[i + 1]
            genomes_dir = os.path.abspath(genomes_dir) + "/"
            kmc_out_dir = genomes_dir + "kmc_output/"
        elif (args[i] == "-k" or args[i] == "--kmc_out_dir"):
            kmc_out_dir = args[i + 1]
            kmc_out_dir = os.path.abspath(kmc_out_dir) + "/"
        elif (args[i] == "-s" or args[i] == "--kmr_size"):
            kmr_size = args[i + 1]
        elif(args[i] == "-c"or args[i] == "--csv_path"):
            csv_path = args[i + 1]
            csv_path = os.path.abspath(csv_path)
        elif (args[i] == "-i" or args[i] == "--dictionary_name"):
            dictionary_name = args[i + 1]
            dictionary_name = os.path.abspath(dictionary_name)
        elif (args[i] == "-m" or args[i] == "--model_name"):
            model_name = args[i + 1]
            model_name = os.path.abspath(model_name)
        elif (args[i] == "-n" or args[i] == "--genome_name"):
            genome_name = args[i + 1]



    # exitting the script if the required parameters were not passed in
    if (genome_name == "" or csv_path == ""):
        print("Error: required parameter not entered (csv_parh (-c) or genome_name (-n))")
        sys.exit()


    return genomes_dir, kmc_out_dir, genome_name, kmr_size, csv_path, model_name, dictionary_name



# returns a string of all the options for the script if the script was called with -h or --help
def helpOption():
    s = "-g --genomes_dir:\tthe directory containing the fasta files to train the model, the KMC executable package, and the kmc.sh script. The default is './'"
    s+= "\n-k --kmc_out_dir:\tthe directory in which the output files of KMC are stored"
    s+= "\n-s --kmr_size:\tthe size of the k-mer to run KMC with. The default is 10"
    s+= "\n-m --model_name:\tthe name of the stored Ct value prediction model. The default is 'ct_model.sav'"
    s+= "\n-i --dictionary_name:\tthe name of the dictionary { k-mer : column number } used to create the DataFrame that the model was trained on(must be a .pkl file). The default is 'kmr_dictionary.pkl'"
    s+= "\n-n --genome_name:\tthe genome fasta file to predict the Ct value of (must be named <genome_id>.fasta). There is no default for this option."
    s+= "\n-c --csv_path:\tthe path to the metadata .csv file with the genome_id and testing instrument of the genome There is no default for this option."
    return s


# runs KMC on the genome file to compute a list of the unique k-mers to be used by the model to predict the Ct value
# parameters:
#    genomes_dir: the directory containing genome_name (the name of the genome file in .fasta format)
#    kmr_size: the size of the k-mers to be used when running KMC (should be the same size as used to create the model)
def runKMC(genomes_dir, genome_name, kmr_size, kmc_out_dir):
    # run KMC on genome_name
    out_file = genome_name.replace(".fasta", "_kmc")

    # running KMC:
    os.chdir(genomes_dir)
    cmd = 'kmc.sh' + ' ' + str(kmr_size) + ' ' + genome_name + ' ' + out_file + ' ' + kmc_out_dir
    os.system(cmd)

    # deleting the temprorary files created by kmc:
    os.system('rm ' + out_file + '.kmc_pre')
    os.system('rm ' + out_file + '.kmc_suf')

    # returning the name of the output file of KMC containing the frequency of unique 10-mers
    return genome_name.strip(".fasta") + "_kmc." + str(kmr_size) + ".kmrs"


# creates a 1D numpy array from the k-mer counts of the input genome with the same features as the matrix the model was trained on
# parameters:
#    kmc_out_dir: the directory containing kmr_output_file (the output file of KMC with the list of unique k-mers)
#    kmr_dictionary: the dictionary of k-mer : column number used to construct the matrix the model was trained on
#    csv_path: the path to the file containing the instruments and MCoV-ids
def createRow(kmc_out_dir, kmr_output_file, kmr_dictionary, kmr_size, csv_path):
    # initalizing the row to have all 0s:
    row = []
    for i in range(len(kmr_dictionary)):
        row.append(0)

    # reads in the output file of kmc one k-mer at a time and updates the corresponding spot in row
    file = open((kmc_out_dir + kmr_output_file), "r")
    for aline in file:
        # getting the k-mer and k-mer frequency for every line:
        kmr = aline.split()[0]
        frequency = aline.split()[1]
        # getting the column number of the k-mer from the dictionary:
        col_num = kmr_dictionary[kmr]
        # updating the right index with the frequency of the k-mer:
        row[col_num] = frequency

    # adding the instrument:
    genome_id = kmr_output_file.strip("_kmc." + str(kmr_size) + ".kmrs")
    instrument = getInstrument(csv_path, genome_id) # getting the instrument corresponding to the MCoV-id

    for i in range(3):
        row.insert(0, 0) # adding 3 columns initialized to 0 for the instrument

    if (instrument == "ALINITY"):
        row[0] = 1
    elif (instrument == "PANTHER"):
        row[1] = 1
    elif (instrument == "CEPHEID"):
        row[2] = 1

    # returning the row as a 1D numpy array
    return np.array(row).reshape(1, -1)


# returns the instrument corresponding to the passed in MCoV-id
# parameters:
#    csv_path: the directory to the file containing information on the MCoV-ids and the instrument
#    mcovid: the MCoV-id for which to find the instruments
# returns: the instrument
def getInstrument(csv_path, genome_id):
    data_labels = pd.read_csv(csv_path)
    for index, row in data_labels.iterrows():
        if (row["Genome ID"] == genome_id):
            return row["Instrument"]



# main function
# predicts the Ct value of a genome
#  runs KMC on the genome to get a list of unique k-mers
#  creates a 1D numpy array with the same features as the matrix the model was trained on
#  predicts the Ct value of the array using the stored model
def main(argv):
    # current working directory:
    start_dir = os.getcwd()

    args = sys.argv
    # reads in parameters passed in by user through the command line or setting paramters to default values
    genomes_dir, kmc_out_dir, genome_name, kmr_size, csv_path, model_name, dictionary_name = parseParams(args, start_dir)

    # running KMC on the file containing the genome to predict
    kmc_output_file_name = runKMC(genomes_dir, genome_name, kmr_size, kmc_out_dir)
    print("--predictCt.py-- ran KMC on genome file")

    # opening the dictionary:
    kmr_dictionary = pickle.load(open(dictionary_name, "rb"))
    # creating a 1D numoy array with the same features as the matrix the model was trained on:
    row = createRow(kmc_out_dir, kmc_output_file_name, kmr_dictionary, kmr_size, csv_path)
    print("--predictCt.py-- created numpy array from k-mer counts")

    # opening the model:
    model = pickle.load(open(model_name, 'rb'))

    # predicting the ct value of the row:
    ct_prediction = model.predict(row)
    print("--predictCt.py-- predicted Ct value of genome")

    # printing the Ct value prediction
    print("\n\nGenome: ", genome_name.strip(".fasta"), "  Ct value prediction:  ",  str(ct_prediction[0]))


# if this is the script called by python, run main function
if __name__ == '__main__':
	main(sys.argv)

