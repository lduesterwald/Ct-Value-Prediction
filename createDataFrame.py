import sys
import os
import numpy as np
import pandas as pd
import pickle

# this function parses any parameters passed in through the command line or sets them to a default value
# parameters:
#    args: the list of arguments passed in through the command line
# returns: genomes_dir, kmc_out_dir, kmr_size, csv_path, output_dir, df_name, dictionary_name
def parseParams(args):
    # setting default values for parameters:
    genomes_dir = "./" # (-g) the directory with genomes as .fasta files and KMC executibles
    kmc_out_dir = genomes_dir + "kmc_output" # (-k) the directory to which to store the output files of KMC with k-mer counts
    kmr_size = 10 # (-s) the size k-mer to run KMC with
    output_dir = genomes_dir + "output" # (-o) the directory in which to store the final DataFrame and the dictionary used to create the DataFrame
    df_name = "kmr_df.csv" # (-d) the name to store the k-mer DataFrame created by this script as
    dictionary_name = "kmr_dictionary.pkl" # (-i) the name to store the dictionary (k-mer : column number) created by the script as

    # required parameter:
    csv_path = "" # (-c) required, the path to the .csv file containg the genome_id and instrument for all files in genomes_dir

    # parsing any parameters passed in through the command line
    for i in range(len(args)):
        if(args[i] == "-h" or args[i] == "--help"):
            print(helpOption())
            sys.exit()
        if (i == len(args) - 1):
            break
        elif (args[i] == "-g" or args[i] == "--genomes_dir"):
            genomes_dir = args[i + 1]
            if (genomes_dir.endswith("/") == False):
                genomes_dir = genomes_dir + "/"
            kmc_out_dir = genomes_dir + "kmc_output"
        elif (args[i] == "-k" or args[i] == "--kmc_out_dir"):
            kmc_out_dir = args[i + 1]
        elif (args[i] == "-s" or args[i] == "--kmr_size"):
            kmr_size = int(args[i + 1])
        elif (args[i] == "-c" or args[i] == "--csv_path"):
            csv_path = args[i + 1]
        elif (args[i] == "-o" or args[i] == "--output_dir"):
            output_dir = args[i + 1]
        elif (args[i] == "-d" or args[i] == "--df_name"):
            df_name = args[i + 1]
        elif (args[i] == "-i" or args[i] == "--dictionary_name"):
            dictionary_name = args[i + 1]


    # exitting the script if the required parameter (csv_path) was not passed in
    if (csv_path == ""):
        print("Error: csv_path (-c) (required parameter) not entered")
        sys.exit()

    return genomes_dir, kmc_out_dir, kmr_size, csv_path, output_dir, df_name, dictionary_name



# returns a string of all the options for the script if the script was called with -h or --help
def helpOption():
    s = "-g --genomes_dir:\tthe directory containing the fasta files to train the model, the KMC executable package, and the kmc.sh script. The default is './'"
    s+= "\n-k --kmc_out_dir:\tthe directory for storing the output files created by KMC. The default is ~/genomes_dir/kmc_out_dir"
    s+= "\n-s --kmr_size:\tthe size of the k-mer to run KMC with. The default is 10"
    s+= "\n-c --csv_path:\tthe path to the metadata .csv file with the genome_id, testing instrument, and Ct value of all genomes. There is no default for this option."
    s+= "\n-o --output_dir:\tthe directory for storing the k-mer DataFrame created by the script (will be created if it does not exist). The default is  ~/<genomes_dir>/output/"
    s+= "\n-d --df_name:\tthe name that the k-mer DataFrame created by the script will be stored as (must be a .csv file). The default is  'kmr_df.csv'"
    s+= "\n-i --dictionary_name:\tthe name that the dictionary { k-mer : column number } used to create the DataFrame will be stored as (must be a .pkl file). The default is 'kmr_dictionary.pkl'"
    return s


# concatenates all .fasta files in a directory into one file (to be used in runKMCConcat)
# parameters:
#    dir: the directory in which to concatenate all .fasta all_files
# returns: 'concat_file', the name of the concattenated file
def concatFiles(dir):
    os.chdir(dir)
    os.system("touch concat_file")

    # iterates through all .fasta files in the directory to add them all to the concat_file
    for filename in os.scandir(dir):
        if (filename.path.endswith(".fasta")):
            with open(filename) as f:
                file_contents = (f.read() + '\n')
            f.close()
            # writing the contents of the current .fasta file to the concat_file
            with open ('concat_file', 'a') as f:
                f.write(file_contents)
            f.close()

    return 'concat_file'



# runs KMC on the concattenated file with all the genomes
# this generates a list of all unique k-mers across all genomes
# parameters:
#    kmr_size: the size of the k-mer with which to run KMC
#    genomes_dir: the directory containing concat_file
#    concat_file_name: the name of the concattenated file with which to run KMC
# returns: the name of the file outputted by KMC containing all unique k-mers and their frequencies
def runKMCConcat(kmr_size, genomes_dir, concat_file_name):
    os.chdir(genomes_dir)
    cmd = './kmc.sh' + ' ' + str(kmr_size) + ' ' +  concat_file_name + ' concat_KMC ' + genomes_dir
    os.system(cmd)

    # Deleting the extra files produced by KMC
    os.system("rm concat_KMC.kmc_pre")
    os.system("rm concat_KMC.kmc_suf")

    #returning the name of the file containing all unique k-mers across all genomes in geomes_dir
    return "concat_KMC." + str(kmr_size) + ".kmrs"



# creates a dictionary of k-mer : column number to be used when filling in the DataFrame
# parameters:
#    kmr_size: the size of the k-mers used when running KMC
#    genomes_dir: the directory in which all_kmrs_file is
#    all_kmrs_file: the output of kmc on the concattenated file, contains all k-mers across all genomes
# returns: the dictionary of k-mer : column number for all k-mers
def createDictionary(kmr_size, genomes_dir, all_kmrs_file):
    kmr_dictionary = {}

    # iterrates through all lines all_kmrs_file (each line is a unique k-mer) and adds that k-mer to the dictionary
    # incriments col_num which keeps track of the column number to which the k-mer will correspond
    os.chdir(genomes_dir)
    file = open(all_kmrs_file)
    col_num = 0
    for aline in file:
        values = aline.split() # values: [k-mer, frequency]
        kmr_dictionary.update({values[0]: col_num}) #dictionary is k_mer : column number
        col_num+= 1
    file.close()

    return kmr_dictionary


# initializes a dataFrame with the correct number of columns and rows and initializes the values to 0
# creates a column for every k-mer and the instrument, genome_id, and ct value and a row for every genome in the directory
# parameters:
#    genomes_dir: the directory containing the genomes that will be added to the df
#    kmr_dictionary: the dictionary of k-mer : column index used to find the right columns to change the frequency of
# returns: the DataFrame initialized to 0s and transposed for more efficient access in memory
def initializeDf(genomes_dir, kmr_dictionary):
    # number of columns:
    N = len(kmr_dictionary) + 5 # creating 5 extra columns: 3 for instruments, 1 for genome_id and 1 for ct value
    # number of rows (number of genomes in the directory):
    M = 0
    for filename in os.scandir(genomes_dir):
        if (filename.path.endswith(".fasta")):
            M = M + 1

    df = pd.DataFrame(np.zeros((M, N), dtype=int)) # creating the DataFrame and initializing it to 0s
    # update first row (genome_id) to be strings:
    df[0] = df[0].astype(str)
    # transposing the DataFrame for more efficient memory access
    df = df.transpose()

    return df




# reads in the k-mers from the output files of KMC
# fills in the DataFrame with the frequencies of each k-mer and the Ct value and instrument of each genome
# parameters:
#    kmr_df: the DataFrame to add the k-mer frequencies to
#    kmr_dictionary: the dictionary of k-mer : column index used to find the right columns to change the frequency of
#    csv_path: the path to the file containing information on the genome_id, instrument, and Ct value of every genome
#    kmc_out_dir: the directory containing the outputs of KMC (the files with k-mer frequency)
#    kmr_size: the size of k-mers used
# returns: the updated DataFrame with k-mer frequencies, instrument, and Ct value filled in
def fillDf(kmr_df, kmr_dictionary, csv_path, kmc_out_dir, kmr_size):
    num_rows = len(kmr_df.index) # the number of rows (k-mers + instrument, Ct value, and genome_id)

    index = 0 # index counts the genomes that are being read in

    # iterates through all files in the kmc_out_dir and reads in the k-mers to the initialized DataFrmae
    os.chdir(kmc_out_dir)
    for filename in os.scandir(kmc_out_dir):
        if (filename.path.startswith(kmc_out_dir + "/MCoV-")): # checking if the file is in the correct format
            print("  Processing file:   ", filename, " (", index, ")")

            kmc_file = open(filename)

            # Getting the genome_id from the file name:
            genome_id = str(filename).replace("<DirEntry '", "")
            genome_id = genome_id.replace(("_kmc." + str(kmr_size) + ".kmrs'>"), "")

            # getting the Instrument and Ct value of the genome:
            genome_info = getInfo(csv_path, genome_id)

            if (genome_info == None): # the genome_id was not found in the metadata file
                # deleting the column of the DataFrame corresponding to the genome that was not found
                kmr_df = kmr_df.drop(kmr_df.columns[[0]], axis=1)
            else:
                kmr_df.at[0, index] = genome_id # adding the genome id to the first row of the current column

                # parsing the k-mer output file to get the frequency of every k-mer in the current genome and updating the DataFrmae:
                for aline in kmc_file:
                    values = aline.split() #values[0] = k-mer, values[1] = frequency

                    # dict_encoding is the row that corresponds to the current k-mer from the dictionary
                    dict_encoding = kmr_dictionary[values[0]] + 1 # adding 1 due to the genome_id in the first row
                    kmr_df.at[dict_encoding, index] = values[1] # updating the DataFrame with the frequency of the current k-mer

                ins = genome_info[0]
                ct = genome_info[1]
                # adding the Ct value of the current genome
                kmr_df.at[(num_rows - 1), index]  = ct # adding Ct value to the last row

                # adding instrument of the current genome:
                # getting the correct row of the genome:
                if (ins == "ALINITY"):
                    row = num_rows - 4
                elif (ins == "PANTHER"):
                    row = num_rows - 3
                elif (ins == "CEPHEID"):
                    row = num_rows - 2
                kmr_df.at[row, index] = 1 # updating the correct instrument row to be 1

            kmc_file.close()
            index = index + 1

    return kmr_df



# gets the instrument and ct value of the given genome_id from the csv file
# gets the Ct value and instrument of the a specific genome
# parameters:
#    csv_path: the path to the metadata csv file containing the information about the genome
#    genome_id: the genome id of the genome for which to find the Instrument and Ct value
# returns: a list of the instrument and Ct value or None if the genome_id was not in the csv file
def getInfo(csv_path, genome_id):
    csv_file = pd.read_csv(csv_path)
    for index, row in csv_file.iterrows():
        if (row["genome_id"] == genome_id):
            instrument =  row["INSTRUMENT"]
            ct = row["ct_value"]
            return [instrument, ct]
    return None # the genome_id was not found in the metadata file


# Transposes the DataFrame and renames the columns
# parameters: df: the DataFrame to transpose
# returns: the transposed DataFrame with column titles
def transposeDf(df):
    df = df.transpose()
    num_cols = len(df.columns)

    col_names = ["genome_id"] # creating a list of the column names
    kmr_index = 0
    for i in range(1, (num_cols) - 4): # leaving 4 columns at the end: 3 for the instruments and 1 for the Ct value
        col_names.append(kmr_index)
        kmr_index = kmr_index + 1

    # Appending the instruments and the Ct value
    col_names.append("alinity")
    col_names.append("panther")
    col_names.append("cepheid")
    col_names.append("ct_value")

    # updating the columns of the DataFrame:
    df.columns = col_names

    return df



# stores the DataFrame as a csv file
# parameters:
#    df: the DataFrame to be stpred
#    df_dir: the directory in which to store the DataFrame
#    df_name: the name to store the DataFrame as
def storeDataFrame(df, df_dir, df_name):
    if (os.path.isfile(df_dir) == False):
        os.system("mkdir " + df_dir) # creating the output directory if it does not already exist
    os.chdir(df_dir)
    df.to_csv(df_name)




# main function
# creates and stores a DataFrame with columns of every k-mer across all genomes, the instrument (one-hot encoded),
#  the Ct value (label), and the genome_id
# every genome is one row in the DataFrame
def main(argv):
    args = sys.argv
    # reads in parameters passed in by user through the command line or setting paramters to default values
    genomes_dir, kmc_out_dir, kmr_size, csv_path, output_dir, df_name, dictionary_name = parseParams(args)

    # concatenates all genomes in the directory into one file
    concat_file_name = concatFiles(genomes_dir)
    print("--createDataFrame.py-- concattenated genome files")

    # runs KMC on the concattenated file to generate a list of all k-mers across all genomes
    all_kmrs_file = runKMCConcat(kmr_size, genomes_dir, concat_file_name)
    print("--createDataFrame.py-- ran KMC on the concattenated file")

    # creates a dictionary of k-mer : column number to be used when creating the DataFrame
    kmr_dict = createDictionary(kmr_size, genomes_dir, all_kmrs_file)
    print("--createDataFrame.py-- created dictionary of k-mer : column number")

    # initializes the columns of the DataFrame to include every k-mer (using kmr_dict)
    init_df = initializeDf(genomes_dir, kmr_dict)
    print("--createDataFrame.py-- initialized DataFrame")

    # reads in the output files of KMC and adds a row to the DataFrame for every genome with the frequency of every k-mer (column)
    kmr_df = fillDf(init_df, kmr_dict, csv_path, kmc_out_dir, kmr_size)
    print("--createDataFrame.py-- filled in DataFrame with the frequency of every k-mer and the instrument and Ct value")

    # transposes the DataFrame and adds column titles
    kmr_df = transposeDf(kmr_df)
    print("--createDataFrame.py-- transposed the DataFrame and added column ttiles")

    # stores the DataFrame as a .csv file
    storeDataFrame(kmr_df, output_dir, df_name)

    # stores the dictionary as a .pkl file
    os.chdir(output_dir)
    pickle.dump(kmr_dict, open(dictionary_name, "wb"))
    print("--createDataFrame.py-- stored DataFrame as '", df_name, "' and dictionary as '", dictionary_name, "'  in  ", output_dir)



# if this is the script called by python, run main function
if __name__ == '__main__':
	main(sys.argv)
