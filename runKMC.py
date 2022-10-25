import sys
import os

# this function parses any parameters passed in through the command line or sets them to a default value
# parameters:
#    args: the list of arguments passed in through the command line
# returns: genomes_dir, kmc_out_dir, kmr_size
def parseParams(args):
    # setting default values for parameters:
    genomes_dir = "./" # (-g) the directory with genomes as .fasta files
    kmc_out_dir = genomes_dir + "kmc_output/" # (-k) the directory for storing the output files of KMC with k-mer counts
    kmr_size = 10 # (-s) the size k-mer to run KMC with

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
            kmc_out_dir = genomes_dir + "kmc_output"
        elif (args[i] == "-k" or args[i] == "--kmc_out_dir"):
            kmc_out_dir = args[i + 1]
            kmc_out_dir = os.path.abspath(kmc_out_dir) + "/"
        elif (args[i] == "-s" or args[i] == "--kmr_size"):
            kmr_size = args[i + 1]


    return genomes_dir, kmc_out_dir, kmr_size


# returns a string of all the options for the script if the script was called with -h or --help
def helpOption():
    s = "-g --genomes_dir:\tthe directory containing the fasta files to train the model, the KMC executable package, and the kmc.sh script. The default is './'"
    s+= "\n-k --kmc_out_dir:\tthe directory for storing the output files created by KMC. The default is genomes_dir/kmc_out_dir"
    s+= "\n-s --kmr_size:\tthe size of the k-mer to run KMC with. The default is 10"
    return s


# runs KMC on all .fasta files in a directory and puts the outputs into an output sub-directory
# parameters:
#    genomes_dir: the directory containing the genomes as .fasta all_files
#    kmc_out_dir: the directory to which the output of KMC will be stored
#    kmr_size: the size of the k-mers with which to run KMC
def runKMC(genomes_dir, kmc_out_dir, kmr_size):
    os.chdir(genomes_dir)
    os.system("mkdir " + kmc_out_dir) # Create sub directory for kmc_out_dir

    # iterate over all files in genomes_dir
    #  for each .fasta file, runs kmc.sh (runs kmc and kmc_dump)
    for filename in os.scandir(genomes_dir):
        if (filename.path.endswith(".fasta")):
            mcov = str(filename).strip("<DirEntry '.fasta'>")
            input_file = mcov + ".fasta"
            out_file = mcov + "_kmc"

            cmd = 'kmc.sh' + ' ' + str(kmr_size) + ' ' + input_file + ' ' + out_file + ' ' + kmc_out_dir
            os.system(cmd) #runs kmc.sh for the current file



# deletes all temporary files that end in .suf or .pre (generated by KMC)
# moves the output files from KMC into the output directory
# parameters:
#    genomes_dir: the directory containing the genomes as .fasta all_files
#    kmc_out_dir: the directory to which the output of KMC will be stored
def cleanDir(genomes_dir, kmc_out_dir):
    for filename in os.scandir(genomes_dir):
        if (filename.path.endswith("suf") or filename.path.endswith("pre")):
            os.system("rm " + filename.path) #removes .pre and .suf files from KMC
        elif (filename.path.endswith("kmrs")):
            cmd = "mv " + filename.path + " " + kmc_out_dir
            os.system(cmd)



# main function
# generates k-mer lists for every genome in a directory by running KMC
def main(argv):
    args = sys.argv
    # reads in parameters passed in by user through the command line or setting paramters to default values
    genomes_dir, kmc_out_dir, kmr_size = parseParams(args)

    #runs kmc.sh (kmc and kmc_dump) for every file in genomes_dir directory
    runKMC(genomes_dir, kmc_out_dir, kmr_size)
    print("--runKMC.py-- finished running KMC")

    #calls cleanDir which moves all files to an output directory and deletes extra kmc files .suf and .pre
    cleanDir(genomes_dir, kmc_out_dir)



# if this is the script called by python, run main function
if __name__ == '__main__':
	main(sys.argv)
