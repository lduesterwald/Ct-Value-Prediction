import sys
import os
import pandas as pd
import numpy as np
import pickle

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

# For the model evaulation:
import math
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error


# this function parses any parameters passed in through the command line or sets them to a default value
# parameters:
#    args: the list of arguments passed in through the command line
# returns: output_dir, df_name, dictionary_name, model_name, output_file_name, num_features, test_size, nt, td, rs
def parseParams(args):
    # setting default values for parameters:
    output_dir = "./output/" # (-o) the directory where dataFrame is stored and where the model will be stored
    df_name = "kmr_df.csv" # (-d) the name that of the k-mer dataFrame (in output_dir)
    model_name = "ct_model.sav" # (-m) the name to store the Ct value prediction model as (in output_dir)
    output_file_name = "output_file_trainModel" # (-f) the name of the file to which to write the results
    test_size = 0.2 # (-ts) the test size for the train test split of the data
    # parameters for the model, default values are the optimal ones found during prameter tuning
    num_trees = 400 # (-nt) number of trees
    tree_depth = None # (-td) tree depth
    row_subsampling = 0.25 # (-rs) row subsampling


    # parsing any parameters passed in through the command line
    for i in range(len(args)):
        if(args[i] == "-h" or args[i] == "--help"):
            print(helpOption())
            sys.exit()
        if (i == len(args) - 1):
            break
        elif (args[i] == "-o" or args[i] == "--output_dir"):
            output_dir = args[i + 1]
            output_dir = os.path.abspath(output_dir) + "/"
        elif (args[i] == "-d" or args[i] == "--df_name"):
            df_name = args[i + 1]
        elif (args[i] == "-m" or args[i] == "--model_name"):
            model_name = args[i + 1]
        elif (args[i] == "-f" or args[i] == "--output_file_name"):
            output_file_name = args[i + 1]
        elif( args[i] == "-ts" or args[i] == "--test_size"):
            test_size = args[i + 1]
        elif( args[i] == "-nt" or args[i] == "--num_trees"):
            num_trees = args[i + 1]
        elif( args[i] == "-td" or args[i] == "--tree_depth"):
            td = args[i + 1]
        elif( args[i] == "-rs" or args[i] == "--row_subsampling"):
            rs = args[i + 1]

    return output_dir, df_name, model_name, output_file_name, test_size, num_trees, tree_depth, row_subsampling




# returns a string of all the options for the script if the script was called with -h or --help
def helpOption():
    s = "-o --output_dir:\tthe directory where dataFrame is stored and where the model will be stored. The default is './output/'"
    s+= "\n-d --df_name:\tthe name of the k-mer DataFrame to train the model on (must be a .csv file). The default is 'kmr_df.csv'"
    s+= "\n-m --model_name:\tthe name to store the Ct value prediction model created by the script as (must be a .sav file). The default is 'ct_model.sav'"
    s+= "\n-f --output_file_name:\tthe name of the output file for this script. This file will be created in the <output_dir> directory. The default is 'output_file_trainModel'"
    s+= "\n-ts --test_size:\tthe size of the test set to be used in the train_test_split during model training and evaluation. The default is 0.2"
    s+="\n-nt --num_trees:\tthe 'n_estimators' (number of trees) parameter in the Random Forest regressor. The default is 400"
    s+="\n-td --tree_depth:\tthe 'max_depth' (tree depth) parameter in the Random Forest regressor. The default is None"
    s+="\n-rs --row_subsampling:\tthe 'max_samples' parameter in the Random Forest regressor. the default is 0.25"
    return s


# opens the DataFrame (stored as a .csv file)
# parameters:
#    df_dir: the directory in which the DataFrame is stored
#    df_name: the name that the DataFrame is stored as
# returns: the DataFrame
def openDataFrame(df_dir, df_name, start_dir):
    os.chdir(df_dir)
    df_path = df_dir + df_name
    df = pd.read_csv(df_path)
    df = df.dropna()   #dropping any Nans
    os.chdir(start_dir)

    return df


# splits the DataFrame into a train set and a test set and converts them to numpy array
# parameter:
#    df: the DataFrame to be split into a train and test set
#    ts: the test set size (as a decimal)
# returns: the train set, train labels, test set, and test labels
def splitDf(df, ts):
    mat = pd.DataFrame(df).to_numpy()
    train_set, test_set = train_test_split(mat, test_size=ts, random_state=42)

    r, c = train_set.shape
    train_labels = train_set[:,(c-1)].copy() # copying the Ct values into train_labels
    #dropping the columns for the ct_value (c-1), index (1), and genome_id (0)
    train_set = np.delete(train_set, [(c - 1), 1, 0], 1)

    r, c = test_set.shape
    test_labels = test_set[:,(c - 1)].copy() # copying the Ct values into test_labels
    #dropping the columns for the ct_value (c-1), index (1), and genome_id (0)
    test_set = np.delete(test_set, [(c - 1), 1, 0], 1)

    return train_set, train_labels, test_set, test_labels


# evaluates the model on the test set
# calculates the R2 and RMSE using the model's predictions on the test set
# calls getAccuracyWithinIntervals to compute the model's accuracy within a set of intervals
# writes the results to a file
# parameters:
#    model: the fitted  model to evaluate
#    test set: the test set to predict
#    test_labels: the true values for the test set
#    output_file_path: the path to the file to which to write the results
def evaluateModel(model, test_set, test_labels, output_file_path):
    # predicting the test set:
    predictions = model.predict(test_set)
    print("--trainMode.py-- predicted test set")

    mse = mean_squared_error(test_labels, predictions)
    rmse = math.sqrt(mse)
    r2 = r2_score(test_labels, predictions)

    #Writing the R2 and RMSE to the output file
    f = open((output_file_path), "a")
    s = "\nModel Accuracy: \nR2: " + str(r2) + "\nRMSE: " + str(rmse)
    f.write(s)
    f.close()

    # calculating the accuracy within intervals
    getAccuracyWithinIntervals(predictions, test_labels, output_file_path)


# computes model's accuracy within intervals as % of times prediction was within interval of actual and writes the results to a file
# parameters:
#    predictions and true_values: the model's predictions and the true_values it was predicting
#    output_file_path: the path to the output file to which to write the results
def getAccuracyWithinIntervals(predictions, true_values, output_file_path):
    intervals = [6, 5, 4, 3, 2, 1] # the intervals within which to calculate the model' accuracy

    # accuracies will keep track of the total times that the model's prediction was within the interval of the true value
    accuracies = []
    for i in range(len(intervals)):
        accuracies.append(0)

    # iterating through all the predictions and actuals and updates the accuracies list
    for i in range(len(predictions)):
        true = true_values.item(i)
        predicted = predictions.item(i)
        # iterates through all the intervals to calculate whether the prediction was within that interval of the true value
        for j in range(len(intervals)):
           intv = intervals[j]
           if (abs(true - predicted) <= intv): #checking if the true was within the interval of the prediction
               accuracies[j] = accuracies[j] + 1


    f = open((output_file_path), "a")
    f.write("\n\nAccuracy Within Intervals:\n")

    for i in range(len(accuracies)):
        # dividing the number of times the prediction was within the interval of the true value by the total number of predictions
        accuracies[i] = accuracies[i] / len(predictions)
        #Writing the accuracy of the model within the interval to the output file
        f.write("within " + str(intervals[i]) + ":  " + str(accuracies[i]) + "\n")

    f.close()


# main function
# trains a model on the k-mer dataFrame created by createDataFrame.py
# evaluates the model's accuracy (r2 score, RMSE, accuracy within intervals)
# prints a list of the top features (k-mers) of the model
def main(argv):
    # current working directory:
    start_dir = os.getcwd()

    args = sys.argv
    # reads in parameters passed in by user through the command line or setting paramters to default values
    output_dir, df_name, model_name, output_file_name, test_size, num_trees, tree_depth, row_subsampling = parseParams(args)


    # opening the DataFrame
    df = openDataFrame(output_dir, df_name, start_dir)
    print("--trainModel.py-- opened DataFrame")

    #splitting the DataFrame into a train and test set
    train_set, train_labels, test_set, test_labels = splitDf(df, test_size)
    print("--trainModel.py-- split the DataFrame into train and test sets")

    model =  RandomForestRegressor(n_estimators=num_trees, max_depth=tree_depth, random_state=42, max_samples=row_subsampling)

    # training the model:
    model.fit(train_set, train_labels)
    print("--trainModel.py-- trained the model")

    # creating the output file:
    os.chdir(output_dir)
    os.system("touch " + output_file_name)
    os.chdir(start_dir)

    # evaluating the model's accuracy and writing the results to a file
    evaluateModel(model, test_set, test_labels, (output_dir + output_file_name))


    # storing the model to the output_dir
    os.chdir(output_dir)
    pickle.dump(model, open(model_name, "wb"))
    print("--trainModel.py-- stored results as '", output_file_name, "' in  '", output_dir, "'")


# if this is the script called by python, run main function
if __name__ == '__main__':
	main(sys.argv)

