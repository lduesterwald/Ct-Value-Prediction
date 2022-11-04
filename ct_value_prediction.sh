#!/bin/bash
# runs the Ct value prediction script pipeline

# sets the arguments to the appropriate variables
while getopts g:k:s:c:d:i:m:f:e:r:p:u:n: flag
do
    case "${flag}" in
        g) genomes_dir=${OPTARG};;
        k) kmc_out_dir=${OPTARG};;
        s) kmr_size=${OPTARG};;
        c) csv_path=${OPTARG};;
        d) df_name=${OPTARG};;
        i) dictionary_name=${OPTARG};;
        m) model_name=${OPTARG};;
        f) output_file_name=${OPTARG};;
        e) test_size=${OPTARG};;
        r) num_trees=${OPTARG};;
        p) tree_depth=${OPTARG};;
        u) row_sub=${OPTARG};;
        n) genome_name=${OPTARG};;
    esac
done


# checks if any required arguments are null
if [ -z "$csv_path" ] || [ -z "$genome_name" ]
then
    space="                              "
    echo usage: "$0 [-g genomes_dir] [-k kmc_out_dir] [-s kmr_size]"
    echo "$space[-c csv_path - required] [-d df_name] [-i dictionary_name]"
    echo "$space[-m model_name] [-f output_file_name] [-e test_size]"
    echo "$space[-r num_trees] [-p tree_depth] [-u row_subsampling]"
    echo "$space[-n genome_name - required]"
    echo " "
    echo "one or more required arguments missing: "
    echo "    -c: csv_path"
    echo "    -n: genome_name"
    exit 1
fi


# getting the directory where the bash script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# checking all the arguments and creating the commands for each script:
base="python3 "
base+="$SCRIPT_DIR"
base+="/"
c1="$base"
c1+="runKMC.py"
c2="$base"
c2+="createDataFrame.py"
c3="$base"
c3+="trainModel.py"
c4="$base"
c4+="predictCt.py"


if ! [ -z "$genomes_dir" ]; then c1+=" -g $genomes_dir"; c2+=" -g $genomes_dir"; c4+=" -g $genomes_dir"; fi
if ! [ -z "$kmc_out_dir" ]; then c1+=" -k $kmc_out_dir"; c2+=" -k $kmc_out_dir"; c4+=" -k $kmc_out_dir"; fi
if ! [ -z "$kmr_size" ]; then c1+=" -s $kmr_size"; c4+=" -s $kmr_size"; fi
c2+=" -c $csv_path"; c4+=" -c $csv_path"
if ! [ -z "$df_name" ]; then c2+=" -d $df_name"; c3+=" -d $df_name"; fi
if ! [ -z "$dictionary_name" ]; then c2+=" -i $dictionary_name"; c3+=" -i $dictionary_name"; c4+=" -i $dictionary_name"; fi
if ! [ -z "$model_name" ]; then c3+=" -m $model_name"; c4+=" -m $model_name"; fi
if ! [ -z "$output_file_name" ]; then c3+=" -f $output_file_name"; fi
if ! [ -z "$test_size" ]; then c3+=" -ts $test_size"; fi
if ! [ -z "$num_trees" ]; then c3+=" -nt $num_trees"; fi
if ! [ -z "$tree_depth" ]; then c3+=" -td $tree_depth"; fi
if ! [ -z "$row_sub" ]; then c3+=" -rs $row_sub"; fi
c4+=" -n $genome_name"


# running the scripts:
echo "$c1"
$c1
echo "$c2"
$c2
echo "$c3"
$c3
echo "$c4"
$c4

if [ -z "$output_file_name" ]
then
    echo "stored all output as output_file_trainModel"
else
    echo "stored all output as $output_file_name"
fi

