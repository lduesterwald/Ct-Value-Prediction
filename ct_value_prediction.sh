#!/bin/bash
# runs the Ct value prediction script pipeline

# sets the arguments to the appropriate variables
while getopts g:k:s:c:o:d:i:m:f:t:e:r:p:u:n: flag
do
    case "${flag}" in
        g) genomes_dir=${OPTARG};;
        k) kmc_out_dir=${OPTARG};;
        s) kmr_size=${OPTARG};;
        c) csv_path=${OPTARG};;
        o) output_dir=${OPTARG};;
        d) df_name=${OPTARG};;
        i) dictionary_name=${OPTARG};;
        m) model_name=${OPTARG};;
        f) output_file_name=${OPTARG};;
        t) num_features=${OPTARG};;
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
    echo "$space[-c csv_path - required] [-o output_dir] [-d df_name]"
    echo "$space[-i dictionary_name] [-m model_name] [-f output_file_name] "
    echo "$space[-t num_features] [-e test_size] [-r num_trees]"
    echo "$space[-p tree_depth] [-u row_subsampling] [-n genome_name - required]"
    echo " "
    echo "one or more required arguments missing: "
    echo "    -c: csv_path"
    echo "    -n: genome_name"
    exit 1
fi


# checking all the arguments and creating the commands for each script:
c1="python3 runKMC.py"
c2="python3 createDataFrame.py"
c3="python3 trainModel.py"
c4="python3 predictCt.py"

if ! [ -z "$genomes_dir" ]; then c1+=" -g $genomes_dir"; c2+=" -g $genomes_dir"; c4+=" -g $genomes_dir"; fi
if ! [ -z "$kmc_out_dir" ]; then c1+=" -k $kmc_out_dir"; c2+=" -k $kmc_out_dir"; fi
if ! [ -z "$kmr_size" ]; then c1+=" -s $kmr_size"; c4+=" -s $kmr_size"; fi
c2+=" -c $csv_path"; c4+=" -c $csv_path"
if ! [ -z "$output_dir" ]; then c2+=" -o $output_dir"; c3+=" -o $output_dir"; c4+=" -o $output_dir"; fi
if ! [ -z "$df_name" ]; then c2+=" -d $df_name"; c3+=" -d $df_name"; fi
if ! [ -z "$dictionary_name" ]; then c2+=" -i $dictionary_name"; c3+=" -i $dictionary_name"; c4+=" -i $dictionary_name"; fi
if ! [ -z "$model_name" ]; then c3+=" -m $model_name"; c4+=" -m $model_name"; fi
if ! [ -z "$output_file_name" ]; then c3+=" -f $output_file_name"; fi
if ! [ -z "$num_features" ]; then c3+=" -t $num_features"; fi
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

echo "stored all output as $output_file_name in $output_dir"
