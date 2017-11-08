#!/usr/bin/env bash

# Helper script to install dependencies for the RBNS pipeline

# See README.md for more information

cd "$(dirname "$0")"

conda_path=$(which conda)
os_name=$(uname)

if [ ! -f $conda_path ]; then
	echo "conda binary not found. Conda is required for RBNS_pipeline functionality. See the README file."
	exit 1
fi

if [ $os_name == "Linux" ]; then
	echo "Attempting to create shapeware environment..."
	conda env create -n rbns_pipeline --file conda/linux.yml

else
	echo "OS $os_name not supported."

fi


echo "Do you want to download the example files? (y/n)"

read dl_option

curl_path=$(command -v curl)
wget_path=$(command -v wget)

if [ $dl_option == 'y' ] || [ $dl_option == 'Y' ]; then

	if [ ! -f $curl_path ]; then
		echo "curl binary not found. curl or wget is required to download files."
		
	elif [ ! -f $wget_path ]; then
		echo "wget binary not found. curl or wget is required to download files."
		exit 1
	else
		example_data/download_files.sh
		echo "You should now be able to test the RBNS_pipeline with the following commands:"
		echo "source activate rbns_pipeline"
		echo "python src/RBNS_main.py test_data/settings.RBFOX3.json"

	fi
fi

