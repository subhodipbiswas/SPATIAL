# SPATIAL

This is the GitHub repository corresponding to the algorithm **S**warm-based s**pa**atial meme**ti**c **al**gorithm (**SPATIAL**). Here, we apply SPATIAL to solve the problem of school boundary formation (also called school redistricting).

## Installation

The code is written in Python3.6 and the experiments were run on a machine using Ubuntu 18.04 LTS. You can follow the commands below for setting up your project.

### Setting up virtual environment
Assuming you have Python3, set up a virtual environment
```
pip install virtualenv
virtualenv -p python3 venv
```

### Activate the environment
Always make sure to activate it before running any python script of this project
```
source venv/bin/activate
```

## Package installation
Install the required packages contained in the file requirements.txt. This is a one-time thing, make sure the virtual environment is activated before performing this step.
```
pip install -r requirements.txt
```

Note that some geospatial packages in Python require dependencies like [GDAL](https://gdal.org/) to be already installed.

### Navigate to the 'src' folder containing the source scripts
```
cd ./src
```

## Run the code

You can simulate all the experiments using the following command:
```
make SPATIAL
```

OR


Simulate experiments for
1. Elementary school
```
./run_algo.py -s ES
```
2. Middle school
```
./run_algo.py -s MS
```
3. High school
```
./run_algo.py -s HS
```


### Deactivate the environment
Deactivate it before exiting the project
```
deactivate
```

### Data
The geospatial data used here is of [LCPS](https://www.lcps.org/) school district for the school year 2019-20. The data has been pre-processed for usage and may not accurately represent the policies of LCPS.

## Citation
If you use this data/code for your work, please consider citing the paper:
```
@inproceedings{10.1145/3397536.3422265,
author = {Biswas, Subhodip and Chen, Fanglan and Chen, Zhiqian and Lu, Chang-Tien and Ramakrishnan, Naren},
title = {Incorporating Domain Knowledge into Memetic Algorithms for Solving Spatial Optimization Problems},
year = {2020},
isbn = {9781450380195},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3397536.3422265},
doi = {10.1145/3397536.3422265},
booktitle = {Proceedings of the 28th International Conference on Advances in Geographic Information Systems},
pages = {25–35},
numpages = {11},
keywords = {Metaheuristic, domain knowledge, spatial optimization},
location = {Seattle, WA, USA},
series = {SIGSPATIAL '20}
}
```
## Help
Should you have queries, feel free to send an email to subhodip [at] cs.vt.edu
