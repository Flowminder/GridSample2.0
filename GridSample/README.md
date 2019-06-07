# GRID SAMPLE SET-UP INSTRUCTIONS 

## Introduction

The GridSample2.0 code is a mixture of Python and R ( for the gridEZ algorithm).

The code obtained from the GitHub repository contains a demo called `run_demo.py` and the files required to run it (population rasters,shapefiles, etc), other than two datasets which are too large to be hosted on GitHub:

- **WorldPop dataset**:  https://www.worldpop.org/geodata/summary?id=85
- **GHS-SMOD**: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_POP_GLOBE_R2016A/

Users should download these datasets and add them to the `GridSample2.0/GridSample/SampleBuilder/test_files/` directory.

It is provided "dockerized" to help set up the required environment easily, though the user may wish to set an alternative one up themselves.

## Build using docker

You will require some knowledge of docker, and to have docker installed (see www.docker.com).

The instructions given here assume the GIT repository is pulled into a dir called `GridSample_standalone`.

Go to GridSample_standalone

Change the 3 paths in the `docker-compose.yml` file to match the locations where you have put the `src` dir, the `test_files` dir, and where you wish the results to be located.

Build the docker container:

```
	docker-compose build
```
	
Start the container:

```
	docker-compose up -d
```
	
Check it is running (Look for a container with a name like gridsample_grid_sample_standalone_1 when listing all running containers):

```
	docker ps
```	
	
## Run the demo app

Enter the docker container in interactive mode (replacing gridsample_grid_sample_standalone_1 in the below code with the name of your container):

```
	docker exec -it gridsample_grid_sample_standalone_1 bash
```
	
Run the demo:

```
	python /src/run_demo.py
```
	
	
## Other notes

The two main source files will print out debug messages if required.

In sample_builder_utils.python and/or Mainfx_phase2_v1.py set:


```
	enable_debug = True
```
	
There is also an option to keep a number of intermediate files created by the algorithm which might be useful.
In  Mainfx_phase2_v1.py set:


```
	delete_intermediate_outputs = False
```
	






