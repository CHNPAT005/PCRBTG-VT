# The Epps effect under alternative sampling schemes

## Authors
* Patrick Chang
* Etienne Pienaar
* Tim Gebbie

## Link to resources:

Dataset DOI: [10.25375/uct.12315092](https://zivahub.uct.ac.za/articles/Using_the_Epps_effect_to_detect_discrete_data_generating_processes_Dataset/12315092/1)

## Steps for Replication:
- Change directories for all the files under [/Scripts](https://github.com/CHNPAT005/PCRBTG-VT/tree/main/Scripts). Currently the directories are set as: `cd("/Users/patrickchang1/PCEPTG-VT")`. Change this to where you have stored the repository. 

	- Run [/Scripts/Simulation.jl](https://github.com/CHNPAT005/PCRBTG-VT/blob/main/Scripts/Simulation.jl) and [/Scripts/MoreVT.jl](https://github.com/CHNPAT005/PCRBTG-VT/blob/main/Scripts/MoreVT.jl) to reproduce the simulation work
	
 - To reproduce the Empirical analysis - download the processed dataset from ZivaHub and put the csv files into the folder `/Real Data`.
 	- Run the scripts that start with `Emp` to reproduce the empirical work

- We have included the plots under `/Plots` and Computed results under `/Computed Data` if one does not wish to re-run everything.

## Using the functions for other purposes:

We have built a suite of useful functions that can be used for other applications. For a demonstration on how to use them, please visit: https://github.com/CHNPAT005/PCEPTG-MSC.