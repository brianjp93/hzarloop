# Installation for ACISS Users
The hzarloop and clusterfit programs can be run on the University of Oregon's computing cluster using as many nodes as you can get your hands on.  These installation instructions are for installing on ACISS specifically.  However, if you have access to another computer cluster, these instructions may have close analogues for how to install for your cluster.

### Install hzar for R
To do this, you must first be ssh'd into ACISS.

`$ ssh username@aciss.uoregon.edu`

Then load up R.  I used 3.2.3, but later versions will probably work just as well. (maybe stick to 3.2.3 if you want to be safe)  Install hzar in an interactive R session.  I had to use the http servers because ACISS didn't let me download from the https servers.

`$ module load R/3.2.3`

We also need to make a directory for HZAR to be installed.

`$ mkdir ~/Rpackages`

Start up an R session

`$ R`

Then once you are in an interactive R session, type...

`> install.packages(hzar, lib="~/Rpackages")`

Use any of the http servers to install hzar.

Then we need to create a .Renviron file in the home directory so that R knows where to look for the hzar module.

`$ touch ~/.Renviron`

edit the file and add the line

`R_LIBS=~/Rpackages`

