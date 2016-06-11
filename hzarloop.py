"""
hzarloop.py

April 25, 2016
- Being developed in python 3 on ubuntu 16.04
- A substantial edit to Madeline's original script "hzar_loop.py"
    which also fit clines to loci data.  Runs an R script to
    make the fits.
"""

import os
import subprocess
import csv
import time
import sys

__author__ = "Brian Perrett"


class Hzarloop:

    # include modelIII
    rawrscript = """args = commandArgs(trailingOnly=TRUE)
        lociname = args[1]
        num = args[2]
        csvstring = args[3]

        loci<-read.csv(text='{csvtext}', header=T, sep=',')
        library(hzar)
        # library(hzar, lib.loc="/home7/perrett/Rpackages")
        # Set chain length. default = 1e5
        chainLength={chainlength}

        # set seeds for dif models
        # make random?
        # mainSeed= list(A=c(596,528,124,978,544,99),
        #     B=c(528,124,978,544,99,596),
        #     C=c(124,978,544,99,596,528))
        
        mainSeed= list(A=c(sample(1:1000, 6)),
            B=c(sample(1:1000, 6)),
            C=c(sample(1:1000, 6)))


        # execute parallel mode
        # if(require(doMC)){{
          # registerDoMC()
          # }} else {{
          registerDoSEQ();
          # }}

        if(length(apropos('^mim$',ignore.case=FALSE)) == 0 ||
            !is.list(mim) ) mim <- list()
        mim[[lociname]] <- list();
        ##Space to hold the observed data
        mim[[lociname]]$obs <- list();
        ##Space to hold the models to fit
        mim[[lociname]]$models <- list();
        ##Space to hold the compiled fit requests
        mim[[lociname]]$fitRs <- list();
        ##Space to hold the output data chains
        mim[[lociname]]$runs <- list();
        ##Space to hold the analysed data
        mim[[lociname]]$analysis <- list(); 

        #Load data
        print("distance")
        print(loci$distance)
        print("lociname")
        print(loci[[lociname]])
        print("num")
        print(loci[[num]])
        mim[[lociname]]$obs <- hzar.doMolecularData1DPops(loci$distance, loci[[lociname]],loci[[num]]);

        #Make function to load models
        mim.loadmodel <- function(scaling,tails, id=paste(scaling,tails,sep='.'))
        mim[[lociname]]$models[[id]] <<- hzar.makeCline1DFreq(mim[[lociname]]$obs, scaling, tails)
        #load sigmoid model
        mim.loadmodel('free' ,'none','modelII');
        #load three-parts model
        mim.loadmodel('free' ,'both','modelIII');
        
        # constrain distances to -30 to 50
        mim[[lociname]]$models <- sapply(mim[[lociname]]$models, hzar.model.addBoxReq, -30, 50, simplify=FALSE)
        mim[[lociname]]$models <- sapply(mim[[lociname]]$models, hzar.model.addCenterRange, -1, 1, simplify=FALSE)


        #Compile models to prepare for fitting
        # Brian: Maybe change the next line to include distance restrictions (-30, 50) ?? 
        mim[[lociname]]$fitRs$init <- sapply(mim[[lociname]]$models, hzar.first.fitRequest.old.ML, obsData=mim[[lociname]]$obs, verbose=FALSE, simplify=FALSE)
        # Define chainlength, burnin, and seed
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$chainLength <- chainLength;
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$burnin <- chainLength %/% 10;
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$seed[[1]] <- mainSeed$B

        mim[[lociname]]$fitRs$init$modelIII$mcmcParam$chainLength <- chainLength;
        mim[[lociname]]$fitRs$init$modelIII$mcmcParam$burnin <- chainLength %/% 10;
        mim[[lociname]]$fitRs$init$modelIII$mcmcParam$seed[[1]] <- mainSeed$B 

        #Run model for initial chain
        mim[[lociname]]$runs$init <- list()
        mim[[lociname]]$runs$init$modelII <-hzar.doFit(mim[[lociname]]$fitRs$init$modelII)
        mim[[lociname]]$runs$init$modelIII <- hzar.doFit(mim[[lociname]]$fitRs$init$modelIII)

        ##Compile a new set of fit requests using the initial chains
        mim[[lociname]]$fitRs$chains <- lapply(mim[[lociname]]$runs$init, hzar.next.fitRequest)
        ## Replicate each fit request 3 times, keeping the original
        ## seeds while switching to a new seed channel.
        mim[[lociname]]$fitRs$chains <- hzar.multiFitRequest(mim[[lociname]]$fitRs$chains,each=3,baseSeed=NULL)

        ## runif(9,-30,600) center for modelII, modelIII    
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['center']= 0

        ## runif(9,0,630) width for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['width']= 20

        ## runif(6,0,1) pMin for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['pMin']= 0.2
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['pMin']= 0.2 
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['pMin']= 0.2 
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['pMin']= 0.2 
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['pMin']= 0.2 
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['pMin']= 0.2

        ## runif(6,0,1) pMax for modelII, modelIII 
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['pMax']= 0.8
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['pMax']= 0.8
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['pMax']= 0.8 
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['pMax']=  0.8 
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['pMax']=  0.8 
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['pMax']=  0.8

        ## runif(3,0,630) deltaL for modelIII
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['deltaL']= 20
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['deltaL']= 20 
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['deltaL']= 20

        ## runif(3,0,1) tauL for modelIII 
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['tauL']= runif(1, 0.0, 1.0)
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['tauL']= runif(1, 0.0, 1.0)
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['tauL']= runif(1, 0.0, 1.0)

        ## runif(3,0,630) deltaR for modelIII
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['deltaR']=60 
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['deltaR']=60 
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['deltaR']=60

        ## runif(3,0,1) tauR for modelIII
        mim[[lociname]]$fitRs$chains[[4]]$modelParam$init['tauR']= runif(1, 0.0, 1.0)
        mim[[lociname]]$fitRs$chains[[5]]$modelParam$init['tauR']= runif(1, 0.0, 1.0)
        mim[[lociname]]$fitRs$chains[[6]]$modelParam$init['tauR']= runif(1, 0.0, 1.0)

        ## Run a chain of 3 runs for every fit request
        mim[[lociname]]$runs$chains <-  hzar.doChain.multi(mim[[lociname]]$fitRs$chains, doPar=TRUE, inOrder=FALSE,count=2)

        # Create null model
        mim[[lociname]]$analysis$initDGs <- list( nullModel = hzar.dataGroup.null(mim[[lociname]]$obs))

        ## Create a model data group (hzar.dataGroup object) for each
        ## model from the initial runs.
        mim[[lociname]]$analysis$initDGs$modelII <- hzar.dataGroup.add(mim[[lociname]]$runs$init$modelII)
        mim[[lociname]]$analysis$initDGs$modelIII <- hzar.dataGroup.add(mim[[lociname]]$runs$init$modelIII)


        ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
        ## just created, copying the naming scheme (nullModel,
        ## modelII, modelIII).
        mim[[lociname]]$analysis$oDG <- hzar.make.obsDataGroup(mim[[lociname]]$analysis$initDGs)
        mim[[lociname]]$analysis$oDG <- hzar.copyModelLabels(mim[[lociname]]$analysis$initDGs, mim[[lociname]]$analysis$oDG)

        ## Convert all runs to hzar.dataGroup objects, adding them to
        ## the hzar.obsDataGroup object.
        mim[[lociname]]$analysis$oDG <- hzar.make.obsDataGroup(lapply(mim[[lociname]]$runs$chains,hzar.dataGroup.add),mim[[lociname]]$analysis$oDG);

        ## Do model selection based on the AICc scores
        print(mim[[lociname]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(mim[[lociname]]$analysis$oDG));

        ## Print out the model with the minimum AICc score
        print(mim[[lociname]]$analysis$model.name <-rownames(mim[[lociname]]$analysis$AICcTable)[[ which.min(mim[[lociname]]$analysis$AICcTable$AICc )]])

        ## Extract the hzar.dataGroup object for the selected model
        mim[[lociname]]$analysis$model.selected <-mim[[lociname]]$analysis$oDG$data.groups[[mim[[lociname]]$analysis$model.name]]

        ## Look at the variation in parameters for the selected model
        # print(hzar.getLLCutParam(mim[[lociname]]$analysis$model.selected, names(mim[[lociname]]$analysis$model.selected$data.param)));
        if (mim[[lociname]]$analysis$model.name == "nullModel") {{
            print("hello")
        }} else {{
            print(hzar.getLLCutParam(mim[[lociname]]$analysis$model.selected, names(mim[[lociname]]$analysis$model.selected$data.param)))
        }}
        

        ## Print the maximum likelihood cline for the selected model
        print(hzar.get.ML.cline(mim[[lociname]]$analysis$model.selected))

        # hzar.plot.cline(mim[[lociname]]$analysis$model.selected);"""  # It's a RAWR script. Hah. But actually it's a raw R script.

    # no modelIII
    rawrscript1 = """args = commandArgs(trailingOnly=TRUE)
        lociname = args[1]
        num = args[2]
        csvstring = args[3]

        loci<-read.csv(text='{csvtext}', header=T, sep=',')
        # library(hzar)
        library(hzar, lib.loc="/home7/perrett/Rpackages")
        # Set chain length. default = 1e5
        chainLength={chainlength}

        # set seeds for dif models
        # make random?
        # mainSeed= list(A=c(596,528,124,978,544,99),
        #     B=c(528,124,978,544,99,596),
        #     C=c(124,978,544,99,596,528))

        mainSeed= list(A=c(sample(1:1000, 6)),
            B=c(sample(1:1000, 6)),
            C=c(sample(1:1000, 6)))


        # execute parallel mode
        # if(require(doMC)){{
          # registerDoMC()
          # }} else {{
          registerDoSEQ();
          # }}

        if(length(apropos('^mim$',ignore.case=FALSE)) == 0 ||
            !is.list(mim) ) mim <- list()
        mim[[lociname]] <- list();
        ##Space to hold the observed data
        mim[[lociname]]$obs <- list();
        ##Space to hold the models to fit
        mim[[lociname]]$models <- list();
        ##Space to hold the compiled fit requests
        mim[[lociname]]$fitRs <- list();
        ##Space to hold the output data chains
        mim[[lociname]]$runs <- list();
        ##Space to hold the analysed data
        mim[[lociname]]$analysis <- list();

        #Load data
        print("distance")
        print(loci$distance)
        print("lociname")
        print(loci[[lociname]])
        print("num")
        print(loci[[num]])
        mim[[lociname]]$obs <- hzar.doMolecularData1DPops(loci$distance, loci[[lociname]],loci[[num]]);

        #Make function to load models
        mim.loadmodel <- function(scaling,tails, id=paste(scaling,tails,sep='.'))
        mim[[lociname]]$models[[id]] <<- hzar.makeCline1DFreq(mim[[lociname]]$obs, scaling, tails)
        #load sigmoid model
        mim.loadmodel('free' ,'none','modelII');

        # constrain distances to -30 to 50
        mim[[lociname]]$models <- sapply(mim[[lociname]]$models, hzar.model.addBoxReq, -30, 50, simplify=FALSE)
        mim[[lociname]]$models <- sapply(mim[[lociname]]$models, hzar.model.addCenterRange, -1, 1, simplify=FALSE)


        #Compile models to prepare for fitting
        # Brian: Maybe change the next line to include distance restrictions (-30, 50) ??
        mim[[lociname]]$fitRs$init <- sapply(mim[[lociname]]$models, hzar.first.fitRequest.old.ML, obsData=mim[[lociname]]$obs, verbose=FALSE, simplify=FALSE)
        # Define chainlength, burnin, and seed
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$chainLength <- chainLength;
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$burnin <- chainLength %/% 10;
        mim[[lociname]]$fitRs$init$modelII$mcmcParam$seed[[1]] <- mainSeed$B

        #Run model for initial chain
        mim[[lociname]]$runs$init <- list()
        mim[[lociname]]$runs$init$modelII <-hzar.doFit(mim[[lociname]]$fitRs$init$modelII)

        ##Compile a new set of fit requests using the initial chains
        mim[[lociname]]$fitRs$chains <- lapply(mim[[lociname]]$runs$init, hzar.next.fitRequest)
        ## Replicate each fit request 3 times, keeping the original
        ## seeds while switching to a new seed channel.
        mim[[lociname]]$fitRs$chains <- hzar.multiFitRequest(mim[[lociname]]$fitRs$chains,each=3,baseSeed=NULL)

        ## runif(9,-30,600) center for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['center']= 0
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['center']= 0

        ## runif(9,0,630) width for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['width']= 20
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['width']= 20

        ## runif(6,0,1) pMin for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['pMin']= 0.2
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['pMin']= 0.2
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['pMin']= 0.2

        ## runif(6,0,1) pMax for modelII, modelIII
        mim[[lociname]]$fitRs$chains[[1]]$modelParam$init['pMax']= 0.8
        mim[[lociname]]$fitRs$chains[[2]]$modelParam$init['pMax']= 0.8
        mim[[lociname]]$fitRs$chains[[3]]$modelParam$init['pMax']= 0.8

        ## Run a chain of 3 runs for every fit request
        mim[[lociname]]$runs$chains <-  hzar.doChain.multi(mim[[lociname]]$fitRs$chains, doPar=TRUE, inOrder=FALSE,count=2)

        # Create null model
        mim[[lociname]]$analysis$initDGs <- list( nullModel = hzar.dataGroup.null(mim[[lociname]]$obs))

        ## Create a model data group (hzar.dataGroup object) for each
        ## model from the initial runs.
        mim[[lociname]]$analysis$initDGs$modelII <- hzar.dataGroup.add(mim[[lociname]]$runs$init$modelII)


        ## Create a hzar.obsDataGroup object from the four hzar.dataGroup
        ## just created, copying the naming scheme (nullModel,
        ## modelII, modelIII).
        mim[[lociname]]$analysis$oDG <- hzar.make.obsDataGroup(mim[[lociname]]$analysis$initDGs)
        mim[[lociname]]$analysis$oDG <- hzar.copyModelLabels(mim[[lociname]]$analysis$initDGs, mim[[lociname]]$analysis$oDG)

        ## Convert all runs to hzar.dataGroup objects, adding them to
        ## the hzar.obsDataGroup object.
        mim[[lociname]]$analysis$oDG <- hzar.make.obsDataGroup(lapply(mim[[lociname]]$runs$chains,hzar.dataGroup.add),mim[[lociname]]$analysis$oDG);

        ## Do model selection based on the AICc scores
        print(mim[[lociname]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(mim[[lociname]]$analysis$oDG));

        ## Print out the model with the minimum AICc score
        print(mim[[lociname]]$analysis$model.name <-rownames(mim[[lociname]]$analysis$AICcTable)[[ which.min(mim[[lociname]]$analysis$AICcTable$AICc )]])

        ## Extract the hzar.dataGroup object for the selected model
        mim[[lociname]]$analysis$model.selected <-mim[[lociname]]$analysis$oDG$data.groups[[mim[[lociname]]$analysis$model.name]]

        ## Look at the variation in parameters for the selected model
        # print(hzar.getLLCutParam(mim[[lociname]]$analysis$model.selected, names(mim[[lociname]]$analysis$model.selected$data.param)));
        if (mim[[lociname]]$analysis$model.name == "nullModel") {{
            print("hello")
        }} else {{
            print(hzar.getLLCutParam(mim[[lociname]]$analysis$model.selected, names(mim[[lociname]]$analysis$model.selected$data.param)))
        }}

        ## Print the maximum likelihood cline for the selected model
        print(hzar.get.ML.cline(mim[[lociname]]$analysis$model.selected))

        # hzar.plot.cline(mim[[lociname]]$analysis$model.selected);"""  # It's a RAWR script. Hah. But actually it's a raw R script.

    headers = [
        "ID",
        "ModelII_AICc",
        "ModelIII_AICc",
        "nullModel_AICc",
        "best_model",
        "best_model_loglike",
        "center",
        "center2LLLow",
        "center2LLHigh",
        "width",
        "width2LLLow",
        "width2LLHigh",
        "slope",
        "slope2LLLow",
        "slope2LLHigh",
        "Pmin",
        "Pmin2LLLow",
        "Pmin2LLHigh",
        "Pmax",
        "Pmax2LLLow",
        "Pmax2LLHigh",
        "Pmax-Pmin",
        "Pmax-PminLLLow",
        "Pmax-PminLLHigh",
        "deltaL",
        "deltaL2LLLow",
        "deltaL2LLHigh",
        "tauL",
        "tauL2LLLow",
        "tauL2LLHigh",
        "deltaR",
        "deltaR2LLLow",
        "deltaR2LLHigh",
        "tauR",
        "tauR2LLLow",
        "tauR2LLHigh",
        "lg",
        "read_pos",
        "scaffold_id",
        "lg_position"
    ]
    conversion = {
                "pMin2LLLow": "Pmin2LLLow",
                "width2LLLow": "width2LLLow",
                "pMax2LLHigh": "Pmax2LLHigh",
                "center": "center",
                "pMax": "Pmax",
                "pMin2LLHigh": "Pmin2LLHigh",
                "pMin": "Pmin",
                "center2LLHigh": "center2LLHigh",
                "width": "width",
                "center2LLLow": "center2LLLow",
                "width2LLHigh": "width2LLHigh",
                "logLike": "best_model_loglike",
                "pMax2LLLow": "Pmax2LLLow",
                "best_model": "best_model",
                "modelII": "ModelII_AICc",
                "nullModel": "nullModel_AICc",
                "modelIII": "ModelIII_AICc"
                }

    def __init__(self, inputfile, outputcsv="", locinames=None, pickAlleles=False, iteration=1, modelIII=True, chainlength="1e5"):

        self.inputfile = os.path.abspath(inputfile).replace("\\", "/")                     # In case someone decides to run this on windows
        self.outputcsv = outputcsv
        self.locinames = self.getAlleleList(locinames) if locinames is not None else None
        self.pickAlleles(pickAlleles) if pickAlleles is not False else False
        self.iteration = iteration
        self.rawrscript = Hzarloop.rawrscript1 if not modelIII else Hzarloop.rawrscript
        self.chainlength = chainlength

    def pickAlleles(self, outputfilename="alleleList.txt"):
        """
        parse input spreadsheet to find which of the 2 alleles are positive on
            average for negative distances.  Write to a file.
        A column header called "classification" can be used to tell this program
            which of the two alleles to use for plotting.  Putting a "1" in the
            csv column dictates whether to add the data to the sum or not.
            - If this header is not available, all of the a and b alleles on the
                negative side of the transect (aka negative distances) will be summed,
                and the greater of the 2 will be used.
        """
        classification = False 
        with open(self.inputfile, "r") as f:                                                   # open csv data file object as f
            reader = csv.reader(f)                                                             # csv reader object
            i = 0                                                                              # keeping track of index so we know when we are at the header
            for row in reader:
                i += 1
                if i == 1:                                                                     # if we are looking at the header...
                    header = row
                    complements = [[j, j+1] for j in range(len(header)) if ".a" in header[j].lower()]  # list comprehension.  gives a list of tuples of complement alleles.
                    if "classification" in header:
                        classification = True
                        clIndex = header.index("classification")
                    distIndex = header.index("distance")                                       # get index of distance column
                    colsums = [[0, 0] for j in range(len(complements))]                        # create list of tuples of 0's like [[0,0],[0,0],...,[0,0]]
                else:                                                                          # aka if not header
                    if classification:
                        cl = int(row[clIndex])
                        if cl == 1:
                            for j, indices in enumerate(complements):                          # sum the frequencies for each .a and .b allele separately
                                colsums[j][0] += float(row[indices[0]])
                                colsums[j][1] += float(row[indices[1]])
                    else:
                        dist = float(row[distIndex])                                           # get distance
                        if dist < 0:                                                           # only look at negative distances
                            for j, indices in enumerate(complements):                          # sum the frequencies for each .a and .b allele separately
                                colsums[j][0] += float(row[indices[0]])
                                colsums[j][1] += float(row[indices[1]])
        with open(outputfilename, "w") as f:
            for i, s in enumerate(colsums):
                if s[0] > s[1]:
                    # print(header[complements[i][0]])
                    f.write("{}\n".format(header[complements[i][0]]))
                else:
                    # print(header[complements[i][1]])
                    f.write("{}\n".format(header[complements[i][1]]))
        self.locinames = self.getAlleleList(outputfilename)

    def getAlleleList(self, allelenamefile):
        """
        reads in the alleles to fit in the allelenamefile file.
        """
        with open(allelenamefile, "r") as f:
            alleles = ["X{}".format(line.strip()) for line in f if len(line.strip()) != 0]  # Add X to the beginning to match R headers
            return alleles

    def serialFit(self):
        with open(self.outputcsv, "w") as f:                                            # create file object f
            writer = csv.DictWriter(f, fieldnames=self.headers, lineterminator="\n")    # writer object for writing to csv
            writer.writeheader()                                                        # writes the self.headers to file
        for i, loci in enumerate(self.locinames):                                       # "i" is the index, "loci" is the name
            cparams = {}  # converted params to match headers of output                 # param names don't match headers in csv, need to convert using conversion dictionary
            loci = "{}{}".format(loci[:-1], loci[-1].lower())
            rloci = loci.replace("_", "")
            print(loci)                                                                 # For debugging, can be commented out
            # print(rloci)
            rnum = "{}n".format(rloci[:-1])
            locivals = loci[:-2].split("_")
            suffix = loci[-2:]
            headersplit = ["ID", "lg", "read_pos", "scaffold_id", "lg_position"]
            for i, val in enumerate(locivals):
                if i == 0:
                    cparams[headersplit[i]] = "{}{}".format(val, suffix)
                else:
                    cparams[headersplit[i]] = val
            lociname = locivals[0]
            # print(lociname)
            num = "{}n".format(loci[:-1])                                               # Replace last character of name with "n", as was done in Madeline's script
            # print(num)
            with open(self.inputfile, "r") as f:
                csvstring = f.read()
                csvstring = csvstring.replace(".A", ".a").replace(".B", ".b").replace(".N", ".n")
                csvstring = csvstring.replace(loci[1:-2], rloci[1:-2])
            # rscript = self.rawrscript.format(loci, num, csvstring)                         # create R script for specific loci
            # csvstring = csvstring.replace(loci, lociname)
            with open('rfile.R', "w", encoding='utf-8') as f:                                # create file object f
                f.write(self.rawrscript.format(**{"csvtext": csvstring, "chainlength": self.chainlength}))  # write R Script
            # cmd1 = "Rscript hzar_loop.R"                                                   # Defining the command to be run in shell
            cmd1 = "Rscript --vanilla rfile.R {} {}".format(rloci, rnum)
            sp = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)                  # using Popen so that the output can be saved to python variable
            # sp = subprocess.Popen(cmd1, shell=True)
            data = str(sp.stdout.read()).split("\\n")                                        # saving output to "data" variable
            params = self.extractParams(data)                                                # extract information from data
            # print(params)
            for p in params:                                                                                              # Convert params to correct names
                if p in self.conversion:
                    cparams[self.conversion[p]] = params[p]
                else:
                    cparams[p] = params[p]
            if "Pmax" in cparams:                                                                                         # nullModel won't have Pmax or any of the following parameters
                cparams["Pmax-Pmin"] = str(round(float(cparams["Pmax"]) - float(cparams["Pmin"]), 7))                     # define and calculate Pmax - Pmin
                cparams["slope"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width"]), 7))                   # define and calculate slope
                cparams["Pmax-PminLLLow"] = str(round(float(cparams["Pmax2LLLow"]) - float(cparams["Pmin2LLHigh"]), 7))   # pmax_low - pmin_max = min delta p
                cparams["Pmax-PminLLHigh"] = str(round(float(cparams["Pmax2LLHigh"]) - float(cparams["Pmin2LLLow"]), 7))  # pmax_high - pmin_low = max delta p
                cparams["slope2LLLow"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width2LLHigh"]), 7))      # min delta p / width = min slope
                cparams["slope2LLHigh"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width2LLLow"]), 7))      # max delta p / with = max slope
            else:                                                                                                         # if Pmax isn't in cparams, we have a nullModel
                cparams["slope"] = "0"                                                                                    # nullModel implies slope = 0
            with open(self.outputcsv, "a") as f:                                                                          # file object f
                writer = csv.DictWriter(f, fieldnames=self.headers, lineterminator="\n")                                  # csv writer object "writer"
                writer.writerow(cparams)                                                                                  # write row of data to csv

    @staticmethod
    def makeFit(locinames, rawrscript, csvstring, rank):
        """
        to be used in conjunction with clusterFit.
        Given a list of locinames, and the csv data file this function will
            generate the R script to fit the data and return information
            extracted from the R output in a python dictionary to be written
            to a spreadsheet as the data is collected.
        rawrscript is the raw R script to be modified and generated.
        """
        cparamslist = []
        for i, loci in enumerate(locinames):                                            # "i" is the index, "loci" is the name
            # print("Fitting ID {}".format(loci))                                                                 # For debugging, can be commented out
            cparams = {}  # converted params to match headers of output                 # param names don't match headers in csv, need to convert using conversion dictionary
            loci = "{}{}".format(loci[:-1], loci[-1].lower())
            rloci = loci.replace("_", "")
            print(loci)                                                                 # For debugging, can be commented out
            # print(rloci)
            rnum = "{}n".format(rloci[:-1])
            locivals = loci[:-2].split("_")
            suffix = loci[-2:]
            headersplit = ["ID", "lg", "read_pos", "scaffold_id", "lg_position"]
            for i, val in enumerate(locivals):
                if i == 0:
                    cparams[headersplit[i]] = "{}{}".format(val, suffix)
                else:
                    cparams[headersplit[i]] = val
            lociname = locivals[0]
            num = "{}n".format(loci[:-1])                                               # Replace last character of name with "n", as was done in Madeline's script
            cmd1 = ["Rscript", "--vanilla", "rfile.R", rloci, rnum]
            start = time.time()
            # sp = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)             # using Popen so that the output can be saved to python variable
            # subprocess.Popen(cmd1)
            print("Fitting ID {}".format(loci))
            sp = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
            data = str(sp.stdout.read()).split("\\n")                                   # saving output to "data" variable
            # print(data)
            print("Executing R script took {} seconds.".format(time.time() - start))
            # os.remove(rname)
            print("Extracting data for loci {}.".format(loci))
            params = Hzarloop.extractParams(data)                                           # extract information from data
            # print(params)
            # print("nullModel: {}".format(params["nullModel"]))
            # print("modelII:   {}".format(params["modelII"]))
            # print("modelIII:  {}".format(params["modelIII"]))
            condition1 = (float(params["nullModel"]) - 10) < float(params["modelII"])
            condition2 = (float(params["nullModel"]) - 10) < float(params["modelIII"]) if "modelIII" in params else True
            if condition1 and condition2:
                cparams["best_model"] = "nullModel"
                cparams["slope"] = "0"
                cparams["nullModel_AICc"] = params["nullModel"]
                cparams["ModelII_AICc"] = params["modelII"]
                cparams["best_model_loglike"] = params["logLike"]
                if "modelIII" in params:
                    cparams["ModelIII_AICc"] = params["modelIII"]
            else:
                for p in params:                                                                                              # Convert params to correct names
                    if p in Hzarloop.conversion:
                        cparams[Hzarloop.conversion[p]] = params[p]
                    else:
                        cparams[p] = params[p]
                if "Pmax" in cparams:                                                                                         # nullModel won't have Pmax or any of the following parameters
                    cparams["Pmax-Pmin"] = str(round(float(cparams["Pmax"]) - float(cparams["Pmin"]), 7))                     # define and calculate Pmax - Pmin
                    cparams["slope"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width"]), 7))                   # define and calculate slope
                    cparams["Pmax-PminLLLow"] = str(round(float(cparams["Pmax2LLLow"]) - float(cparams["Pmin2LLHigh"]), 7))   # pmax_low - pmin_max = min delta p
                    cparams["Pmax-PminLLHigh"] = str(round(float(cparams["Pmax2LLHigh"]) - float(cparams["Pmin2LLLow"]), 7))  # pmax_high - pmin_low = max delta p
                    cparams["slope2LLLow"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width2LLHigh"]), 7))      # min delta p / width = min slope
                    cparams["slope2LLHigh"] = str(round(float(cparams["Pmax-Pmin"]) / float(cparams["width2LLLow"]), 7))      # max delta p / with = max slope
                else:                                                                                                         # if Pmax isn't in cparams, we have a nullModel
                    cparams["slope"] = "0"                                                                                    # nullModel implies slope = 0
            cparamslist.append(cparams)
        # print("cparamslist: {}".format(cparamslist))
        written = False
        while not written:
            try:
                with open("done.txt", "a") as f:
                    f.write("{}\n".format(loci[1:]))
                written = True
            except:
                pass
        return cparamslist

    @staticmethod
    def extractParams(data):
        """
        given the data printed from the R file, return the necessary
            parameters and values for each model.
        data - list of strings (1 for each new line)

        TODO - check that parsing modelIII and nullModel works.
        """
        # print("Exctracting data from string.")
        params = {}                                                 # The dictionary to be built on and returned
        i = 0                                                       # an index for which line we are at
        wcount = -1
        # for line in data:                                           # checking output for debugging purposes
        #     print(line)
        while i < len(data):                                        # we want to iterate over the entire file
            wcount += 1
            if wcount % 50 == 0:
                print("Running through while loop # {}.".format(wcount))
            line = data[i].strip()                                  # remove spaces at the beginning and end
            # print("line 464: {}".format(line))
            if "logLike" in line:
                i += 1
                # print(line)
                params["logLike"] = data[i].strip().split()[1]      # logLike value is on the next line.  Split on spaces and retrieve the second value.  (Python uses 0-based indexing)
            if "AICc" not in line:                                  # We don't care about any of the data before seeing "AICc"
                i += 1
            else:
                moveon = False
                while not moveon:                                   # The next 3 lines contain data for the models and their scores
                    i += 1
                    if "[1]" in data[i]:                            # if we have less than 3 models, break out of our 3 loop forloop
                        i -= 1
                        moveon = True
                    else:
                        model, score = data[i].split()
                        params[model] = score                       # Add model and score to params dictionary
                i += 1
                # print("bestmodel: {}.".format(data[i]))
                bestmodel = data[i].split()[1].replace('"', '')     # Next line declares best model
                paramnames = []
                paramvalues = []
                params["best_model"] = bestmodel                    # Add best_model to params dictionary
                if bestmodel in ["modelII", "modelIII"]:            # Handle data for modelII or modelIII.  They are formatted similarly.
                    for j in range(3):
                        if j == 2 and bestmodel != "modelIII":
                            i += 1
                        i += 1
                        paramnames += data[i].strip().split()
                        i += 1
                        paramvalues += data[i].strip().split()[1:]
                    if bestmodel == "modelIII":
                        i += 2
                        paramnames += data[i].strip().split()
                        i += 1
                        paramvalues += data[i].strip().split()[1:]
                    # print("paramnames: {}".format(paramnames))
                    # print("paramvalues: {}".format(paramvalues))
                    for j, p in enumerate(paramnames):
                        params[p] = paramvalues[j]
                else:  # best model is nullModel - no data to extract
                    i += 1
        return params


def test():
    """
    Used for debugging
    """
    print(sys.version)
    start = time.time()
    hzar = Hzarloop("ref7big.csv", "small_output/ref7model3.csv", "ref7small.txt")
    hzar.serialFit()
    end = time.time()
    print("time: {}.".format(end - start))


def testsingle():
    """
    Used for debugging
    """
    hzar = Hzarloop("hzar_10loci.csv", "small_output/singleacisstest.csv", "locinames_short.txt")
    hzar.serialFit()


def testPickAlleles():
    fname = sys.argv[1]
    output = sys.argv[2]
    hzar = Hzarloop(fname, pickAlleles=output)


if __name__ == '__main__':
    # test()
    # testsingle()
    testPickAlleles()
