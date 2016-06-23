#' Run the CONFETTI model
#'
#' @param pars Model parameters of CONFETTI. Either a vector of length 23 or a named vector of variable length. The names have to match the following list. When the parameter is unnamed the order has to correspond to the following list.
#'
#'\enumerate{
#'    \item \bold{metaSR} - Species richness of the metacommunity
#'    \item \bold{metaCV} - Coefficient of variation (sd/mean) of the metacommunity lognormal abundance distribution
#' \item \bold{m} - Immigration rate (Interval [0,1])
#' \item \bold{rmax} - Radius of interactions neighborhood (Interval [0,20])
#' \item \bold{aRec} - Parameter for recruitment competition (Interval [0,0.1])
#' \item \bold{disp.m} - Mean dispersal distance of species
#' \item \bold{disp.cv} - Coefficient of variation (sd/mean) of dispersal distance among species
#' \item \bold{CNDD.m} - Mean conspecific negative density dependence (CNDD)
#' \item \bold{CNDD.cv} - Coefficient of variation (sd/mean) of CNDD among species
#' \item \bold{pRec.m} - mean recruitment probability without competition in interval [0.001; 1]
#' \item \bold{pRec.cv} - Coefficient of variation (sd/mean) of recruitment probability among species
#' \item \bold{trade1.CNDD.pRec} - mode for trade-off between CNDD and recruitment probability. There are three different possible modes (parameterizations) for each of the three trade-offs:
#' \itemize {
#'    \item 0: no trade-off
#'    \item 1: linear: \eqn{y = a + b * x}
#'    \item 2: logarithmic: \eqn{y = a + b * log(x)}
#'    \item 3: exponential: \eqn{y = a + b * exp(c * x)}
#' }
#' \item \bold{a.CNDD.pRec} - parameters for CNDD-recruitment trade-off
#' \item \bold{b.CNDD.pRec}
#' \item \bold{c.CNDD.pRec}
#' \item \bold{trade2.CNDD.abund} - mode for trade-off between CNDD and metacommunity relative abundance:
#' \item \bold{a.CNDD.abund} - parameters for CNDD-metacommunity abundance trade-off
#' \item \bold{b.CNDD.abund}
#' \item \bold{c.CNDD.abund}
#' \item \bold{trade3.disp.pRec} - mode for trade-off between dispersal and recruitment probability:
#' \item \bold{a.disp.pRec} - parameters for dispersal-recruitment trade-off
#' \item \bold{b.disp.pRec}
#' \item \bold{c.disp.pRec}
#' }
#' @param nRep Number of replicate runs with the same parameter set
#' @param nGen Number of generations simulated. This means the model is simulates nGen*nTrees birth-death events
#' @param nSteps.out Number of time steps with model output. If nSteps.out == 1 there is only model output from the last simulated time step. If nSteps.out > 1 there is output over time. The output time interval is calculated as nGen/nSteps.out. The initial condition is included in the output. Temporal output works only on combination with avg = F.
#' @param nTrees Number of trees in the local community
#' @param Xext Size of the local community in x-direction
#' @param Yext Size of the local community in y-direction
#' @param dist.max Maximum distance (in meters) over which the spatial patterns F(r) and PCF(r) are calculated.
#' @param bin.width The distance bin-width used for the calculation of the spatial patterns F(r) and PCF(r),
#' @param meta.SAD The abundance distribution of the metacommunity: (0) uniform distribution, (1) log-normal distribution with coefficient of variation metaCV.
#' @param avg Logical variable: Should the non-temporal model output be averaged over the replicate runs? This works only if nSteps.out == 1.
#'
#' @return The model simulates virtual tree censuses with coordinates and species identities for each tree. Several summary statistics are calculated from the census data. The summary statistics are calculated for each of the nSteps.out time steps. The model returns a list with the following elements:
#' \itemize{
#'    \item \bold{census} - a dataframe with x,y coordinates and species ID for every tree in the simulated community. The census represents the last simulated time step and is only provided when nRep = 1 and avg = F.
#'    \item \bold{abundance} - a vector (or matrix) with the number of individuals for every species. Only provided when nRep = 1 and avg = F.
#'    \item \bold{species} - a dataframe with species properties: (1) the species ID, (2) the relative abundance in the metacommunity, (3) the mean dispersal distance, (4) the species-specifc conspecific negative density dependence (CNDD), (5) the species-specific recruitment probability without competition. Only provided when nRep = 1 and avg = F.
#'    \item \bold{generations} - Number of generations for which model output is calculated.
#'    \item \bold{nSpecies} - Number of species
#'    \item \bold{Shannon} - Shannon diversity index \eqn{H = - \sum log(pi) * pi}
#'    \item \bold{Simpson} - Simpson diversity index \eqn{S = 1 - \sum 1/pi^2}, where pi is the relative abundance of species i. This essentially estimates the probability of randomly drawing two individuals from different species (probability of interspecific encounter)
#'    \item \bold{SAD} - species abundance distribution with logarithmic abundance classes (1, 2-3, 4-7, 8-15, ..., >=2048 individuals)
#'    \item \bold{Area} - sampling areas for species-area relationship (SAR) in m^2
#'    \item \bold{SAR} - species-area relationship: average species numbers in quadrats of different sizes
#'    \item \bold{radius} - neighbourhood radii for calculation of F(r) and PCF(r) functions
#'    \item \bold{Fr} - Proportion of conspecific neighbors: This estimates the probability that two trees seperated by distance r belong to the same species.
#'    \item \bold{PCF} - pair-correlation function. Measure of spatial aggregation of regularity across scales. A value of 1 indicates a random pattern, PCF > 1 indicates aggregation and PCF < 1 regularity (=repulsion of individuals)
#' }
#'
#' When nRep > 1 or avg = T only the list with the summary statistics is returned, but no output on trees or single species.
#'
#' @examples
#' #Run CONFETTI with standard parameters
#' confetti.run()
#'
#' # Define your own parameter vector
#' parvec <- c(500, 2.0, 0.01, 20, 0.005, 55, 0.2, 2.0, 0.2, 1.0, 0.0,
#'             0, 1.0, 0.0, 0.0, 0, 1.0, 0.0, 0, 1.0, 0.0, 0.0)
#' out1 <- confetti.run(pars = parvec)
#'
#' # plot model output
#' plot(0:11, out1$SAD,type="b",xlab="log2(Abundance)", ylab="No. of species")
#' plot(out1$Area, out1$SAR, type="b", log="xy", xlab="Area [m2]", ylab="No. of species")
#'
#' @details There are some constraints on parameter values that override the parameter settings. Species traits, i.e. recruitment probability, dispersal distance, and CNDD are simulated from probability distributions. To avoid biologically unrealistiv values a few constraints are implemented: Recruitment probability is restricted to the interval [0.001; 1]; For dispersal distance there is a lower bound of 0.1 m, and for CNDD the lower bound is 1, which means intraspecific competition equals interspedific competition. This lower bound excludes conspecific positive density dependence.
#'
#' There can be only one trade-off that involves CNDD. If the parameter settings imply a trade-off between CNDD and recruitment as well as between CNDD and metacommunity abundance, the second one is ignored and just the trade-off with recruitment probability is simulated.
#'
#' @references May, F.; Huth, A. & Wiegand, T. (2015) Moving beyond abundance distributions: neutral theory and spatial patterns in a tropical forest Proceedings of the Royal Society of London B: Biological Sciences, 282, 20141657
#'
#' @references May, F.; Wiegand, T.; Lehmann, S. & Huth, A. (2016) Do abundance distributions and species aggregation correctly predict macroecological biodiversity patterns in tropical forests? Global Ecology and Biogeography, 25, 575-585
#'
#' @useDynLib confettiRbasic
#' @importFrom Rcpp sourceCpp

confetti.run <- function(pars = c(metaSR   = 100,
                                  metaCV   = 1.0,
                                  m        = 0.1,
                                  rmax     = 10.0,
                                  aRec     = 0.005,
                                  disp.m   = 30.0,
                                  disp.cv  = 0.0,
                                  CNDD.m   = 1.0,
                                  CNDD.cv  = 0.0,
                                  pRec.m   = 1.0,
                                  pRec.cv  = 0.0,
                                  trade1.CNDD.pRec = 0,
                                  a.CNDD.pRec = 0.0,
                                  b.CNDD.pRec = 0.0,
                                  c.CNDD.pRec = 0.0,
                                  trade2.CNDD.abund = 0,
                                  a.CNDD.abund = 0.0,
                                  b.CNDD.abund = 0.0,
                                  c.CNDD.abund = 0.0,
                                  trade3.disp.pRec = 0,
                                  a.disp.pRec = 0.0,
                                  b.disp.pRec = 0.0,
                                  c.disp.pRec = 0.0),
                         nRep = 1,
                         nGen = 100,
                         nSteps.out = 1,
                         nTrees = 10000,
                         Xext = 500,
                         Yext = 500,
                         dist.max = 100,
                         bin.width = 1,
                         meta.SAD = 1,
                         avg = FALSE
                         )
{
   my.pars <- c(metaSR   = 100,
                metaCV   = 1.0,
                m        = 0.1,
                rmax     = 10.0,
                aRec     = 0.005,
                disp.m   = 30.0,
                disp.cv  = 0.0,
                CNDD.m   = 1.0,
                CNDD.cv  = 0.0,
                pRec.m   = 1.0,
                pRec.cv  = 0.0,
                trade1.CNDD.pRec = 0,
                a.CNDD.pRec = 0.0,
                b.CNDD.pRec = 0.0,
                c.CNDD.pRec = 0.0,
                trade2.CNDD.abund = 0,
                a.CNDD.abund = 0.0,
                b.CNDD.abund = 0.0,
                c.CNDD.abund = 0.0,
                trade3.disp.pRec = 0,
                a.disp.pRec = 0.0,
                b.disp.pRec = 0.0,
                c.disp.pRec = 0.0)

   if (is.null(names(pars)) & length(pars) != 22){
      return("Error: pars needs to be a named vector or of length 22")
   }

   if (!is.null(names(pars))){
      unknown.names <- pars[!names(pars) %in% names(my.pars)]
      if (length(unknown.names) >= 1)
         return(paste("Error:",names(unknown.names),"is no valid model parameter"))
   }

   if (is.data.frame(pars)){
      pars.val <- as.numeric(pars)
      names(pars.val) <- names(pars)
      pars <- pars.val
   }

   if (length(pars) == 23)
      my.pars[1:23] <- pars
   else
      my.pars[names(pars)] <- pars

   if (my.pars["trade1.CNDD.pRec"] > 0 & my.pars["trade2.CNDD.abund"] > 0){
      cat("Warning: species-specific CNDD can only be defined by one trade-off.\nOnly the trade-off between CNDD and recruitment is considered")
      my.pars["trade2.CNDD.abund"] <- 0
   }

   if (nSteps.out == 0)
      nSteps.out <- 1

   #if(nGen <= 0)
   #   nGen <- 1

   if (nRep > 1 && nSteps.out > 1){
      cat("Warning: temporal output is ignored if nRep > 1.\nEither use nRep = 1 or nSteps.out = 1")
      nSteps.out <- 1
   }

   if (nRep == 1 && avg == TRUE){
      print("Warning: When nRep == 1, avg is set to FALSE")
      avg <- FALSE
   }

   # create matrix with replicate parameter sets
   pars.mat <- matrix(rep(my.pars,times=nRep), ncol=length(my.pars), byrow=T)

   # apply the model to each line of the model output
   out1 <- apply(pars.mat, MARGIN=1, FUN=EvalConfetti,
                 ngen = nGen,
                 nsteps_out = nSteps.out,
                 ntrees = nTrees,
                 xext = Xext,
                 yext = Yext,
                 rmax = dist.max,
                 bw1  = bin.width,
                 metaSAD = meta.SAD
                 )

   # extract the output
   nSpecies <- as.numeric(unlist(sapply(out1,"[","nSpecies")))
   Shannon <- as.numeric(unlist(sapply(out1,"[","Shannon")))
   Simpson <- as.numeric(unlist(sapply(out1,"[","Simpson")))

   SAD.list <- sapply(out1,"[","SAD")
   SAR.list <- sapply(out1,"[","SAR")
   Fr.list <- sapply(out1,"[","Fr")
   PCF.list <- sapply(out1,"[","PCF")

   Area_m2  <- out1[[1]]$Area_m2
   radius   <- out1[[1]]$radius

   if (nSteps.out == 1) {
      SAD      <- drop(matrix(unlist(SAD.list), ncol=length(SAD.list[[1]]), byrow=T))
      SAR      <- drop(matrix(unlist(SAR.list), ncol=length(SAR.list[[1]]), byrow=T))
      Fr       <- drop(matrix(unlist(Fr.list), ncol=length(PCF.list[[1]]), byrow=T))
      PCF      <- drop(matrix(unlist(PCF.list), ncol=length(PCF.list[[1]]), byrow=T))
   } else {
      SAD      <- SAD.list[[1]]
      SAR      <- SAR.list[[1]]
      Fr       <- Fr.list[[1]]
      PCF      <- PCF.list[[1]]
   }

   if (nRep == 1){
      census <- out1[[1]]$Trees
      abundance <- out1[[1]]$Abundance
      species <- out1[[1]]$Species
      generations <- out1[[1]]$Generations

      return(list(census = census,
                  abundance = abundance,
                  species = species,
                  generations = generations,
                  nSpecies = nSpecies,
                  Shannon = Shannon,
                  Simpson = Simpson,
                  SAD = SAD,
                  Area_m2 = Area_m2,
                  SAR = SAR,
                  radius = radius,
                  Fr = Fr,
                  PCF = PCF))
   } else {

      if (avg == T) {
         nSpecies <- mean(nSpecies)
         Shannon <- mean(Shannon)
         Simpson <- mean(Simpson)
         SAD <- colMeans(SAD)
         SAR <- colMeans(SAR)
         Fr <- colMeans(Fr)
         PCF <- colMeans(PCF)
      }

      return(list( nSpecies = nSpecies,
                   Shannon = Shannon,
                   Simpson = Simpson,
                   SAD = SAD,
                   Area_m2 = Area_m2,
                   SAR = SAR,
                   radius = radius,
                   Fr = Fr,
                   PCF = PCF))

   } # if avg = F
}
