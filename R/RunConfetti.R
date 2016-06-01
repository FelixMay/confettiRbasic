#' Run the CONFETTI model
#'
#' @param pars A vector of eight numbers that include the model parameters for CONFETTI. The vector can be named but need not. Note that the order of parameters matters.
#'
#'\enumerate{
#'    \item \bold{metaSR} - Species richness of the metacommunity
#'    \item \bold{metaCV} - Coefficient of variation (sd/mean) of the metacommunity lognormal abundance distribution
#' \item \bold{m} - Immigration rate (Interval [0,1])
#' \item \bold{rmax} - Radius of interactions neighborhood (Interval [0,20])
#' \item \bold{aRec} - Parameter for recruitment competition (Interval [0,0.1])
#' \item \bold{disp} - Mean dispersal distance
#' \item \bold{CNDD.m} - Mean conspecific negative density dependence (CNDD)
#' \item \bold{CNDD.cv} - Coefficient of variation (sd/mean) of CNDD among species
#' }
#' @param nRep Number of replicate runs with the same parameter set
#' @param nGen Number of generations simulated. This means the model is simulates nGen*nTrees birth-death events
#' @param nTrees Number of trees in the local community
#' @param Xext Size of the local community in x-direction
#' @param Yext Size of the local community in y-direction
#' @param avg logical variable: Should the model output be averaged over the replicate runs
#'
#' @return The model simulates virtual tree censuses with coordinates and species identities for each tree. The census data is only provided as output when nRep = 1 and avg = F. In this case the output is s list with three elements:
#' \enumerate{
#'    \item \bold{census} - a dataframe with x,y coordinates and species ID for every tree in the simulated community
#'    \item \bold{abundance} - a names vector with the number of individuals for every species
#'    \item \bold{summary.stats} a list with eigth elements that includes non-spatial and spatial summary statistics of the simulated community. The list elements are:
#' \enumerate{
#'    \item \bold{nSpecies} - Number of species
#'    \item \bold{Shannon} - Shannon diversity index \eqn{H = - \sum log(pi) * pi}
#'    \item \bold{Simpson} - Simpson diversity index \eqn{S = 1 - \sum 1/pi^2}, where pi is the relative abundance of species i. This essentially estimates the probability of randomly drawing two individuals from different species (probability of interspecific encounter)
#'    \item \bold{SAD} - species abundance distribution with logarithmic abundance classes (1, 2-3, 4-7, 8-15, ..., >=2048 individuals)
#'    \item \bold{Area} - sampling areas for species area relationship in m^2
#'    \item \bold{SAR} - species-area relationship: average species numbers in quadrats of different sizes
#'    \item \bold{Fr} - Proportion of conspecific neighbors: This estimates the probability that two trees seperated by distance r belong to the same species.
#'    \item \bold{PCF} - pair-correlation function. Measure of spatial aggregation of regularity across scales. A value of 1 indicates a random pattern, PCF > 1 indicates aggregation and PCF < 1 regularity (=repulsion of individuals)
#' }
#' }
#' When nRep > 1 or avg = T only the list with the summary statistics is returned.
#'
#' @examples
#' #Run CONFETTI with standard parameters
#' confetti.run()
#'
#' # Define your own parameter vector
#' parvec <- c(500, 2.0, 0.01, 20, 0.005, 55, 2.0, 0.2)
#' out1 <- confetti.run(pars = parvec)
#'
#' # plot model output
#' plot(0:11, out1$SAD,type="b",xlab="log2(Abundance)", ylab="No. of species")
#' plot(out1$Area, out1$SAR, type="b", log="xy", xlab="Area [m2]", ylab="No. of species")
#'
#' @references May, F.; Huth, A. & Wiegand, T. (2015) Moving beyond abundance distributions: neutral theory and spatial patterns in a tropical forest Proceedings of the Royal Society of London B: Biological Sciences, 282, 20141657
#'
#' @references May, F.; Wiegand, T.; Lehmann, S. & Huth, A. (2016) Do abundance distributions and species aggregation correctly predict macroecological biodiversity patterns in tropical forests? Global Ecology and Biogeography, 25, 575-585
#'
#' @useDynLib confettiRbasic
#' @importFrom Rcpp sourceCpp


confetti.run <- function(pars = c(metaSR = 300,
                                  metaCV = 1.0,
                                  m      = 0.1,
                                  rmax   = 10,
                                  aRec   = 0.005,
                                  disp   = 30,
                                  CNDD.m = 1.0,
                                  CNDD.cv = 0.0),
                         nRep = 1,
                         nGen = 100,
                         nTrees = 10000,
                         Xext = 500,
                         Yext = 500,
                         avg = FALSE)
{
   # create matrix with replicate parameter sets
   pars <- as.numeric(pars)
   pars.mat <- matrix(rep(pars,times=nRep), ncol=length(pars), byrow=T)

   # apply the model to each line of the model output
   out1 <- apply(pars.mat, MARGIN=1, FUN=EvalConfetti,
                 ngen = nGen,
                 ntrees = nTrees,
                 xext = Xext,
                 yext = Yext)

   # extract the output
   nSpecies <- as.numeric(unlist(sapply(out1,"[","nSpecies")))
   Shannon <- as.numeric(unlist(sapply(out1,"[","Shannon")))
   Simpson <- as.numeric(unlist(sapply(out1,"[","Simpson")))

   SAD.list <- sapply(out1,"[","SAD")
   SAR.list <- sapply(out1,"[","SAR")
   Fr.list <- sapply(out1,"[","Fr")
   PCF.list <- sapply(out1,"[","PCF")

   Area_m2  <- out1[[1]]$Area_m2

   SAD      <- matrix(unlist(SAD.list),ncol=length(SAD.list[[1]]),byrow=T)
   SAR      <- matrix(unlist(SAR.list),ncol=length(SAR.list[[1]]),byrow=T)
   Fr       <- matrix(unlist(Fr.list),ncol=length(PCF.list[[1]]),byrow=T)
   PCF      <- matrix(unlist(PCF.list),ncol=length(PCF.list[[1]]),byrow=T)

   summary.stats <- list()

   # if necessary calculate the average patterns
   if (avg == TRUE){
      summary.stats$nSpecies <- mean(nSpecies)
      summary.stats$Shannon <- mean(Shannon)
      summary.stats$Simpson <- mean(Simpson)
      summary.stats$SAD <- colMeans(SAD)
      summary.stats$Area <- Area_m2
      summary.stats$SAR <- colMeans(SAR)
      summary.stats$Fr <- colMeans(Fr)
      summary.stats$PCF <- colMeans(PCF)

      return( summary.stats = summary.stats)

   } else {
      if (nRep == 1){
         census <- data.frame(X = out1[[1]]$X,
                              Y = out1[[1]]$Y,
                              SpecID = factor(out1[[1]]$SpecID))
         abundance <- table(out1[[1]]$SpecID)
      }
      summary.stats$nSpecies <- nSpecies
      summary.stats$Shannon <- Shannon
      summary.stats$Shannon <- Simpson
      summary.stats$SAD <- SAD
      summary.stats$Area <- Area_m2
      summary.stats$SAR <- SAR
      summary.stats$Fr <- Fr
      summary.stats$PCF <- PCF
      return(list(census = census, abundance = abundance, summary.stats = summary.stats))
   }
}
