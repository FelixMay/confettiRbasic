# confettiRbasic
An R interface for the spatial simulation model CONFETTI

## Installing `confettiRbasic`

```r
library(devtools)
install_github("FelixMay/confettiRbasic", build_vignettes = T)
library(confettiRbasic)
```

## Getting started with `confettiRbasic`

There are several vignettes included in the package that introduce the usage of the package.

```r
help(package = "confettiRbasic")
browseVignettes("confettiRbasic")
```

For full reference see the published articles using the CONFETTI model

## Contributing to the development of `confettiRBasic`

When you plant to adapt the code of `confettiRbasic` please use GitHub to track your changes so that all developments on the package are available. To do that please follow this protocol:

1) Fork the repo to your local GitHub account

2) Clone your forked version of the repo to your machine

`git clone git@github.com:your_user_name/confettiRbasic.git`

3) Link your local repo back to this repository

`git remote add upstream git@github.com:FelixMay/confettiRbasic.git`

4) Create a branch for your changes

`git branch new_version`

5) Checkout your branch

`git checkout new_version`

6) Make your commits on that branch and when you are done push it to your
forked copy of the repo

`git push origin new_function`

7) Submit a pull request on the GitHub website by going to your forked copy
of the repo and clicking on the pull request button 

8) After your changes are merged with master you'll want to merge that
update to master with your copies as well. 

```
git pull upstream master
git push origin master
# delete your branch as its no longer needed
git branch -d new_function
```

Before your start work on the project in the future you'll want to repeat
step 8 so that your version of the repo does not become out-of-sync
with the main repository. 

## Terms and conditions

When you publish work using the CONFETTI model, please cite the following papers:

1. May, F., Huth, A. & Wiegand, T. (2015). Moving beyond abundance distributions: neutral theory and spatial patterns in a tropical forest.
Proceedings of the Royal Society of London B: Biological Sciences, 282, 20141657.

2. May, F., Wiegand, T., Lehmann, S. & Huth, A. (2016). Do abundance distributions and species aggregation correctly predict macroecological
biodiversity patterns in tropical forests? Global Ecology and Biogeography, 25, 575–585.



