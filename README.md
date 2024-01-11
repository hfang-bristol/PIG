# [An R package `PIG`](https://github.com/hfang-bristol/PIG)

## @ Overview

> The `PIG` is an R package enabling genetics-led and network-driven drug target prioritisation for glaucoma.

## @ Installation

### 1. Install R

Please install R (version 4.3.2 or above); see https://cran.r-project.org

If installed on `Ubuntu` (assuming you have a `ROOT (sudo)` privilege), please do so below

```ruby
sudo su
# here enter your password

wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.3.2.tar.gz
tar xvfz R-4.3.2.tar.gz
cd ~/R-4.3.2
./configure
make
make check
make install
R # start R
```

### 2. Install R packages

```ruby
R # start R

# if the package 'BiocManager' not installed, please do so
if(!("BiocManager" %in% rownames(installed.packages()))) install.packages("BiocManager")

# first, install basic packages: remotes, tidyverse
BiocManager::install(c('remotes','tidyverse'), dependencies=T)

# then, install the package 'PIG' (now hosted at github)
BiocManager::install("hfang-bristol/PIG", dependencies=T, force=T)

# check the package 'PIG' successfully installed
library(help=PIG)
```


## @ Contact

> Please drop [email](mailto:fh12355@rjh.com.cn) for bug reports or other enquiries.


