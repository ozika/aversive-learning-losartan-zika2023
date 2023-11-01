<img src ="img/medicine.png" width="100">


# Description

This repo is associated with the following paper: Zika, O., Appel, J., Klinge, C., Shkreli, L., Browning, M., Wiech, K. and Reinecke, A. (2023), *Reduction of aversive learning rates in Pavlovian conditioning by angiotensin II antagonist losartan.*

# Citation
Please cite this code and data. 


```
Zika, O., Appel, J., Klinge, C., Shkreli, L., Browning, M., Wiech, K. and Reinecke, A. (2023). Reduction of aversive learning rates in Pavlovian conditioning by angiotensin II antagonist losartan. 
```


It contains **data** and **scripts** to reproduce the results in the main text.

# Instructions

The main analysis scripts are organised in a single `r-markdown` notebook. The computational environment can be reproduced using the `renv` package for R environments.

**Requirements**
- R 3.6.0+
- RStudio (or other r-markdown editor)
- set up SSH key with github.com (see [instructions here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account))

**Scripts tested on**

Operating systems: `Ubuntu 18.04.6 LTS`,  `macOS Catalina 10.15.7`, `macOS Monterey 12.3.1`
R versions: `3.6.0 (2019-04-26)`, `4.0.0 (2020-04-24)`, `4.1.2 (2021-11-01)`

**Expected runtime**
- about 10 minutes on a normal computer

**Instructions**
1. Glone git repo (this will create a folder with all scripts and data called `aversive-learning-losartan-zika2023`)

```bash
git clone  git@github.com:ozika/aversive-learning-losartan-zika2023.git
```

2. In RStudio open `scripts/main_analyses.Rmd` 
3. Run. The script should load the `renv.lock` file and use it to download and install all necessary packages.
4. To check that all worked, the script should create a `output/figures/` folder and print the paper figures into it.

# Licence
Code: GNU GPLv3 (see LICENSE.md)  
Data: CC BY-NC 4.0 International license [see here](https://creativecommons.org/licenses/by-nc/4.0/)