##hctsa
###Highly comparative time-series analysis code repository

This repository draws on the work published as:

B. D. Fulcher, M. A. Little, N. S. Jones (2013) [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). J. Roy. Soc. Interface. 10, 83.

As part of this work, hundreds of pieces of time-series analysis code were developed that produce thousands of features summarizing structural properties of time series.
An accompanying web resource for this project is [comp-engine time series](www.comp-engine.org/timeseries), which allows users to compare (and more recently upload) thousands of diverse types of time-series analysis code and time-series data.
This code repository is a more in-depth accompaniment to this work, that helps the user set up a mySQL database through Matlab and perform highly comparative time-series analysis on custom datasets using the latest set of code and input settings.

To use the repository, please make a private fork of it, clone it to your local machine, and then set an [upstream remote](https://help.github.com/articles/fork-a-repo/#step-3-configure-git-to-sync-your-fork-with-the-original-spoon-knife-repository) to keep it synchronized with the main repository e.g., using the following code:
```
git remote add upstream git://github.com/SystemsAndSignalsGroup/hctsa.git
```
(make sure that you have [generated an ssh key](https://help.github.com/articles/generating-ssh-keys/) and associated it with your github account).

You can then update to the latest stable version of the repository by pulling the master branch to your local repository:
```
git pull upstream master
```

We recommend keeping the repository without any changes so that the latest version can be pulled from the upstream repository.
However, details on how to merge the latest version of the repository with the local changes in your fork can be found [here](https://help.github.com/articles/syncing-a-fork/).

Any feedback is hugely helpful ([email me](mailto:ben.d.fulcher@gmail.com)) and, in particular, any improvements to the code would be _much_ appreciated in the form of [issues](https://github.com/SystemsAndSignalsGroup/hctsa/issues) or [pull requests](https://help.github.com/articles/using-pull-requests/).

Many thanks go to Romesh Abeysuriya for helping with the database set-up process and install scripts.
