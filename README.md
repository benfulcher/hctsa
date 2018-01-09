# *hctsa*, highly comparative time-series analysis

*hctsa* is a software package for running highly comparative time-series analysis using [Matlab](www.mathworks.com/products/matlab/) (full support for versions R2014b or later; for use in python cf. [pyopy](https://github.com/strawlab/pyopy)).

The software provides a code framework that allows thousands of time-series analysis features to be extracted from time series (or a time-series dataset), as well as tools for normalizing and clustering the data, producing low-dimensional representations of the data, identifying discriminating features between different classes of time series, learning multivariate classification models using large sets of time-series features, finding nearest matches to a time series of interest, and a range of other visualizations and analyses.

If you use this software, please read and cite these (open access) articles:

* B.D. Fulcher and N.S. Jones. [_hctsa_: A computational framework for automated time-series phenotyping using massive feature extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6). *Cell Systems* **5**, 527 (2017).
* B.D. Fulcher, M.A. Little, N.S. Jones [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). *J. Roy. Soc. Interface* **10**, 83 (2013).

Any feedback, as [email](mailto:ben.d.fulcher@gmail.com), [github issues](https://github.com/benfulcher/hctsa/issues) or [pull requests](https://help.github.com/articles/using-pull-requests/), is much appreciated.

### Getting started
&#x1F4D6; &#x1F4D6;
Comprehensive documentation
&#x1F4D6; &#x1F4D6;
for *hctsa* is [on gitbook](https://www.gitbook.com/book/benfulcher/highly-comparative-time-series-analysis-manual/details), which can be read online or downloaded as a pdf, epub, or mobi file.

### Downloading the repository

For users unfamiliar with git, the current version of the repository can be downloaded by simply clicking the green *Download .zip* button.

It is recommended to use the repository with git.
For this, please [make a fork](https://help.github.com/articles/fork-a-repo/) of it, clone it to your local machine, and then set an [upstream remote](https://help.github.com/articles/fork-a-repo/#step-3-configure-git-to-sync-your-fork-with-the-original-spoon-knife-repository) to keep it synchronized with the main repository e.g., using the following code:
```
git remote add upstream git://github.com/benfulcher/hctsa.git
```
(make sure that you have [generated an ssh key](https://help.github.com/articles/generating-ssh-keys/) and associated it with your github account).

You can then update to the latest stable version of the repository by pulling the master branch to your local repository:
```
git pull upstream master
```

For analyzing specific datasets, we recommend working outside of the repository so that incremental updates can be pulled from the upstream repository.
Details on how to merge the latest version of the repository with the local changes in your fork can be found [here](https://help.github.com/articles/syncing-a-fork/).

## *hctsa* licenses

### Internal licenses

There are two licenses applied to the core parts of the repository:

1. Sections of the repository required to compute features from time-series data is licensed as [GNU General Public License version 3](http://www.gnu.org/licenses/gpl-3.0.en.html).

2. Sections implementing the framework for running *hctsa* analyses and visualizations is licensed as the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).

To use *hctsa* for commercial applications, please contact [Ben Fulcher](ben.d.fulcher@gmail.com).

A range of external code packages are provided in the **Toolboxes** directory of the repository, and each have their own associated license (as outlined below).

### External packages and dependencies

The following [Matlab toolboxes](https://mathworks.com/programs/nrd/matlab-toolbox-price-request.html?ref=ggl&s_eid=ppc_18665571802&q=matlab%20toolboxes%20price) are used by *hctsa* and are required for full functionality of the software.
In the case that some toolboxes are unavailable, the *hctsa* software can still be used, but only a reduced set of time-series features will be computed.

1. Statistics Toolbox
2. Signal Processing Toolbox
3. Curve Fitting Toolbox
4. System Identification Toolbox
5. Wavelet Toolbox
6. Econometrics Toolbox

---

The following time-series analysis packages are provided with the software (in the **Toolboxes** directory), and are used by our main feature extraction algorithms to compute meaningful structural features from time series:

* [*TISEAN* package for nonlinear time-series analysis, version 3.0.1](http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html) (GPL license).
* [*TSTOOL* package for nonlinear time-series analysis, version 1.2](http://www.dpi.physik.uni-goettingen.de/tstool/) (GPL license).
* Joseph T. Lizier's [Java Information Dynamics Toolkit (JIDT)](https://github.com/jlizier/jidt) for studying information-theoretic measures of computation in complex systems, version 1.3 (GPL license).
* Time-series analysis code developed by [Michael Small](http://staffhome.ecm.uwa.edu.au/~00027830/code.html) (unlicensed).
* Max Little's [Time-series analysis code](http://www.maxlittle.net/software/index.php) (GPL license).
* Sample Entropy code from [Physionet](http://www.physionet.org/faq.shtml#license) (GPL license).
* [*ARFIT* Toolbox for AR model estimation](http://climate-dynamics.org/software/#arfit) (unlicensed).
* [*gpml* Toolbox for Gaussian Process regression model estimation, version 3.5](http://www.gaussianprocess.org/gpml/code/matlab/doc/) (FreeBSD license).
* Danilo P. Mandic's [delay vector variance code](http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm) (GPL license).
* [Cross Recurrence Plot Toolbox](http://tocsy.pik-potsdam.de/CRPtoolbox/) (GPL license)
* Zoubin Ghahramani's [Hidden Markov Model (HMM) code](http://mlg.eng.cam.ac.uk/zoubin/software.html) (MIT license).
* Danny Kaplan's Code for embedding statistics (GPL license).
* Two-dimensional histogram code from Matlab Central (BSD license).
* Various histogram and entropy code by Rudy Moddemeijer (unlicensed).

## Publications

See the following publications for examples of *hctsa* use:
* ***Implementation paper introducing the hctsa package, with applications to high throughput phenotyping of C. Elegans and Drosophila movement time series*** &#x1F4D7; : B. D. Fulcher & N. S. Jones. _hctsa_: A Computational Framework for Automated Time-Series Phenotyping Using Massive Feature Extraction. *Cell Systems* **5**, 527 (2017). [Link](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6).
* ***Introduction to feature-based time-series analysis*** &#x1F4D7; : B. D. Fulcher. Feature-based time-series analysis. *arXiv* 1709.08055 (2017) [Link](https://arxiv.org/abs/1709.08055).
* ***Application to fMRI data*** &#x1F4D7; : S. S. Sethi, V. Zerbi, N. Wenderoth, A. Fornito, B. D. Fulcher. Structural connectome topology relates to regional BOLD signal dynamics in the mouse brain. *Chaos* **27**, 047405 (2017). [Link](http://aip.scitation.org/doi/10.1063/1.4979281), [preprint](http://biorxiv.org/lookup/doi/10.1101/085514).
* ***Application to time-series data mining*** &#x1F4D7; : B. D. Fulcher & N. S. Jones. Highly comparative feature-based time-series classification. *IEEE Trans. Knowl. Data Eng.* **26**, 3026 (2014). [Link](http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=6786425).
* ***Application to fetal heart rate time series*** &#x1F4D7; : B. D. Fulcher, A. E. Georgieva, C. W. G. Redman, N. S. Jones. Highly comparative fetal heart rate analysis. *34th Ann. Int. Conf. IEEE EMBC* 3135 (2012). [Link](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6346629).
* ***Original paper, showing that the behavior of thousands of time-series methods on thousands of different time series can provide structure to the interdisciplinary time-series analysis literature*** &#x1F4D7; : B. D. Fulcher, M. A. Little, N. S. Jones. Highly comparative time-series analysis: the empirical structure of time series and their methods. *J. Roy. Soc. Interface* **10**, 20130048 (2013). [Link](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full).

## Acknowledgements

Many thanks go to [Romesh Abeysuriya](https://github.com/RomeshA) for helping with the mySQL database set-up and install scripts, and [Santi Villalba](https://github.com/sdvillal) for lots of helpful feedback and advice on the software.

## Related resources

### _hctsa_ datasets
There are a range of open datasets with pre-computed _hctsa_ features.
(If you have data to share and host, let me know and I'll add it to this list):
* [1000 empirical time series](https://figshare.com/articles/1000_Empirical_Time_series/5436136)
* [C. elegans movement speed data](https://figshare.com/articles/Highly_comparative_time-series_analysis_of_Caenorhabditis_elegans_movement_speed/3863559) and associated [analysis code](https://github.com/benfulcher/hctsa_phenotypingWorm).
* [Drosophila movement speed](https://figshare.com/articles/Highly_comparative_time-series_analysis_of_Drosophila_melanogaster_movement_speed/3863553) and associated [analysis code](https://github.com/benfulcher/hctsa_phenotypingFly).

### Comp-Engine Time Series

An accompanying web resource for this project is [Comp-Engine Time Series](http://www.comp-engine.org/timeseries), which allows users to compare thousands of diverse types of time-series analysis code and time-series data.
This website also allows you to download large volumes of univariate time-series data.
Note that the code files on Comp-Engine Time Series are based on an early implementation and rarely match with the updated features and functions contained in this repository.

### `pyopy`

This excellent repository allows users to run *hctsa* software from within python: [pyopy](https://github.com/strawlab/pyopy).

### `tsfresh`

Native python time-series code to extract hundreds of time-series features, with in-built feature filtering, is [tsfresh](https://github.com/blue-yonder/tsfresh).

### `tscompdata` and `tsfeatures`

These R packages are by [Rob Hyndman](https://twitter.com/robjhyndman).
The first, [`tscompdata`](https://github.com/robjhyndman/tscompdata), makes available existing collections of time-series data for analysis.
The second, [`tsfeatures`](https://github.com/robjhyndman/tsfeatures), includes implementations of a range of time-series features.

### `pyunicorn`

A python-based nonlinear time-series analysis and complex systems code package, [pyunicorn](http://scitation.aip.org/content/aip/journal/chaos/25/11/10.1063/1.4934554).
