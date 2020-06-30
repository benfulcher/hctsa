# :wavy_dash: *hctsa* :wavy_dash:: highly comparative time-series analysis

*hctsa* is a software package for running highly comparative time-series analysis using [Matlab](https://www.mathworks.com/products/matlab/) (full support for versions R2018b or later).

The software provides a code framework that enables the extraction of thousands of time-series features from a time series (or a time-series dataset).
It also provides a range of tools for visualizing and analyzing the resulting time-series feature matrix, including:
1. Normalizing and clustering the data,
2. Producing low-dimensional representations of the data,
3. Identifying and interpreting discriminating features between different classes of time series,
4. Learning multivariate classification models.

**Feel free to [email me](mailto:ben.d.fulcher@gmail.com) for help with real-world applications of _hctsa_** :nerd_face:

### Acknowledgement :+1:
If you use this software, please read and cite these open-access articles:

* B.D. Fulcher and N.S. Jones. [_hctsa_: A computational framework for automated time-series phenotyping using massive feature extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6). *Cell Systems* **5**, 527 (2017).
* B.D. Fulcher, M.A. Little, N.S. Jones. [Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full). *J. Roy. Soc. Interface* **10**, 83 (2013).

Feedback, as [email](mailto:ben.d.fulcher@gmail.com), [github issues](https://github.com/benfulcher/hctsa/issues) or [pull requests](https://help.github.com/articles/using-pull-requests/), is much appreciated.

**For commercial use of *hctsa*, including licensing and consulting, contact [Engine Analytics](http://www.engineanalytics.org/).**

## Getting Started :blush:

### Documentation &#x1F4D6;
__Comprehensive documentation__ for *hctsa*, from getting started through to more advanced analyses is on [gitbook](https://hctsa-users.gitbook.io/hctsa-manual).

### Downloading the repository :arrow_down:

For users unfamiliar with git, the current version of the repository can be downloaded by simply clicking the green _Clone or download_ button, and then clicking _Download .zip_.

It is recommended to use the repository with git.
For this, please [make a fork](https://help.github.com/articles/fork-a-repo/) of it, clone it to your local machine, and then set an [upstream remote](https://help.github.com/articles/fork-a-repo/#step-3-configure-git-to-sync-your-fork-with-the-original-spoon-knife-repository) to keep it synchronized with the main repository e.g., using the following code:
```
git remote add upstream git://github.com/benfulcher/hctsa.git
```
(make sure that you have [generated an ssh key](https://help.github.com/articles/generating-ssh-keys/) and associated it with your Github account).

You can then update to the latest stable version of the repository by pulling the master branch to your local repository:
```
git pull upstream master
```

For analyzing specific datasets, we recommend working outside of the repository so that incremental updates can be pulled from the upstream repository.
Details on how to merge the latest version of the repository with the local changes in your fork can be found [here](https://help.github.com/articles/syncing-a-fork/).

## Related resources

### _CompEngine_ :collision:

[_CompEngine_](http://www.comp-engine.org) is an accompanying web resource for this project.
It is a self-organizing database of time-series data that allows users to upload, explore, and compare thousands of diverse types of time-series data.
This vast and growing collection of time-series data can also be downloaded.
You can read more about it in our [&#x1F4D9;preprint](https://arxiv.org/abs/1905.01042).

### _catch22_ :two::two:

Is over 7000 just a few too many features for your application?
Do you not have access to a Matlab license?
_catch22_ has all of your faux-rhetorical questions covered.
This reduced set of 22 features, determined through a combination of classification performance and mutual redundancy as explained in [this paper](https://arxiv.org/abs/1901.10200v2), is available [here](https://github.com/chlubba/catch22) as an efficiently coded C implementation with wrappers for python and R.

### _hctsa_ datasets and example workflows :floppy_disk:
There are a range of open datasets with pre-computed _hctsa_ features, as well as some examples of _hctsa_ workflows.
* [_C. elegans_ movement speed data](https://figshare.com/articles/Highly_comparative_time-series_analysis_of_Caenorhabditis_elegans_movement_speed/3863559) and associated [analysis code](https://github.com/benfulcher/hctsa_phenotypingWorm).
* [Drosophila movement speed](https://figshare.com/articles/Highly_comparative_time-series_analysis_of_Drosophila_melanogaster_movement_speed/3863553) and associated [analysis code](https://github.com/benfulcher/hctsa_phenotypingFly).
* [1000 empirical time series](https://figshare.com/articles/1000_Empirical_Time_series/5436136)

(If you have data to share and host, let me know and I'll add it to this list)

### Running _hctsa_ on a cluster :computer:

Matlab code for computing features for an initialized `HCTSA.mat` file, by distributing the computation across a large number of cluster jobs (using pbs or slurm schedulers) is [here](https://github.com/benfulcher/distributed_hctsa).

## Publications :closed_book:

Here we provide a list of publications that have used _hctsa_.

### Our publications

Where journal articles (&#x1F4D7;) are not open access, we also provide a link to the preprint (&#x1F4D9;).
Links to Github code repositories (:octocat:) are provided where appropriate.

#### The development of *hctsa* and other resources for feature-based time-series analysis

See the following publications for details of how the highly-comparative approach to time-series analysis has developed since our initial publication in 2013:

* **We reduced *hctsa* down to a reduced set of 22 efficiently coded features**.
  * [&#x1F4D7; CAnonical Time-series CHaracteristics](https://doi.org/10.1007/s10618-019-00647-x). *Data Mining and Knowledge Discovery* **33**, 1821 (2019).
  * [:octocat: Code](https://github.com/chlubba/catch22).
* **We developed a software package for highly-comparative time-series analysis, *hctsa*** (includes applications to high throughput phenotyping of _C. Elegans_ and Drosophila movement time series).
  * [&#x1F4D7; _hctsa_: A Computational Framework for Automated Time-Series Phenotyping Using Massive Feature Extraction](http://www.cell.com/cell-systems/fulltext/S2405-4712\(17\)30438-6).   *Cell Systems* **5**, 527 (2017).
  * [:octocat: Code (fly)](https://github.com/benfulcher/hctsa_phenotypingFly) &
  [:octocat: Code (worm)](https://github.com/benfulcher/hctsa_phenotypingWorm).
* **Introduction to using features for time-series analysis**
  * [&#x1F4D7; Feature-based time-series analysis](https://www.crcpress.com/Feature-Engineering-for-Machine-Learning-and-Data-Analytics/Dong-Liu/p/book/9781138744387).
    *Feature Engineering for Machine Learning and Data Analytics*, CRC Press (2018).
  * [&#x1F4D9; Preprint](https://arxiv.org/abs/1709.08055).
* **The behavior of thousands of time-series methods on thousands of different time series can be used to organize an interdisciplinary time-series analysis literature**
  * [&#x1F4D7; Highly comparative time-series analysis: the empirical structure of time series and their methods](http://rsif.royalsocietypublishing.org/content/10/83/20130048.full).
    *J. Roy. Soc. Interface* (2013).

#### Applications of *hctsa*
We have used *hctsa* to:
* **Distinguish targeted perturbations to mouse fMRI dynamics**
  * [&#x1F4D7; Cortical Excitation:Inhibition Imbalance Causes Abnormal Brain Network Dynamics as Observed in Neurodevelopmental Disorders](https://academic.oup.com/cercor/advance-article/doi/10.1093/cercor/bhaa084/5823074). *Cerebral Cortex* (2020).
  * [:octocat: Code](https://github.com/benfulcher/hctsa_DREADD).
* **Connect structural brain connectivity to fMRI dynamics**
  * [&#x1F4D7; Structural connectome topology relates to regional BOLD signal dynamics in the mouse brain](http://aip.scitation.org/doi/10.1063/1.4979281). *Chaos* (2017).
  * [&#x1F4D9; Preprint](http://biorxiv.org/lookup/doi/10.1101/085514).
* **Distinguish time-series patterns for data-mining applications**
  * [&#x1F4D7; Highly comparative feature-based time-series classification](http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=6786425). *IEEE Trans. Knowl. Data Eng.* (2014).
* **Classify babies with low blood pH from fetal heart rate time series**
  * [&#x1F4D7; Highly comparative fetal heart rate analysis](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6346629). *34th Ann. Int. Conf. IEEE EMBC* (2012).

### Other Publications

_hctsa_ has been used by others in neuroscience, engineering, and biomedicine.
An updated list of publications using _hctsa_ is on this [wiki page](https://github.com/benfulcher/hctsa/wiki/Publications-using-hctsa).

## *hctsa* licenses

### Internal licenses

There are two licenses applied to the core parts of the repository:

1. The framework for running *hctsa* analyses and visualizations is licensed as the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).
A license for commercial use is available from [Engine Analytics](http://www.engineanalytics.org/).

2. Code for computing features from time-series data is licensed as [GNU General Public License version 3](http://www.gnu.org/licenses/gpl-3.0.en.html).

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
* Sample Entropy code from [Physionet](https://archive.physionet.org/faq.shtml#license) (GPL license).
* [*ARFIT* Toolbox for AR model estimation](http://climate-dynamics.org/software/#arfit) (unlicensed).
* [*gpml* Toolbox for Gaussian Process regression model estimation, version 3.5](http://www.gaussianprocess.org/gpml/code/matlab/doc/) (FreeBSD license).
* Danilo P. Mandic's [delay vector variance code](http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm) (GPL license).
* [Cross Recurrence Plot Toolbox](http://tocsy.pik-potsdam.de/CRPtoolbox/) (GPL license)
* Zoubin Ghahramani's [Hidden Markov Model (HMM) code](http://mlg.eng.cam.ac.uk/zoubin/software.html) (MIT license).
* Danny Kaplan's Code for embedding statistics (GPL license).
* Two-dimensional histogram code from Matlab Central (BSD license).
* Various histogram and entropy code by Rudy Moddemeijer (unlicensed).

## Other time-series analysis resources

Other good resources for time-series analysis, e.g., in other programming languages (python and R) are listed here.

### `pyopy`

This excellent repository allows users to run *hctsa* software from within python: [`pyopy`](https://github.com/strawlab/pyopy).

### `hctsaAnalysisPython`

Some beginner-level python code for analyzing the results of _hctsa_ calculations is [here](https://github.com/benfulcher/hctsaAnalysisPython).

### Generating time-series data from synthetic models

A Matlab repository for generating time-series data from diverse model systems is [here](https://github.com/benfulcher/TimeSeriesGeneration).

### `tsfresh`

Native python time-series code to extract hundreds of time-series features, with in-built feature filtering, is [tsfresh](https://github.com/blue-yonder/tsfresh); cf. [their paper](https://www.sciencedirect.com/science/article/pii/S0925231218304843).

### `tscompdata` and `tsfeatures`

These R packages are by [Rob Hyndman](https://twitter.com/robjhyndman).
The first, [`tscompdata`](https://github.com/robjhyndman/tscompdata), makes available existing collections of time-series data for analysis.
The second, [`tsfeatures`](https://github.com/robjhyndman/tsfeatures), includes implementations of a range of time-series features.

### `TSFEL`

[TSFEL](https://github.com/fraunhoferportugal/tsfel), 'Time Series Feature Extraction Library', is a python package with implementations of 60 simple time-series features (with unit tests).

### `Khiva`

[Khiva](https://github.com/shapelets/khiva) is an open-source library of efficient algorithms to analyse time series in GPU and CPU.

### `pyunicorn`

A python-based nonlinear time-series analysis and complex systems code package, [pyunicorn](http://scitation.aip.org/content/aip/journal/chaos/25/11/10.1063/1.4934554).

### `TSFuse` (python)
[TSFuse](https://github.com/arnedb/tsfuse) can extract features from multivariate time series.

## Acknowledgements :wave:

Many thanks go to [Romesh Abeysuriya](https://github.com/RomeshA) for helping with the mySQL database set-up and install scripts, and [Santi Villalba](https://github.com/sdvillal) for lots of helpful feedback and advice on the software.
