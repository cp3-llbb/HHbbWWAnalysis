# HH Machine Learning 

Software to do multi classification for HH->bbWW analysis

These scripts are used to make hyperparameter scans with Talos and learning on Keras

All scripts are **python3**

## Getting Started

This software is intended to work on Ingrid/Manneback 

### Prerequisites

Modules you will need to load
```
module load root/6.12.04-sl7_gcc73 boost/1.66.0_sl7_gcc73 gcc/gcc-7.3.0-sl7_amd64 python/python36_sl7_gcc73 slurm/slurm_utils 

```

### Installing required python packages 

Below are the required packages that can be installed with pip.

If you do not have sysadmin rights, do not forget to use ``` pipi install --user ...  ```

- [Tensorflow](https://www.tensorflow.org/install/pip) (neural networks learning)
    -> Already on Ingrid, do not install local verison 
- [Keras](https://pypi.org/project/Keras/) (wraper around Tensorflow)
    -> Already on Ingrid, do not install local verison 
- [Talos](https://pypi.org/project/talos/) (hyperparameter scans)
    -> Special branch already available in this directory (not need to install)
- [Root_numpy](https://pypi.org/project/root-numpy/) (From ROOT trees to numpy arrays)
    -> Already on Ingrid, do not install local verison 
- [Seaborn](https://pypi.org/project/seaborn/) (Data Visualization)
- [Numpy](https://pypi.org/project/numpy/) (Data manipulation)
- [Pandas](https://pypi.org/project/pandas/) (Useful to manipulate numpy arrays altogether)
- [Astetik](https://pypi.org/project/astetik/) (Simplified templates of seaborn required by Talos)
- [Enlighten](https://pypi.org/project/enlighten/) (Practical process bar)
- [Scipy](https://pypi.org/project/scipy/) (Data processing)
- [PrettyTable](https://pypi.org/project/PrettyTable/) (table printing)

## Usual workflow

Most of the tweaks are done in one file : parameters.py

They will be described in details in the next subsections

Then we will detail the usual workflow of the hyperparameter scans

### Configuration script

parameters.py contains the global information required by all the scripts, all the variables are accessed via `parameters.something`

They will be decribed in the following 
- paths : where the script will be running, produce the output and models
- ratios : 
    - classic : one part goes for training, one for evaluation of the model (used in hyperparameter scans) and one for producing an output for testing
        A check is done to make sure the total accounts to 1
    - cross-validation : see section below
- Slurm parameters : will be provided to submit_on_slurm.py
- Names : includes 
    - physics config : config YAML (see below), lumi and eras dict
    - categories : which categories to take from config yaml, nodes for classes and channels
    - suffix : used to generate the mask and scaler (see explanation below)
    - cache : to cache the data (see later)
    - json files : contain xsec and event weight sum information
    - resume : name of model to be retrained (very rare)
    - output batch size : for producing the output (just goes faster)
    - split_name : split the output root file per tag or sample name (see later)
- evaluation criterion : to select the best model (based on val_loss or evaluation error, latter better)
- Callbacks for learning :
    - early_stopping_params : see https://keras.io/api/callbacks/early_stopping/
    - reduceLR_params : see https://keras.io/api/callbacks/reduce_lr_on_plateau/
- Scan dictionary : 
    Keys are the names of the parameters we want to scan
    Values are the possible combinations (must always be a list, even for single item)
    Repetition : number of times one hyperpameter set needs to be used (almost all the time : 1)
- Variables (can use any ROOT tricks):
    - cut : cut for data importation
    - weights : what branch to use for sampling weights in the learning
    - inputs : list of branches to be used as training variables
    - outputs : list of branches to be used as training targets
        Note : for branches that are not in tree but will be added later (eg tag) : use $string ($ will be removed after)
    - other_variables : other variables you want to keep in the tree but not use as inputs not targets
- make_dtype : this is because we use root_numpy to produce the root files and it does not like '.', '(', ')', '-' or '*'

### Sample dict ###

Yaml dict containing two keys :
- sampleDir : directory to be added in front of all the samples
- sampleDict : dict of all the samples per era
    -> key : era
        -> key : category + channel + node
        -> val : list of samples

*Tip* : A helper `produceSampleList.py` based on analysis bamboo yaml file is available.

### Workflow 

After you have chosen all the parameters, the first try can be on local

Usual command : 
``` 
python HHMachineLearning.py (args) --scan name_of_scan 
```
The args depend on what you have hardcoded in HHMachineLearning.py

*Note* : all the hyperparameter combinations will be run sequentially, this might take time ... 

*Tip* : use one combination only (only lists with one item) and small number of epochs to check everything works

The products a the scripts are :
    - csv file : contains the parameters in the scan, loss, acc and error
    - zip file : contains model architecture+weights, results in the csv, plus other details
    
*Tip* : You can either unzip the .zip and load the json and h5 files with the classic method ([here](https://machinelearningmastery.com/save-load-keras-deep-learning-models/)).
Or you can use the `Restore` method of Talos on the zip archive directly (but you need to submit the preprocessing layer specifically, see code in NeuralNet.py)

To submit on the cluster (using the slurm parameters in parameters.py), `--scan` must be replaced by two arguments
``` 
python HHMachineLearning.py (args) --submit name_of_jobs --split 1
```
`--submit` requires a string as name for the output dir (saved in `slurm`) 
`--split` requires the number of parameters used for each job (almost always 1)

*Note* : if using `--split N`, N! combinations will be used (might be reduncancies between different jobs) --- anyway, you will use 1 almost always
The split parameters will be saved in `split/` it is important that they remain there until the jobs jave finished running. After that they can be removed.

The output and logs will be in `slurm/name_of_jobs`

Now all the zip and csv files will be in the output directory but one needs to find the best one.

The first step is to concatenate the csv, to do that 
```
python HHMachineLearning.py --csv slurm/name_of_jobs/output/
```
This will create a concatenated csv file in model with name `name_of_jobs`, ordered according to the eval_criterion (evaluation error is better)
Note : for classification the F1 score is used and should be ordered in descending order (aka, the higher the better)

The easy way is then to pick the best model in the csv (it is ordered so it is easy), and get the corresponding zip file (also on the same line of the csv)
Let's say the best model is `slurm/name_of_jobs/output/one_of_the_job_output.zip`, to change its name one can use 
```
python Utils.py --zip slurm/name_of_jobs/output/one_of_the_job_output.zip model/my_model.zip
```

*Warning* : just changing the zip name will not work because the content also needs to change name (hence the function in `Utils.py`)

The other option with more details is to use the report option
```
python HHMachineLearning.py --report name_of_jobs
```
(the script automatically looks in `model` and adds the .csv extension, so you should not use it)

The script will then printout the 10 best models (according to the eval_criterion), plots on the console several histograms and produces png files.
The plot definitions are in plot_scans.py, they are [seaborn](https://seaborn.pydata.org/) based 

This will give clues on what parameters are doing better jobs. The zip file can the be dealt with the same way as before.

To then produce the test plots, one can use 
```
python HHMachineLearning.py (args) --model my_model --test
```
(same as csv, no need to specify the directory `model` nor .zip)

This will produce the root output files (split according to `split_name` in parameters.py) on the test set.

The plotting can be done in `Plotting/` (see associated README)

If other files have to be processed, one can use 
```
python HHMachineLearning.py (args) --model my_model --output key
``` 
Where `key` is one of the key in sampleList.py 

*Note* : There can be several keys 

*Warning* : these samples must not have been used in the training, this will cause undetected overfitting

... And that's it !!

#### Resubmission 
If some jobs failed, they can be resubmitted with the command 
```
python HHMachineLearning.py (args) --submit name_of_jobs --split 1 --resubmit /slurm/name_of_jobs/output/
```
The script will check what hyperparameters have been processed and which ones are missing, the corresponding jobs will be in a new directory and need to be moved back to the initial one before the csv concatenation step.

*Warning* : the hyperparameter dict in parameters.py must not change in the meantime !!! (especially number of epochs)
Otherwise the parameters in the csv will have changed. But the slurm parameters and keras callbacks can change at resubmission.

### Preprocessing and training/test split

What has not been dealt with in the previous sections is how the data preparation are handled.

#### Data split 

Depending on the ratios in parameters.py, a boolean mask is generated for each dataset.
    - False -> test set
    - True -> training set
The mask is generated as a npy object based on the suffix in parameters.py.

*Note* : If they do not exist, they will be generated and saved. If they exist they will just be loaded.
Warning : If the data changes, the code will exit with an error because the masks do not fit anymore (either delete them or change suffix).  

*Tip* : the point of the mask is that for each hyperparameter the training and test data will be the same and not randomized at each trial.

#### Preprocessing
Preprocessing is very important in machine learning to give all the features of the training the same importance.
We are using here the [Standard Scaler](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html), the point is to apply :
```
z -> z-mean/std
```
Where mean and std are the mean and standard deviation of the *training* data.

This scaler is saved in a pickle file with suffix in parameters.py as well (same tips and warnings as masks)

The easy way to use it is to transform the training and testing inputs, and do the inverse when saving into root files.
But keeping track of both model and scaler is annoying...

So a custom layer in preprocessing.py incoorporates the mean and std as weights that are then saved in the model. No need to keep track of the scaler anymore when sharing the model.
On the ther side when loading the model, the script must be given so that Keras knows how to handle it (but already included in the machinery here).

### Learning weights

In order to represent in the training the physical significance of the training events, the event weight needs to be used ([doc](https://keras.io/guides/customizing_what_happens_in_fit/#supporting-sampleweight-amp-classweight)).

This is what is given as weight in parameters.py, the issue arises from the negative weights. They can be dealt with in several scenarios
- Use the negative weights as they are : physically makes sens, but if tone batch size contains mostly negative weights the learning will go in the wrong direction.
- Use the absolute value : all statistics, but repressed phase space regions will have more importance in the training.
- Use only the positive values (aka, cut on MC weigt > 0) : reduces the statistics a bit but avois the problem

We use the last option.

This was the events in one sample are correctly weighted wrt to each other. To weigh between samples, the learning weights is set as :
```
learning weight = event weight * Xsec / event weight sum
```
Where these values come from the json files in parameters.py

*Warning* : this is only valid for backgrounds, for signal the Xsec is unknown so better keep only the event weight

On the other side, it is possible that there is less signal statistics than background. To alleviates that, the sum of learning weights is equalized between signal and background.

Eg: learning weights (signal) /= sum(learning weights (signal)) and same for background.

In case of multiclassification (eg, ST, DY, and TT classes) all classes need to have the same sum of learning weights

This is implemented in the beginning of HHMachineLearning.py.

### Cross-validation
[Cross-validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)) is a very useful technique for two reasons :
    - It makes full use of the statistics (reduces the bias when we have a small number of events)
    - Allows all the events to be passed through the model when filling histograms in HEP
In the classical approcah (also called holdout), we train on say ~70% of data and evaluate the performances on the remaining 30% which is unbiased because not seen by the model during the training.
The issue is that these 30% of the dataset are not used for training the network, if the dataset is huge it is not a problem but otherwise we are impairing the training.

To overcome that, cross-validation consists in training several models on subparts of the data : training, evaluation and application. The first one should be the largest and the last one must not be seen at any point by the associated network here. This will be the events evaluated by it only at the analysis level when filling plots.

*Warning* : it implies that each hyperparameter set will be trained several times, take that into account when evaluating computing time.

Example : the data is split in 5 slices and we train 5 models. We will split according to the event number (because it is unique) but any other branch in the tree will work.
For a model i :
    - it will be applied on events for which event_number % 5 == i
    - it will be evaluated on events for which event_number % 5 == (i+1)%5
    - it will be trained on the remaining events
Another possiblity is 3 networks for 6 slices, therefore each network is applied on two data slices

You can also do some fancy stuff like 5 models on 15 slices (10 for learning, 2 for evaluation and 3 for application)

*Warning* : the assumption made here is that each event can only be applied on one model, hence it is required that the application number is equal to number of slices / number of models.
In addition the number of slices must be multiple of the number of models. Several assertions in the code are there to enforce it.

To enable cross-validation, change the boolean `crossvalidation` in `parameters.py` and tweak the numbers below accordingly.
The result is that a mask will be added in the dataframe and will serve to determine which model to train/evaluate.

*Note* :  the mask is now saved in the dataframe, and since it depends on the event number (or another metric in the branches) that is invariant there is no need to save it anymore so no mask is saved as npy file (the scaler will still be).

Apart from that the training will remain the same, except that now several models and csv file will be saved, one for each step with additional string _crossval%d. 

Csv concatenation works the same, however in principle the hyperparameter scan with cross-validation would require two nested loops. This is not implemented (maybe in the future) so right now the user has to select one hyperparameter set that has several submodels (aka with same substring "_dict%d" but different "_crossval%d") - we will use submodels to qualify the different models coming from one cross-validation pass. All these submodels need to be saved in `model/` (possibly using `python Utils.py --zip [] []`).

To testing the command for output becomes 
```
python HHMachineLearning.py (args) --test --model model_crossval0 model_crossval1 model_crossval2 model_crossval3 model_crossval4
```

*Note* : `--model` accepts a list as above and will interpret the models as going in that order ! It is the user responsibility to make sure they are given in the right order otherwise models will not be applied on the correct dataset slice (and no error will be shown). 

Each model will then be evaluated on its application slice (aka the one that should be used at the analysis level) to see what is its behaviour and the output tree will be produced as in the classical approach. 

### Generator
In case there is too much data in the training to put them in the RAM, small chunks can be loaded in turns and trained on.
The advantage is that many threads can be used to generate the training data from root files.

To use the generator, add the flag `--generator`. 

*Note* : no data importation is done before the learning, everything is done on the fly which introduces some complexity

The generator can use both classic holdout or cross-validation but the working differs from before.

*Warning* : the scaler still needs to be computed before the learning can actually take place. If no mask is found (based on suffix in `parameters.py`) it will create one by loopng through the files in the sample list YAML file batch by batch (same variables as for output).

The generator tries to be smart and populate the batches in the same proportion as in the input root files (to avoid a batch being filled only with the same class) and will reweight each class to not imbalance the training. The repartition will be printed out for the user to check everything seems fine. 

*Note* : for this to work well, the batch size needs to be large enough. Given also that the bootleneck comes from the data transfer, large batches will improve computation time, especially for the GPU (the IO bandwidth there is not amazing). 

Keras provides the possibility to use multiprocessing, this is specified by the `workers` variable in `parameters.py`. 0 means the main thread will go fetch a batch from the input root files and then do the learning, sequentially. 1 or more means one thread does the learning and n others will produce batches and put them into queue.

This means that when submitting jobs, additional number of tasks can be asked to accelerate training (see the slurm parameters in `parameters.py`).

Producing output trees will also work in generator mode, except that the workers part is not available (otherwise could not retrieve the inputs to put in output tree). So the production is done in a simple loop (will look later into multiprocessing for that). One root tree will be produced per slice, a simple `hadd` afterwards will do just fine.

#### Holdout mode
This will run when `parameters.crossvalidation = False`.

In this mode, a mask still needs to be computed (with keys 0 for training, 1 for evaluation and 2 for output) : this is done automatically based on `parameters.suffix` (generated or loaded if present).

Apart from that the usual commands will still work

*Note* : given technical issues, the imported dataframe size will be bigger than the batch size (sometimes significantly) and an internal mask will be used to select the events in the corretct set. This will be fixed in the future.

#### Cross-validation

This will run when `parameters.crossvalidation = True`.

Nothing much new compared to usual commands. The output will be produced per slice AND per model (a bit slower then, but `hadd` will also work as fine).


### Cache
The importation from root files can be slow and if the training data is not too big it can be cached (see name in parameters.py).
This is especially useful when modifying the code for rapid testing or evaluation on local. For job submission it is best to use `--nocache`

*Warning* : whenever you change something in sampleList.py, the preprocessing or mask, the cache must be cleared !!!
Otherwise you will still run on the older cache values and not the changes you chose.

*Note* : caching is of course not used if using the generator.

## Authors

* **Florian Bury** -- [Github](https://github.com/FlorianBury)

## Acknowledgments

