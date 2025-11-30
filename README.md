# Preprocessing for Image Mass Cytometry dataset
## Updated: November 29, 2025. By Tina (Yi-Ting) Hsu 
### Background and Rationale
**IMC_preprocess** is a workflow for preprocessing [Imaging Mass Cytometry™](https://www.fluidigm.com/applications/imaging-mass-cytometry) (IMC) dataset. IMC is a spatial modality with single-cell proteiomics expression data as well as spatial context represented by x-y coordinates of all the cells in images. In the analysis work of IMC, we usually receive the dataset in batches due to the generation time (it is commonly 5-6 hours per core). Hence, it costs more and more time on repeating preprocessing steps since we need to change all the parameters everytime for the new dataset. To adress this, this pipeline helps execute all the preprocessing works before cell type annotation in an automatic way which also establishes standardization in the projet. This workflow includes the typical preprocessing steps in the field, such as quality control and clustering, and the specific procedure of the workflow is documented below. Other single cell expression data with changes to the input data format is available as well for the pipeline which is described below in detail. 

The underlying approach to this workflow has 3 modules:
- Batch Effect Checking: Show the batch effect situation by [UMAP plot](https://umap-learn.readthedocs.io/en/latest/). 
- Leiden Clustering: Run Leiden Clustering to make all the cells in clusters/categories.
- Single-cell Protein Expression across Clusters: Produce a protein expression heatmap of each cell across clusters.

<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/DAG.png?raw=true" width="300" height="300">

Package Dependencies:
*You should have Python 3 and Miniconda installed for Linux*
  - nextflow
  - matplotlib = 3.3.3
  - numpy = 1.19.4
  - pandas = 1.1.4
  - seaborn = 0.11.0
  - pip = 20.2.4

### Usage
The execution of the package *requires* a **Linux OS**. If you are using Windows, you can also run it by **wsl**.

1. Clone this file in a directory you prefer. This will download all the files needed for running the pipeline. The output will be stored in the folder created as "output".
```
git clone https://github.com/
```
2. Create an environment with all the dependencies. 
```
conda env create --name nf_imc_preprocessing --file env
```
3. Activate the environment.
```
conda activate nf_imc_preprocessing
```
4. Now run IMC_preprocess by executing the following code. This will create two files in the **output** folder: `batch_umap.png`, `new_adata.h5ad`, `cluster_result.png`.
```
nextflow run workflow.nf
```
5. Generate a DAG from the workflow:
```
snakemake --dag results/heatmap.png | dot -Tsvg > dag.svg
```

### Input
IMC_preprocessing takes in one h5ad file: an anndata-structure file with all quantified spatial data. All the files should be placed in the folder `data/`. 

Inside `data/`, the sample dataset we use here is a segmented and quantified IMC dataset with 4 breast cancer patients represented. All protein expressions and x-y coordinates of per cell should be stored in the h5ad.

H5AD should have the following properties:
*If you are using the spatial data from other technologies, please make sure the data type/format is h5ad and the following properties stored in the correct layers in an anndata object. *

- anndata.X: expression values of all the markers
- anndata.obsm["spatial"]: x-y coordinates of all the cells across iamges.
- anndata.obs["unique_core_id"]: the column in obs with the unique key to the cores or images.

### Output
IMC_preprocessing produces two types of visualizations in png files for interpreting data characteristics in preprocessing step: UMAP plot, expression heatmap plot, and one h5ad file with the clustering results. All the files should be stored in the folder `output/`.

Showcase for the output from IMC_preprocessing 
<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/output/batch_umap.png?raw=true" width="300" height="300">

- **batch_umap.png**: This is a visualization for high dimensional data based on Uniform Manifold Approximation and Projection. Axes do not refer to spatial coordinates. A UMAP plot is included to show that how different are the images, cores or batches. This is a common step in quality control that to evaluate if batch correction is needed or not. If there's no series batch effect across batches, then we can run the clustering with the whole adata for annotation.

------------------------------------------------------

- **new_adata.h5ad**: an anndata object with a new column, leiden in obs, which represents the labels/clusters from Leiden Clustering.

------------------------------------------------------
<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/output/cluster_result.png?raw=true" width="300" height="300">

- **cluster_result.png**: This graph shows the clustering results by the marker expressions ofr each cell across clusters. The objective of the plot is to help finish cell type annotation.For each cluster, we can give the cell type based on the specific higher expressions in the corresponding protein markers/ For example, if most of the cells in Cluster A highly express in CD20 that we may annotate the cluster as B Cell cluster.




11月-29 22:36:34.258 [main] DEBUG nextflow.cli.Launcher - $> nextflow run workflow.nf -with-docker 'tina0420/imc_preprocessing:latest' --work-dir output
11月-29 22:36:34.305 [main] DEBUG nextflow.cli.CmdRun - N E X T F L O W  ~  version 25.10.2
11月-29 22:36:34.318 [main] DEBUG nextflow.plugin.PluginsFacade - Setting up plugin manager > mode=prod; embedded=false; plugins-dir=/Users/Eric/.nextflow/plugins; core-plugins: nf-amazon@3.4.2,nf-azure@1.20.2,nf-cloudcache@0.5.0,nf-codecommit@0.5.0,nf-console@1.3.0,nf-google@1.23.3,nf-k8s@1.2.2,nf-tower@1.17.3,nf-wave@1.16.1
11月-29 22:36:34.331 [main] INFO  o.pf4j.DefaultPluginStatusProvider - Enabled plugins: []
11月-29 22:36:34.331 [main] INFO  o.pf4j.DefaultPluginStatusProvider - Disabled plugins: []
11月-29 22:36:34.333 [main] INFO  org.pf4j.DefaultPluginManager - PF4J version 3.12.0 in 'deployment' mode
11月-29 22:36:34.340 [main] DEBUG nextflow.util.RetryConfig - Missing nextflow session - using default retry config
11月-29 22:36:34.410 [main] INFO  org.pf4j.AbstractPluginManager - No plugins
11月-29 22:36:34.425 [main] DEBUG nextflow.config.ConfigBuilder - Found config local: /Users/Eric/Downloads/IMC_preprocessing/nextflow.config
11月-29 22:36:34.426 [main] DEBUG nextflow.config.ConfigBuilder - Parsing config file: /Users/Eric/Downloads/IMC_preprocessing/nextflow.config
11月-29 22:36:34.440 [main] DEBUG nextflow.config.ConfigBuilder - Applying config profile: `standard`
11月-29 22:36:34.657 [main] DEBUG nextflow.config.ConfigBuilder - Enabling execution in Docker container as requested by command-line option `-with-docker tina0420/imc_preprocessing:latest`
11月-29 22:36:34.660 [main] DEBUG nextflow.plugin.PluginsFacade - Plugins default=[]
11月-29 22:36:34.660 [main] DEBUG nextflow.plugin.PluginsFacade - Plugins resolved requirement=[]
11月-29 22:36:34.667 [main] DEBUG n.secret.LocalSecretsProvider - Secrets store: /Users/Eric/.nextflow/secrets/store.json
11月-29 22:36:34.668 [main] DEBUG nextflow.secret.SecretsLoader - Discovered secrets providers: [nextflow.secret.LocalSecretsProvider@536d97f8] - activable => nextflow.secret.LocalSecretsProvider@536d97f8
11月-29 22:36:34.669 [main] DEBUG nextflow.cli.CmdRun - Applied DSL=2 by global default
11月-29 22:36:34.675 [main] DEBUG nextflow.cli.CmdRun - Launching `workflow.nf` [irreverent_watson] DSL2 - revision: 9e7d406abe
11月-29 22:36:34.705 [main] DEBUG nextflow.Session - Session UUID: b63b14ac-2b29-4937-8280-9c42b00a82e9
11月-29 22:36:34.706 [main] DEBUG nextflow.Session - Run name: irreverent_watson
11月-29 22:36:34.706 [main] DEBUG nextflow.Session - Executor pool size: 8
11月-29 22:36:34.709 [main] DEBUG nextflow.file.FilePorter - File porter settings maxRetries=3; maxTransfers=50; pollTimeout=null
11月-29 22:36:34.712 [main] DEBUG nextflow.util.ThreadPoolBuilder - Creating thread pool 'FileTransfer' minSize=10; maxSize=24; workQueue=LinkedBlockingQueue[-1]; allowCoreThreadTimeout=false
11月-29 22:36:34.723 [main] DEBUG nextflow.cli.CmdRun - 
  Version: 25.10.2 build 10555
  Created: 28-11-2025 19:24 UTC (11:24 PDT)
  System: Mac OS X 15.4.1
  Runtime: Groovy 4.0.28 on Java HotSpot(TM) 64-Bit Server VM 25.0.1+8-LTS-27
  Encoding: UTF-8 (UTF-8)
  Process: 36224@xuyitingdeMacBook-Air.local [127.0.0.1]
  CPUs: 8 - Mem: 16 GB (1.5 GB) - Swap: 0 (0)
11月-29 22:36:34.729 [main] DEBUG nextflow.Session - Work-dir: /Users/Eric/Downloads/IMC_preprocessing/work [Mac OS X]
11月-29 22:36:34.729 [main] DEBUG nextflow.Session - Script base path does not exist or is not a directory: /Users/Eric/Downloads/IMC_preprocessing/bin
11月-29 22:36:34.734 [main] DEBUG nextflow.executor.ExecutorFactory - Extension executors providers=[]
11月-29 22:36:34.737 [main] DEBUG nextflow.Session - Observer factory (v2): LinObserverFactory
11月-29 22:36:34.738 [main] DEBUG nextflow.Session - Observer factory (v2): DefaultObserverFactory
11月-29 22:36:34.753 [main] DEBUG nextflow.cache.CacheFactory - Using Nextflow cache factory: nextflow.cache.DefaultCacheFactory
11月-29 22:36:34.757 [main] DEBUG nextflow.util.CustomThreadPool - Creating default thread pool > poolSize: 9; maxThreads: 1000
11月-29 22:36:34.783 [main] DEBUG nextflow.Session - Session start
11月-29 22:36:34.897 [main] DEBUG nextflow.script.ScriptRunner - > Launching execution
11月-29 22:36:34.935 [main] DEBUG n.script.dsl.ProcessConfigBuilder - Config settings `withName:batch_effect_check_plot` matches process batch_effect_check_plot
11月-29 22:36:34.943 [main] DEBUG nextflow.executor.ExecutorFactory - << taskConfig executor: null
11月-29 22:36:34.944 [main] DEBUG nextflow.executor.ExecutorFactory - >> processorType: 'local'
11月-29 22:36:34.948 [main] DEBUG nextflow.executor.Executor - [warm up] executor > local
11月-29 22:36:34.951 [main] DEBUG n.processor.LocalPollingMonitor - Creating local task monitor for executor 'local' > cpus=8; memory=16 GB; capacity=8; pollInterval=100ms; dumpInterval=5m
11月-29 22:36:34.952 [main] DEBUG n.processor.TaskPollingMonitor - >>> barrier register (monitor: local)
11月-29 22:36:34.962 [main] DEBUG nextflow.processor.TaskProcessor - Creating process 'batch_effect_check_plot': maxForks=0; fair=false; array=0
11月-29 22:36:34.992 [main] DEBUG n.script.dsl.ProcessConfigBuilder - Config settings `withName:run_leiden_clustering` matches process run_leiden_clustering
11月-29 22:36:34.993 [main] DEBUG nextflow.executor.ExecutorFactory - << taskConfig executor: null
11月-29 22:36:34.993 [main] DEBUG nextflow.executor.ExecutorFactory - >> processorType: 'local'
11月-29 22:36:34.993 [main] DEBUG nextflow.processor.TaskProcessor - Creating process 'run_leiden_clustering': maxForks=0; fair=false; array=0
11月-29 22:36:34.994 [main] DEBUG n.script.dsl.ProcessConfigBuilder - Config settings `withName:plot_cluster_heatmap` matches process plot_cluster_heatmap
11月-29 22:36:34.995 [main] DEBUG nextflow.executor.ExecutorFactory - << taskConfig executor: null
11月-29 22:36:34.995 [main] DEBUG nextflow.executor.ExecutorFactory - >> processorType: 'local'
11月-29 22:36:34.995 [main] DEBUG nextflow.processor.TaskProcessor - Creating process 'plot_cluster_heatmap': maxForks=0; fair=false; array=0
11月-29 22:36:34.997 [main] DEBUG nextflow.Session - Workflow process names [dsl2]: run_leiden_clustering, batch_effect_check_plot, plot_cluster_heatmap
11月-29 22:36:34.998 [main] DEBUG nextflow.Session - Igniting dataflow network (3)
11月-29 22:36:34.998 [main] DEBUG nextflow.processor.TaskProcessor - Starting process > batch_effect_check_plot
11月-29 22:36:34.999 [main] DEBUG nextflow.processor.TaskProcessor - Starting process > run_leiden_clustering
11月-29 22:36:34.999 [main] DEBUG nextflow.processor.TaskProcessor - Starting process > plot_cluster_heatmap
11月-29 22:36:34.999 [main] DEBUG nextflow.script.ScriptRunner - Parsed script files:
  Script_88f86be359e9b80c: /Users/Eric/Downloads/IMC_preprocessing/workflow.nf
11月-29 22:36:34.999 [main] DEBUG nextflow.script.ScriptRunner - > Awaiting termination 
11月-29 22:36:34.999 [main] DEBUG nextflow.Session - Session await
11月-29 22:36:35.087 [Task submitter] DEBUG n.executor.local.LocalTaskHandler - Launch cmd line: /bin/bash -ue .command.run
11月-29 22:36:35.088 [Task submitter] INFO  nextflow.Session - [87/da150b] Submitted process > batch_effect_check_plot (Batch effect check plot)
11月-29 22:36:35.094 [Task submitter] DEBUG n.executor.local.LocalTaskHandler - Launch cmd line: /bin/bash -ue .command.run
11月-29 22:36:35.095 [Task submitter] INFO  nextflow.Session - [d1/a17753] Submitted process > run_leiden_clustering (Run leiden clustering)
11月-29 22:36:35.429 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 2; name: run_leiden_clustering (Run leiden clustering); status: COMPLETED; exit: 127; error: -; workDir: /Users/Eric/Downloads/IMC_preprocessing/work/d1/a1775313230cb5722d8c2576841d1f]
11月-29 22:36:35.430 [Task monitor] DEBUG nextflow.util.ThreadPoolBuilder - Creating thread pool 'TaskFinalizer' minSize=10; maxSize=24; workQueue=LinkedBlockingQueue[-1]; allowCoreThreadTimeout=false
11月-29 22:36:35.435 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 1; name: batch_effect_check_plot (Batch effect check plot); status: COMPLETED; exit: 2; error: -; workDir: /Users/Eric/Downloads/IMC_preprocessing/work/87/da150b9254db9af92a6812b9a64456]
11月-29 22:36:35.436 [TaskFinalizer-1] DEBUG nextflow.processor.TaskProcessor - Handling unexpected condition for
  task: name=run_leiden_clustering (Run leiden clustering); work-dir=/Users/Eric/Downloads/IMC_preprocessing/work/d1/a1775313230cb5722d8c2576841d1f
  error [nextflow.exception.ProcessFailedException]: Process `run_leiden_clustering (Run leiden clustering)` terminated with an error exit status (127)
11月-29 22:36:35.437 [TaskFinalizer-2] DEBUG nextflow.processor.TaskProcessor - Handling unexpected condition for
  task: name=batch_effect_check_plot (Batch effect check plot); work-dir=/Users/Eric/Downloads/IMC_preprocessing/work/87/da150b9254db9af92a6812b9a64456
  error [nextflow.exception.ProcessFailedException]: Process `batch_effect_check_plot (Batch effect check plot)` terminated with an error exit status (2)
11月-29 22:36:35.454 [TaskFinalizer-2] ERROR nextflow.processor.TaskProcessor - Error executing process > 'batch_effect_check_plot (Batch effect check plot)'

Caused by:
  Process `batch_effect_check_plot (Batch effect check plot)` terminated with an error exit status (2)


Command executed:

  python batch_effect_check_plot.py     --input biof501_data.h5ad     --figdir .

Command exit status:
  2

Command output:
  (empty)

Command error:
  python: can't open file 'batch_effect_check_plot.py': [Errno 2] No such file or directory

Work dir:
  /Users/Eric/Downloads/IMC_preprocessing/work/87/da150b9254db9af92a6812b9a64456

Container:
  tina0420/imc_preprocessing:latest

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
11月-29 22:36:35.459 [TaskFinalizer-2] DEBUG nextflow.Session - Session aborted -- Cause: Process `batch_effect_check_plot (Batch effect check plot)` terminated with an error exit status (2)
11月-29 22:36:35.459 [main] DEBUG nextflow.Session - Session await > all processes finished
11月-29 22:36:35.499 [TaskFinalizer-2] INFO  nextflow.Nextflow - Error in the workflow!
11月-29 22:36:35.501 [TaskFinalizer-2] DEBUG n.trace.WorkflowStatsObserver - Workflow completed > WorkflowStats[succeededCount=0; failedCount=2; ignoredCount=0; cachedCount=0; pendingCount=0; submittedCount=0; runningCount=0; retriesCount=0; abortedCount=0; succeedDuration=0ms; failedDuration=533ms; cachedDuration=0ms;loadCpus=0; loadMemory=0; peakRunning=2; peakCpus=2; peakMemory=0; ]
11月-29 22:36:35.538 [Task monitor] DEBUG n.processor.TaskPollingMonitor - <<< barrier arrives (monitor: local) - terminating tasks monitor poll loop
11月-29 22:36:35.538 [main] DEBUG nextflow.Session - Session await > all barriers passed
11月-29 22:36:35.540 [main] DEBUG n.trace.WorkflowStatsObserver - Workflow completed > WorkflowStats[succeededCount=0; failedCount=2; ignoredCount=0; cachedCount=0; pendingCount=0; submittedCount=0; runningCount=0; retriesCount=0; abortedCount=0; succeedDuration=0ms; failedDuration=533ms; cachedDuration=0ms;loadCpus=0; loadMemory=0; peakRunning=2; peakCpus=2; peakMemory=0; ]
11月-29 22:36:35.540 [TaskFinalizer-2] DEBUG nextflow.Session - null
java.lang.InterruptedException: null
	at java.base/java.lang.Object.wait0(Native Method)
	at java.base/java.lang.Object.wait(Object.java:389)
	at java.base/java.lang.Thread.join(Thread.java:1887)
	at java.base/java.lang.Thread.join(Thread.java:1963)
	at nextflow.trace.AnsiLogObserver.onFlowComplete(AnsiLogObserver.groovy:520)
	at nextflow.Session$_notifyFlowComplete_lambda39.doCall(Session.groovy:1121)
	at nextflow.Session.notifyEvent(Session.groovy:1151)
	at nextflow.Session.notifyFlowComplete(Session.groovy:1121)
	at nextflow.Session.shutdown0(Session.groovy:779)
	at nextflow.Session.abort(Session.groovy:834)
	at nextflow.Session.fault(Session.groovy:796)
	at nextflow.processor.TaskPollingMonitor.finalizeTask(TaskPollingMonitor.groovy:753)
	at nextflow.processor.TaskPollingMonitor.safeFinalizeTask(TaskPollingMonitor.groovy:736)
	at java.base/jdk.internal.reflect.DirectMethodHandleAccessor.invoke(DirectMethodHandleAccessor.java:104)
	at java.base/java.lang.reflect.Method.invoke(Method.java:565)
	at org.codehaus.groovy.reflection.CachedMethod.invoke(CachedMethod.java:343)
	at groovy.lang.MetaMethod.doMethodInvoke(MetaMethod.java:328)
	at groovy.lang.MetaClassImpl.doInvokeMethod(MetaClassImpl.java:1339)
	at groovy.lang.MetaClassImpl.invokeMethod(MetaClassImpl.java:1094)
	at groovy.lang.MetaClassImpl.invokeMethod(MetaClassImpl.java:1007)
	at org.codehaus.groovy.runtime.InvokerHelper.invokePogoMethod(InvokerHelper.java:645)
	at org.codehaus.groovy.runtime.InvokerHelper.invokeMethod(InvokerHelper.java:628)
	at org.codehaus.groovy.runtime.InvokerHelper.invokeMethodSafe(InvokerHelper.java:82)
	at nextflow.processor.TaskPollingMonitor$_checkTaskStatus_lambda8.doCall(TaskPollingMonitor.groovy:726)
	at java.base/java.util.concurrent.Executors$RunnableAdapter.call(Executors.java:545)
	at java.base/java.util.concurrent.FutureTask.run(FutureTask.java:328)
	at java.base/java.util.concurrent.ThreadPoolExecutor.runWorker(ThreadPoolExecutor.java:1090)
	at java.base/java.util.concurrent.ThreadPoolExecutor$Worker.run(ThreadPoolExecutor.java:614)
	at java.base/java.lang.Thread.run(Thread.java:1474)
11月-29 22:36:35.670 [main] DEBUG nextflow.cache.CacheDB - Closing CacheDB done
11月-29 22:36:35.676 [main] DEBUG nextflow.script.ScriptRunner - > Execution complete -- Goodbye



