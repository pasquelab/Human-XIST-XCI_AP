## Step-by-Step Instructions

Follow these steps to create and activate the conda environment for the XIST loci analysis project.

### Step 1: Create the Conda Environment
- Download the Human-XIST-XCI_AP/SCRIPT/Image Analysis/3D-SIM/environment.yml file

Note: Make sure you have Miniconda or Anaconda installed.

- Run the following command in your Anaconda terminal:

```
conda env create -f PATH/environment.yml
```

This will create a new environment named xist-analysis using the specifications in the environment.yml file.

### Step 2: Activate the Environment
Once the environment is created, activate it using:

```
conda activate xist-analysis
```

### Step 3 (Optional): Deactivate the Environment

```
conda deactivate
```

You are now ready to run the Python scripts for XIST loci analysis.

