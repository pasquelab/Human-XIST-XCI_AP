
# 3D-SIM XIST foci analysis

This section contains the scripts that were used for the quantification of 3D-SIM acquired images of cells (**H9 Naive, EXMCs and TSCs**).

## Imaris-Based Segmentation and Detection

For the segmentation and detection of nuclei and XIST foci, we used **Imaris (Oxford Instruments) version 10.2.0**. The following steps outline the preprocessing workflow:

1. **Download the Parameter File**  
   Retrieve the predefined segmentation parameters from:  
   `Human-XIST-XCI_AP/SCRIPT/Image Analysis/3D-SIM/imarisv10_2_0_cell_segmentation_creationparameteres.icpx`

2. **Import into Imaris**  
   Open Imaris and import the `.icpx` file via the Cell Segmentation module. Adjust the parameters as needed to suit your specific dataset.

3. **Segmentation Configuration**  
   The following settings were used in this study:

   - **Application Type**: 3D  
   - **Cell Type**: Cell Body

   **Nucleus Detection**  
   - Enabled  
   - Estimated Diameter: 2 µm  
   - Background Subtraction: Enabled  
   - Manual Threshold: 91.4578

   **Cell Detection**  
   - Enabled  
   - Estimated Diameter: 10 µm  
   - Background Subtraction: Enabled  
   - Manual Threshold: 20

   **Vesicle Detection**  
   - Enabled  
   - Vesicle Type: Type A  
   - Estimated Diameter: 0.15 µm  
   - Background Subtraction: Enabled  
   - Manual Threshold: 240

   **Segmentation Options**  
   - Split Cells: By Seeds  
   - Split Cells by Masks: Voronoi  
   - Fill Holes: Enabled (for both cell bodies and nuclei)  
   - Erase Nucleus from Cell Body: Enabled

   **Filtering Criteria**  
   - Cell Volume: ≥ 5.85 × 10⁶ voxels  
   - Nucleus Volume: ≥ 108,517 voxels  
   - Cell Seed Quality: ≥ 26.6642

4. **Export Data for Analysis**  
   Export the following CSV files from Imaris for downstream analysis:
   - `Cell_Number_Of_Vesicles_VesicleType=Vesicles_Type_A.csv`: Number of XIST foci per nucleus
   - `Nucleus_Volume.csv`: Segmented nucleus volumes (used for normalization)
   - `Vesicles_Type_A_Position.csv`: 3D centroid coordinates of XIST foci

These exported datasets serve as the input for the Python-based nearest neighbor distance analysis and the R-based visualization and statistical analyses described in the next sections.


## Python-Based Nearest Neighbor Analysis

This Python script calculates the nearest neighbor distances between vesicle positions in 3D space, grouped by cell and image. It outputs both individual neighbor distances and their median values per image between an N-number of neighbors.

### Input

A CSV file containing vesicle position data with the following columns:

- `Original Image Name`
- `CellID`
- `ID`
- `Vesicles Type A Position X`
- `Vesicles Type A Position Y`
- `Vesicles Type A Position Z`

### Parameters

- `file_path`: Path to the input CSV file (default: `/content/Vesicles_Type_A_Position.csv`)
- `neighbors`: Number of nearest neighbors to compute (default: `5`)

### Functionality

#### `calculate_nearest_neighbor(...)`

Calculates the nearest neighbor distances for each vesicle within the same cell and image.

**Parameters:**

- `df`: Input DataFrame
- `filename_column`, `cell_id_column`, `x_column`, `y_column`, `z_column`, `id_column`: Column names
- `neighbors`: Number of neighbors to compute

**Returns:**

- A DataFrame with median nearest neighbor distances for each vesicle.

### Data Processing Steps

1. Load the input `Vesicles_Type_A_Position.csv` file.
2. Filter out rows with missing coordinate values.
3. Group data by `Original Image Name` and `CellID`.
4. Use a KD-Tree to compute nearest neighbor distances.
5. Save:
   - All nearest neighbor distances per foci to `nearest_neighbor_distances.csv`
   - Median distances per image to `median_nearest_neighbor_distances.csv`

### Output

- `nearest_neighbor_distances.csv`: Contains distances to the nearest neighbors for each vesicle.
- `median_nearest_neighbor_distances.csv`: Contains median nearest neighbor distances per image.

## R-Based Data Visualization and Statistical Analysis

### Step 1: Count cells per image per cell type that have "foci out"
#### •	Goal: 
Identify cells where NearestNeighborDistance_5 > 4. We choose the threshold 4 as it the value where all the EXMCs and TSCs still had lots of foci out whereas all naïve cells were within.
•	Groups by image, cell ID, and cell type → determines if that cell has "foci out" (1) or not (0).
•	Keeps only those with foci_out == 1.
•	Then counts how many such cells exist per image and cell type.

### Step 2: Count total unique cells per cell type
•	Get the total number of unique cells per cell type (regardless of "foci out").

### Step 3: Sum the number of "foci out" cells per cell type across all images

### Step 4: Merge with total cells per cell type
•	Combine the count of "foci out" cells with the total cell count per type.

### Step 5: Plot the percentage of "foci out" cells per cell type
•	Create a bar plot showing the percentage of cells with NearestNeighborDistance_5 > 4 per cell type.

### Step 6: Create a vector called test_data_tsc_exmc, combining cell counts from images annotated with "H9 ESC TSCs" and "H9 ESC EXMCs". 
1.	The number of TSC cells.
2.	The non-TSC cells in those images (total_cells - num_cells_total).
3.	The number of EXMC cells.
4.	The non-EXMC cells in those images.
This structure is useful for testing enrichment function which we need to use below. Similarily, create vector for Naïve and TSCs cells.

### Step 7: Perform statistics using enrichment(data) function

#### Input: A 4-element vector (e.g., from test_data_h9naive_exmc or test_data_tsc_exmc).

#### Output:
1.	fisher.test(df) → Fisher's Exact Test result.
2.	chisq.test(df) → Chi-squared Test result.
3.	A data.table formatted for plotting, with calculated percentages of “Yes” and “No” per group:
o	Reshaped using dcast and melt
o	Annotated with genes_factor for barplot faceting

### Step 8: This step reshapes your data and creates a stacked bar plot showing the percentage of Spreaded vs. Non-Spreaded cells.
•	Add a column non_spreaded, the complement of num_cells_total.
•	Reshape the data into long format with:
o	XIST_cloud = "Spreaded" or "Non-Spreaded"
o	cell_count = actual count
•	This enables plotting stacked bars for each condition.
•	Plot- Shows percentage of Spreaded vs. Non-Spreaded cells per cell_type, Stacked bars add to 100% for each cell type, Colors: "red" for Spreaded, "gray" for Non-Spreaded, Minimalist theme for clean visuals

## To calculate the number of foci per cell:
•	group_by(cell_type): Group the data by each cell type.
•	summarise(unique_foci = unique(...)): For each group, extract all unique values of number_of_foci_per_cell.
•	ungroup(): Remove the grouping afterward.

#### Output: Combined violin + box plot with overlayed points. It should produce a nice visualization with:
•	violin plots showing distribution shape and density,
•	narrow boxplots for median,
•	and individual points overlaid.
