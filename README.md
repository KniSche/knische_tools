# KniSche tools

> A selection of custom functions for bioinformatics analysis centred on single cell and spatial next generation sequencing.

## Installation

You can install this list of tools as a package directly from GitHub using pip:

```bash
pip install git+https://github.com/knische/knische_tools
```

## List of Functions
1. `plot_global_spatial` 
	This tool is plot wrapper function around `sc.pl.embedding`, allowing users to plot subset embeddings
	from `adata.obsm["embedding"]` onto the embedding made from the entire set.
	Primary use is solving the problem of visualizing a small subset of cells 
	(e.g., a specific tissue section or cell type) while maintaining the context of the full dataset's embedding
	without requiring the user to manually manage massive coordinate files.
	<details>
	<summary>Usage and Documentation</summary>
	### Prerequisites
	This tool assumes you are using Scanpy and have a global reference embedding in adata.uns
	(e.g., a CSV of X/Y coordinates for the whole atlas). The function looks for `adata.uns['global_spatial']` by default.
		
	This can be set up as follows:
	```python
	import pandas as pd
	import scanpy as sc
	# load the adata	
	adata = sc.read_h5ad("my_subset_data.h5ad")
	
	# Load the Master/Global coordinates (Index must be cell barcodes)
	# Format: DataFrame formatted as [barcodes, x, y], with index=barcodes and columns=[x, y]
	global_df = pd.read_csv("path/to/global_coords.csv", index_col=0)
	
	# Store in .uns
	adata.uns['global_spatial'] = global_df

	```
	### 
	
	### Usage
	```python
	import matplotlib.pyplot as plt
	from knische_tools import plot_global_spatial

	# Simple plot
	plot_global_spatial(
		adata, 
		color=['leiden', 'donor']
	)

	plt.show()
	```
	
	### Additional Documentation
	Plots a subset of data on top of the full atlas coordinates.

	**Parameters:**
	
	* **`adata`** *(AnnData)*: The subset object. Must contain global coords in `adata.uns[uns_key]`.
	* **`color`** *(str or list)*: Keys for annotations of observations/cells or variables/genes (e.g., `'cell_type'` or `['geneA', 'geneB']`).
	* **`basis`** *(str, default='global_spatial')*: The key to store/access the aligned coordinates in `adata.obsm`.
	* **`uns_key`** *(str, default='global_spatial')*: The key in `adata.uns` where the full reference DataFrame is stored.
	* **`size`** *(int, default=10, optional)*: Size of the subset points. If `None`, defaults to Scanpy's automatic sizing.
    * **`**kwargs`**: Additional arguments passed to `sc.pl.embedding` (e.g., `cmap`, `vmax`, `frameon`, `alpha`).
	
	**Background underlay options:**
	* **`background_size`** *(int, default=0.0002)*: Size of the global reference points.
    * **`background_color`** *(str, default='lightgrey')*: Color of the **global reference points**. *Tip: If using a dark figure_color, try a darker grey (e.g., '#333333') for a subtle effect.*
	
	**Style & Theme:**
	* **`figure_color`** *(str, default='#0D0D0D')*: The background color of the figure canvas (95% Black). Set to `'white'` for standard paper figures.
    * **`figure_size`** *(tuple, default=(20,20), optional)*: Manually set the figure size (width, height) in inches, e.g., `(10, 10)`. If `None`, defaults to Scanpy's automatic sizing.
	* **`text_color`** *(str, default='auto')*: Color of axes labels and text. `'auto'` automatically selects white for dark backgrounds and black for light backgrounds.

	### Style example: Switching to "White Paper" Mode
	
	```python
	plot_global_embedding_optimized(
	adata, 
	color=['geneA'], 
	figure_color='white',       # White background
	background_color='gainsboro' # Light grey points
	)
    ```
 
	## Other notes
	- The background layer is automatically rasterized.
			
	</details>
2. 


