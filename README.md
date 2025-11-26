# KniSche tools

> A selection of custom functions for bioinformatics analysis centred on single cell and spatial next generation sequencing.

## Installation

You can install this list of tools as a package directly from GitHub using pip:

```bash
pip install git+[https://github.com/knische/knische_tools.git](https://github.com/knische/knische_tools.git)
```

## List of Functions
1. `plot_global_spatial` 
	This tool is plot wrapper function around `sc.pl.embedding`, allowing users to plot subset embeddings
	from `adata.obsm["embedding"]` onto the embedding made from the entire set.
	Primary use is solving the problem of visualizing a small subset of cells 
	(e.g., a specific tissue section or cell type) while maintaining the context of the full dataset's embedding
	without requiring the user to manually manage massive coordinate files.
	<details>
	<summary></summary>
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
	from knische_spatial import plot_global_spatial

	# Simple plot
	plot_global_spatial(
		adata, 
		color=['cell_type', 'gene_of_interest']
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
	* **`background_color`** *(str, default='lightgrey')*: Color of the global reference points.
	* **`background_size`** *(int, default=5)*: Size of the global reference points.
	* **`subset_size`** *(int, optional)*: Size of the subset points. If `None`, defaults to Scanpy's automatic sizing.
	* **`remove_overlapping_background`** *(bool, default=True)*: If `True`, removes global points that sit directly underneath subset points. Recommended when using transparency/alpha.
	* **`**kwargs`**: Additional arguments passed to `sc.pl.embedding` (e.g., `cmap`, `vmax`, `frameon`, `alpha`).

	## Other notes
	- The background layer is automatically rasterized to keep file sizes small when exporting to PDF/SVG.
			
	</details>
2. 


