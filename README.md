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
	without requiring the user to manually manage massive coordinate files. Later versions will include the actual 
	histology as "background" instead of all points plotted as a raster image. 
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
	
	### Argument documentation
	**Parameters:**
	
	* **`adata`** *(AnnData)*: The subset object. Must contain global coords in `adata.uns[uns_key]`.
	* **`color`** *(str or list)*: Keys for annotations of observations/cells or variables/genes (e.g., `'cell_type'` or `['geneA', 'geneB']`).
	* **`basis`** *(str, default='global_spatial')*: The key to store/access the aligned coordinates in `adata.obsm`.
	* **`uns_key`** *(str, default='global_spatial')*: The key in `adata.uns` where the full reference DataFrame is stored.
	* **`subset`** * (
	* **`size`** *(int, default=10, optional)*: Size of the subset points. If `None`, defaults to Scanpy's automatic sizing.
    	* **`**kwargs`**: Additional arguments passed to `sc.pl.embedding` (e.g., `cmap`, `vmax`, `frameon`, `alpha`).
	
	**Background underlay options:**
	* **`background_size`** *(int, default=0.0002)*: Size of the global reference points.
    	* **`background_color`** *(str, default='lightgrey')*: Color of the **global reference points**. *Tip: If using a dark figure_color, try a darker grey (e.g., '#333333') for a subtle effect.*
	
	**Style & Theme:**
	* **`figure_color`** *(str, default='#0D0D0D')*: The background color of the figure canvas (95% Black). Set to `'white'` for standard paper figures.
    	* **`figure_size`** *(tuple, default=(20,20), optional)*: Manually set the figure size (width, height) in inches, e.g., `(10, 10)`. If `None`, defaults to Scanpy's automatic sizing.
	* **`text_color`** *(str, default='auto')*: Color of axes labels and text. `'auto'` automatically selects white for dark backgrounds and black for light backgrounds.

	**Save image:**
    	* **`save`** *(str, default=None)*: File name to save to, e.g., `'plot.png'`. unlike sc.pl.embedding it will not default to `figures/*`, and the full path needs to be specified.
 
	### Style example: Switching to "White Paper" Mode
	
	```python
	plot_global_spatial(
	adata, 
	color=['geneA'], 
	figure_color='white',       # White background
	background_color='gainsboro' # Light grey points
	)
	plt.show()
    	```
 
	## Other notes
	- The background layer is automatically rasterized.
			
	</details>
2. `plot_global_spatial_rgb`
	This function is scatter plot function independent of `sc.pl.embedding` to plot RGB channels for different 
	gene groups on the global object. In principal it receives similar arguments and plots the `adata.obsm["global_spatial"]`
	embedding. Later versions will include the actual histology as "background" instead of all points plotted as a raster image.	
	<details>Usage and Documentation</summary>
	
	### Argument documentation
	**Parameters:**
	* **`adata`** *(AnnData)*: The subset object containing gene expression. Must contain global coordinates in `adata.uns[uns_key]`.
	* **`r_genes`, `g_genes`, `b_genes`** *(list, default=None)*: Lists of genes to be mapped to the Red, Green, and Blue channels. Features not found in `adata.var_names` are gracefully skipped with a warning.
	* **`basis`** *(str, default='global_spatial')*: The key used to store/access the aligned coordinates in `adata.obsm`.
	* **`uns_key`** *(str, default='global_spatial')*: The key in `adata.uns` where the full reference DataFrame is stored.
	* **`batch_key`** *(str, default=None)*: Column in `adata.obs` used for technical batch correction. Normalization is performed independently within each batch to account for varying signal intensities.
	* **`**kwargs`**: Additional arguments passed to the underlying Matplotlib `ax.scatter` call (e.g., `linewidths`, `marker`, `edgecolors`).

	**RGB Logic & Normalization:**
	* **`gamma`** *(float, default=0.7)*: Power-law brightness correction ($I_{out} = I_{in}^{\gamma}$). Values less than **1.0** brighten the mid-tones, which is essential for making signals visible against dark backgrounds.
	* **`v_max_percentile`** *(int, default=98)*: Percentile used to define the maximum intensity. This protects the visualization from being washed out by high-expressing outliers.
	* **`min_alpha` / `max_alpha`** *(float, default=0.1 / 1.0)*: The range of transparency. Alpha is scaled based on the maximum normalized signal across all three channels.

	**Zoom & Scaling:**
	* **`subset`** *(bool, default=False)*: If `True`, the plot viewport is automatically restricted to the bounding box of the foreground cells.
	* **`auto_size_subset`** *(bool, default=True)*: Dynamically increases both foreground and background spot sizes when `subset=True` to maintain visual density as you "zoom in."
	* **`margin`** *(float, default=0.05)*: Fractional padding (e.g., **5%**) added around the coordinates when the viewport is restricted.
	* **`min_size` / `max_size`** *(int, default=2 / 20)*: The range of spot sizes for foreground cells, scaled by normalized expression intensity.

	**Background Underlay Options:**
	* **`background_size`** *(float, default=None)*: Size of the global reference points. If `None`, size is auto-calculated based on the figure area.
	* **`background_color`** *(str, default='grey')*: Color of the **global reference points**. 
	    > **Tip:** If using a dark `figure_color`, try a darker grey (e.g., `'#333333'`) for a more subtle "atlas" effect.

	**Style & Theme:**
	* **`figure_color`** *(str, default='#0D0D0D')*: The background color of the figure canvas. Set to `'white'` for standard publication-ready figures.
	* **`figure_size`** *(tuple, default=(30, 30))*: Dimensions of the figure in inches. 
	* **`text_color`** *(str, default='auto')*: Color for axes labels and text. `'auto'` automatically toggles between white and black based on the brightness of `figure_color`.

	**Save Image:**
	* **`save`** *(str, default=None)*: File name or full path to save the figure (e.g., `'plots/spatial_rgb.png'`). Unlike standard Scanpy functions, this does not default to a `figures/` directory; the path must be explicit.
	
	### Example:
	
	```python
	plot_global_spatial_rgb(
	adata, 
	r_genes=['geneA']
	g_genes=['geneB', 'geneC'], 
	b_genes=None, 
	figure_color='black',       # Black background (default)
	background_color='grey' # Grey points (default)
	)
	plt.show()
    	```
	</details>
	
## Tutorials & Vignettes
* [**Basic usage & advanced spatial plots**](docs/vignettes/plotGlobalSpatial_tutorial.ipynb): A deep dive into multi-channel visualization, gamma correction, and batch normalization.
	
	
