import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np

def plot_global_spatial(
    adata, 
    color=None, 
    basis="global_spatial", 
    uns_key="global_spatial",
    background_color='lightgrey',
    background_size=0.0002,
    size=10,
    fig_size=(20,20),
    figure_color='#0D0D0D',           # Color of the background / canvas
    text_color='auto',                # 'auto' switches based on figure_color
    **kwargs
):
    """
    Plots the subset data overlaying the global dataset stored in adata.uns.
    
    1. Checks if 'basis' exists in adata.obsm. If not, creates it by slicing 
       adata.uns[uns_key] to match adata.obs_names.
    2. Plots the global background (from uns) and the subset foreground (from obsm).
    
    Parameters:
    -----------
    adata : AnnData
        The subset object. Must contain global coords in adata.uns[uns_key].
    color : str or list
        Variables in adata.obs or adata.var to plot.
    basis : str
        The key to use for plotting. Defaults to "global_spatial".
        The function will look for adata.obsm[basis] or adata.obsm['X_' + basis].
    uns_key : str
        The key in adata.uns where the full dataframe is stored.
    """
    
    # --- STEP 1: Data Integrity Check & Alignment ---
    
    # Check if the global master file exists
    if uns_key not in adata.uns:
        raise KeyError(f"'{uns_key}' not found in adata.uns. Please store the full embedding dataframe there first.")
    
    # Get the master dataframe
    full_df = adata.uns[uns_key]
    
    # Ensure full_df is a DataFrame (needed for index matching)
    if not isinstance(full_df, pd.DataFrame):
        raise TypeError(f"adata.uns['{uns_key}'] must be a Pandas DataFrame with an index matching cell barcodes.")

    # Check if the specific embedding slot is missing or empty
    # We check for both 'basis' and 'X_basis' as Scanpy often prepends 'X_'
    obsm_keys = adata.obsm.keys()
    if basis not in obsm_keys and f"X_{basis}" not in obsm_keys:
        print(f"Generating aligned coordinates for adata.obsm['{basis}']...")
        
        # Robust alignment: Reindex ensures we get the exact rows for this subset
        # .reindex fills missing cells with NaN rather than crashing
        subset_coords = full_df.reindex(adata.obs_names)
        
        # Check for missing data
        if subset_coords.isna().any().any():
            print("Warning: Some cells in the subset were not found in the global reference (NaNs created).")
            
        # Store in obsm (Scanpy expects numpy arrays in obsm)
        adata.obsm[basis] = subset_coords.values

    # --- STEP 2: Plotting ---

    # Prepare background coordinates (Full Atlas) for Matplotlib
    # We take the first two columns of the dataframe in uns
    bg_vals = full_df.values[:, :2]

    # Plot the Subset (Foreground) using Scanpy
    # zorder=2 forces these points to sit ON TOP of the background
    # We only apply fig_size if the user requested it AND they didn't pass an existing axis
    # If text_color is 'auto', choose white for dark backgrounds, black for light.
    if text_color == 'auto':
        # Simple heuristic: if figure_color is dark (starts with #0 or #1 or Black)
        if figure_color.lower() in ['black', '#000000', '#0d0d0d'] or figure_color.startswith(('#0', '#1', '#2')):
            calc_text_color = 'white'
        else:
            calc_text_color = 'black'
    else:
        calc_text_color = text_color

    # --- 2. Setup Plotting Context ---
    rc_params = {
        'figure.facecolor': figure_color,
        'axes.facecolor': figure_color,
        'text.color': calc_text_color,
        'axes.labelcolor': calc_text_color,
        'xtick.color': calc_text_color,
        'ytick.color': calc_text_color,
        'axes.edgecolor': calc_text_color,
        # 'grid.color': calc_text_color # Optional: if you use grids
    }

    # Add custom fig_size if provided
    if fig_size is not None and 'ax' not in kwargs:
        rc_params['figure.fig_size'] = fig_size

    # Apply the context just for this plotting command
    with plt.rc_context(rc_params):
        axes_list = sc.pl.embedding(
            adata, 
            basis=basis, 
            color=color,
            size=size,
            show=False,
            zorder=2, 
            **kwargs
        )

    # Handle Single vs Multiple Panels (Standardize to list)
    if not isinstance(axes_list, list):
        axes_list = [axes_list]

    # Iterate through every axis and paint the background
    for ax in axes_list:
        # Verify it's a plot axis (avoids coloring the colorbar axis)
        if hasattr(ax, 'scatter'):
            ax.scatter(
                bg_vals[:, 0], 
                bg_vals[:, 1], 
                c=background_color, 
                s=background_size,
                zorder=1, # zorder=1 forces this BEHIND the subset
                rasterized=True # Optimization for large backgrounds
            )

    return axes_list
