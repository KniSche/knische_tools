import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import warnings

def plot_global_spatial_rgb(
    adata, 
    basis="global_spatial", 
    uns_key="global_spatial",
    r_genes=None, g_genes=None, b_genes=None, # Default to None for safety
    batch_key=None,
    subset=False,                      
    auto_size_subset=True,             
    margin=0.05,                       
    min_size=1,                        
    max_size=1,                       
    min_alpha=0.1, 
    max_alpha=1.0,
    gamma=0.7,                         
    background_color='grey',
    background_size=None,              
    figure_size=(30, 30),
    figure_color='#0D0D0D',
    v_max_percentile=98,               
    **kwargs
):
    """
    Renders an RGB spatial plot with robust feature-checking. 
    If a gene is missing from the adata, it warns the user and proceeds 
    with a zeroed-out channel for that specific color.
    """

    # --- 1. INITIAL SETUP ---
    save_filename = kwargs.pop('save', None)
    
    if background_size is None:
        area = figure_size[0] * figure_size[1]
        background_size = 0.025 * (area / 100) 

    # --- 2. DATA INTEGRITY & ALIGNMENT ---
    if uns_key not in adata.uns:
        raise KeyError(f"Global reference '{uns_key}' not found in adata.uns.")
    
    full_df = adata.uns[uns_key]
    bg_vals = full_df.values[:, :2]
    
    if basis not in adata.obsm:
        subset_coords = full_df.reindex(adata.obs_names).values[:, :2]
    else:
        subset_coords = adata.obsm[basis][:, :2]

    ## --- 3. ROBUST PER-CHANNEL & PER-GENE SCORING ---
    def get_channel_scores(gene_list, channel_name):
        if gene_list is None or len(gene_list) == 0:
            return np.zeros(adata.n_obs)
            
        found_genes = [g for g in gene_list if g in adata.var_names]
        missing_genes = [g for g in gene_list if g not in adata.var_names]
        
        if missing_genes:
            warnings.warn(f"Genes {missing_genes} not found in adata for {channel_name} channel.")
            
        if not found_genes:
            return np.zeros(adata.n_obs)
            
        # Extract raw expression matrix for these genes
        vals = adata[:, found_genes].X
        if hasattr(vals, "toarray"): vals = vals.toarray()
        
        def scale_vec(v):
            """Core scaling logic: Non-zero 98th percentile Min-Max"""
            if len(v) == 0 or np.all(v == 0): return np.zeros_like(v)
            v_nonzero = v[v > 0]
            if len(v_nonzero) == 0: return np.zeros_like(v)
                
            upper_bound = np.percentile(v_nonzero, v_max_percentile)
            if upper_bound <= 0: upper_bound = v.max()
            
            v_clipped = np.clip(v, 0, upper_bound)
            v_min, v_max = v_clipped.min(), v_clipped.max()
            
            # Return linear 0.0 to 1.0 (Gamma is applied later)
            return (v_clipped - v_min) / (v_max - v_min) if v_max > v_min else np.zeros_like(v)

        def process_matrix(matrix):
            """Applies within-gene scaling, then channel summation"""
            # 1. Scale each gene individually to [0, 1]
            scaled_genes = np.zeros_like(matrix, dtype=float)
            for i in range(matrix.shape[1]):
                scaled_genes[:, i] = scale_vec(matrix[:, i])
                
            # 2. Sum the standardized genes
            summed_scores = scaled_genes.sum(axis=1)
            
            # 3. Scale the final channel sum back to [0, 1] and apply Gamma
            final_channel_norm = scale_vec(summed_scores)
            return np.power(final_channel_norm, gamma)

        # Apply batch correction (if requested)
        if batch_key and batch_key in adata.obs:
            final_vec = np.zeros(adata.n_obs, dtype=float)
            for batch in adata.obs[batch_key].unique():
                mask = (adata.obs[batch_key] == batch).values
                # Slice the matrix for this specific batch and process
                final_vec[mask] = process_matrix(vals[mask, :])
            return final_vec
            
        # Global processing if no batch key is provided
        return process_matrix(vals)

    # Calculate scores with safety checks
    r_norm = get_channel_scores(r_genes, "Red")
    g_norm = get_channel_scores(g_genes, "Green")
    b_norm = get_channel_scores(b_genes, "Blue")

    # --- 4. INTENSITY & ZOOM SCALING ---
    # Driver of size/alpha: The max signal across channels
    # norm_signal = np.max([r_norm, g_norm, b_norm], axis=0)
    norm_signal = r_norm + g_norm + b_norm
    norm_signal[norm_signal > 1] = 1
    
    current_min_s, current_max_s, current_bg_s = min_size, max_size, background_size
    
    if subset and auto_size_subset:
        zoom_factor = np.clip((bg_vals[:, 0].max() - bg_vals[:, 0].min()) / 
                              (subset_coords[:, 0].max() - subset_coords[:, 0].min()), 1.0, 10.0)
        current_min_s *= zoom_factor
        current_max_s *= zoom_factor * 0.5
        current_bg_s *= (zoom_factor * 2)

    # --- 5. RGBA CONSTRUCTION ---
    order = np.argsort(norm_signal)
    
    rgba = np.zeros((len(order), 4))
    rgba[:, 0], rgba[:, 1], rgba[:, 2] = r_norm[order], g_norm[order], b_norm[order]
    
    # Calculate transparency and sizes safely
    rgba[:, 3] = min_alpha + (max_alpha - min_alpha) * norm_signal[order]
    fg_sizes = current_min_s + (current_max_s - current_min_s) * norm_signal[order]
    sorted_coords = subset_coords[order]

    # --- 6. PLOTTING ---
    fig, ax = plt.subplots(figsize=figure_size, facecolor=figure_color)
    ax.set_facecolor(figure_color)
    
    # Background (Atlas)
    ax.scatter(bg_vals[:, 0], bg_vals[:, 1], c=background_color, 
               s=current_bg_s, zorder=1, rasterized=True, edgecolors='none')
    
    # Foreground (Signal)
    ax.scatter(sorted_coords[:, 0], sorted_coords[:, 1], c=rgba, 
               s=fg_sizes, zorder=2, rasterized=True, edgecolors='none', **kwargs)

    if subset:
        x_min, x_max = subset_coords[:, 0].min(), subset_coords[:, 0].max()
        y_min, y_max = subset_coords[:, 1].min(), subset_coords[:, 1].max()
        xr, yr = (x_max - x_min), (y_max - y_min)
        ax.set_xlim(x_min - xr * margin, x_max + xr * margin)
        ax.set_ylim(y_min - yr * margin, y_max + yr * margin)

    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    if save_filename:
        plt.savefig(save_filename, facecolor=figure_color, dpi=300, bbox_inches='tight')
    
    return ax
