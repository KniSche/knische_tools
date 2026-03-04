import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np

def plot_global_spatial(
    adata, 
    color=None, 
    basis="global_spatial",         
    uns_key="global_spatial",          
    subset=False,                      # If True, zooms the view to the foreground data
    auto_size_subset=True,             # Scales spot sizes based on zoom level
    margin=0.1,                        # Fractional padding around the subset zoom
    background_color='grey',           # Color of the global reference points
    background_size=None,              # Auto-calculated if None
    size=None,                         # Foreground spot size
    figure_size=(30.34, 20),           # (Width, Height) in inches
    figure_color='#0D0D0D',            # Canvas background color
    text_color='auto',                 # Switches between white/black based on figure_color
    **kwargs
):
    """
    Overlays subset AnnData results onto a global coordinate reference.
    
    This function handles the alignment of subset data to a master coordinate 
    dataframe (stored in adata.uns), renders a rasterized background for 
    performance, and optionally 'zooms' into the subset area while 
    dynamically scaling spot sizes to maintain visual density.

    Parameters
    ----------
    adata : AnnData
        The subset object containing the cells to be colored/analyzed.
    color : str or list, optional
        Keys in adata.obs or adata.var_names to color by.
    basis : str
        The obsm key for coordinates. If missing, it will be aligned from the master.
    uns_key : str
        The key in adata.uns where the master coordinate DataFrame is stored.
    """

    # --- 1. INITIAL SETUP & SIZE HEURISTICS ---
    save_filename = kwargs.pop('save', None)
    
    # Calculate a baseline background size if not provided (proportional to figure area)
    if background_size is None:
        area = figure_size[0] * figure_size[1]
        background_size = 0.025 * (area / 100) 
        
    if size is None:
        size = 2.0

    # --- 2. DATA ALIGNMENT ---
    if uns_key not in adata.uns:
        raise KeyError(f"Global reference '{uns_key}' not found in adata.uns.")
    
    full_df = adata.uns[uns_key]
    if not isinstance(full_df, pd.DataFrame):
        raise TypeError(f"Global reference must be a Pandas DataFrame.")

    # Align coordinates: If basis is missing from obsm, map it from the master df
    obsm_keys = adata.obsm.keys()
    if basis not in obsm_keys and f"X_{basis}" not in obsm_keys:
        adata.obsm[basis] = full_df.reindex(adata.obs_names).values[:, :2]
    
    # Extract coordinate arrays for scaling calculations
    bg_vals = full_df.values[:, :2]
    actual_basis = basis if basis in adata.obsm else f"X_{basis}"
    fg_vals = adata.obsm[actual_basis][:, :2]

    # --- 3. SMART ZOOM & SCALING LOGIC ---
    current_bg_s = background_size
    current_fg_s = size

    if subset and auto_size_subset:
        # Determine the linear zoom factor based on X-axis range ratio
        g_x_range = bg_vals[:, 0].max() - bg_vals[:, 0].min()
        s_x_range = fg_vals[:, 0].max() - fg_vals[:, 0].min()
        
        if s_x_range > 0:
            zoom_factor = np.clip(g_x_range / s_x_range, 1.0, 10.0)
            
            # Scale sizes: foreground grows to fill space, background grows less to remain subtle
            current_fg_s *= (zoom_factor * 0.5)
            current_bg_s *= (zoom_factor * 2)

    # --- 4. COLOR THEME MANAGEMENT ---
    if text_color == 'auto':
        # Detect if figure background is dark to set appropriate text contrast
        dark_shades = ['black', '#000000', '#0d0d0d']
        is_dark = figure_color.lower() in dark_shades or figure_color.startswith(('#0', '#1', '#2'))
        calc_text_color = 'white' if is_dark else 'black'
    else:
        calc_text_color = text_color

    # Apply aesthetic context
    rc_params = {
        'figure.facecolor': figure_color,
        'axes.facecolor': figure_color,
        'text.color': calc_text_color,
        'axes.labelcolor': calc_text_color,
        'xtick.color': calc_text_color,
        'ytick.color': calc_text_color,
        'axes.edgecolor': calc_text_color,
        'figure.figsize': figure_size
    }

    # --- 5. PLOTTING EXECUTION ---
    with plt.rc_context(rc_params):
        axes_list = sc.pl.embedding(
            adata, 
            basis=basis, 
            color=color,
            size=current_fg_s, 
            show=False,
            zorder=2, 
            **kwargs
        )

    # Ensure axes_list is iterable for multi-panel support
    if not isinstance(axes_list, list):
        axes_list = [axes_list]
    
    for ax in axes_list:
        ax.set_aspect('equal', adjustable='box')

        if hasattr(ax, 'scatter'):
            # Render Background (Rasterized to handle millions of points without lag)
            ax.scatter(
                bg_vals[:, 0], bg_vals[:, 1], 
                c=background_color, 
                s=current_bg_s, 
                zorder=1,
                rasterized=True,
                edgecolors='none'
            )
            
            # Apply coordinate bounding box if subset zoom is requested
            if subset:
                x_min, x_max = fg_vals[:, 0].min(), fg_vals[:, 0].max()
                y_min, y_max = fg_vals[:, 1].min(), fg_vals[:, 1].max()
                
                xr, yr = (x_max - x_min), (y_max - y_min)
                
                ax.set_xlim(x_min - xr * margin, x_max + xr * margin)
                ax.set_ylim(y_min - yr * margin, y_max + yr * margin)

    # Final Save Handling
    if save_filename:
        plt.savefig(save_filename, facecolor=figure_color, bbox_inches='tight', dpi=300)
            
    return axes_list
