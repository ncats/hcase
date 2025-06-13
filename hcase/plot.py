import base64
from io import BytesIO
from PIL import Image, ImageDraw, ImageFont
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFMCS
import seaborn as sns

TRANSPARENT_PNG_BASE64 = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8Xw8AAr8B9FzEBq8AAAAASUVORK5CYII="

def create_custom_linear_bins(freq_matrix, n_bins):
    """
    Create custom bins for frequency:
    - Bin 0: [0, 0]
    - Bin 1: (0, 1]
    - Bin 2: (1, 2)
    - Bin 3: [2, 10)
    - Bin 4: [10, 100)
    - Bin 5: [100, 1000)
    - Bin 6: [1000, 10000)
    - Bin 7: [10000, inf)
    """
    return np.array([0, 1, 2, 10, 100, 1000, 10000, np.inf])

def generate_heatmap_figure(
    agg_combined, plot_type, freq_range, x_centers, y_centers, discrete_colorscale, n_bins,
    special_ref_inchikeys=None, special_unique_ref_inchikeys=None
):
    # Helper to create log-spaced bins, avoids zero by starting from min_bin_val

    # Step 1: Aggregate frequency by bin
    freq_df = agg_combined.groupby(['y_bin', 'x_bin'])["frequency"].sum().unstack(fill_value=0)
    freq_matrix = freq_df.values

    # Build a 2D array of ref_ids matching freq_matrix shape
    refid_df = agg_combined.groupby(['y_bin', 'x_bin'])["ref_id"].first().unstack(fill_value=None)
    refid_matrix = refid_df.values

    # Build a 2D array of base64 images for those ref_ids, but always use blank image for speed
    img_matrix = []
    for row in refid_matrix:
        img_row = []
        for _ in row:
            img_b64 = TRANSPARENT_PNG_BASE64
            img_row.append(img_b64)
        img_matrix.append(img_row)

    # Step 2: Filter bins based on aggregated frequency (after binning)
    if freq_range is not None:
        min_freq, max_freq = freq_range
        mask = (freq_matrix >= min_freq) & (freq_matrix <= max_freq)
    else:
        mask = np.ones_like(freq_matrix, dtype=bool)

    # Step 3: Apply mask to the frequency matrix
    masked_freq_matrix = np.where(mask, freq_matrix, np.nan)

    heatmap_matrix = masked_freq_matrix

    # Step 5: Bin the heatmap values into discrete color bins
    val_min = np.nanmin(masked_freq_matrix)
    val_max = np.nanmax(masked_freq_matrix)

    bins = create_custom_linear_bins(freq_matrix, n_bins)

    # Assign bin 0 to all values that are exactly zero, rest use normal binning
    bin_indices = np.digitize(heatmap_matrix, bins) - 1
    bin_indices = np.where(heatmap_matrix == 0, 0, bin_indices)  # bin 0 for zeros

    # Replace NaNs in heatmap_matrix with bin 0 (lowest bin)
    bin_indices = np.where(np.isnan(heatmap_matrix), 0, bin_indices)

    # Build customdata as [[ [freq, img_b64], ...], ... ]
    customdata = []
    for i in range(len(freq_matrix)):
        row = []
        for j in range(len(freq_matrix[i])):
            freq = freq_matrix[i][j]
            img_b64 = img_matrix[i][j]
            row.append([freq, img_b64])
        customdata.append(row)

    # Step 6: Create Plotly heatmap
    fig = go.Figure(go.Heatmap(
        z=bin_indices.tolist(),
        x=x_centers.tolist(),
        y=y_centers.tolist(),
        colorscale=discrete_colorscale,
        zmin=0,
        zmax=n_bins - 1,
        zauto=False,
        showscale=False,
        customdata=customdata,
        hovertemplate=(
            "dim_1: %{x}<br>"
            "dim_2: %{y}<br>"
            "frequency: %{customdata[0]}<extra></extra>"
        ),
        name=""
    ))

    # Add special markers for ref_inchikeys if provided
    if special_ref_inchikeys is not None:
        # Determine which are unique and which are common
        unique_set = set(special_unique_ref_inchikeys) if special_unique_ref_inchikeys is not None else set()
        all_set = set(special_ref_inchikeys)
        common_set = all_set - unique_set

        # Plot unique (red star)
        if unique_set:
            unique_points = agg_combined[agg_combined['ref_inchikey'].isin(unique_set)]
            if not unique_points.empty:
                unique_x = [x_centers[int(row['x_bin'])] for _, row in unique_points.iterrows()]
                unique_y = [y_centers[int(row['y_bin'])] for _, row in unique_points.iterrows()]
                fig.add_trace(go.Scatter(
                    x=unique_x,
                    y=unique_y,
                    mode='markers',
                    marker=dict(
                        size=24,
                        color='red',
                        symbol='star',
                        opacity=1.0,
                        line=dict(width=3, color='black')
                    ),
                    name='Special Unique InChIKeys',
                    showlegend=True
                ))

        # Plot common (green star)
        if common_set:
            common_points = agg_combined[agg_combined['ref_inchikey'].isin(common_set)]
            if not common_points.empty:
                common_x = [x_centers[int(row['x_bin'])] for _, row in common_points.iterrows()]
                common_y = [y_centers[int(row['y_bin'])] for _, row in common_points.iterrows()]
                fig.add_trace(go.Scatter(
                    x=common_x,
                    y=common_y,
                    mode='markers',
                    marker=dict(
                        size=24,
                        color='green',
                        symbol='star',
                        opacity=1.0,
                        line=dict(width=3, color='black')
                    ),
                    name='Special Common InChIKeys',
                    showlegend=True
                ))

    fig.update_layout(
        title={
            'text': f"Heatmap of {plot_type}",
            'x': 0.5,
            'xanchor': 'center',
        },
        xaxis_title="dim_1",
        yaxis_title="dim_2",
        yaxis_autorange='reversed',
        xaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks='outside'
        ),
        yaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks='outside'
        )
    )

    # Step 7: Return bin edges and ranges
    freq_bins = bins
    min_freq = val_min
    max_freq = val_max

    return fig, freq_bins, min_freq, max_freq


def rgb_tuple_to_str(rgb_tuple):
    r, g, b = [int(255 * x) for x in rgb_tuple]
    return f'rgb({r},{g},{b})'

def add_label_inside_image(img, label, label_height=30):
    width, height = img.size
    draw = ImageDraw.Draw(img)

    try:
        font = ImageFont.truetype("arial.ttf", 14)
    except IOError:
        font = ImageFont.load_default()

    bbox = draw.textbbox((0, 0), label, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    text_x = (width - text_width) // 2
    text_y = height - label_height + (label_height - text_height) // 2  # inside bottom

    draw.text((text_x, text_y), label, fill='black', font=font)
    return img

def get_mcs_match(parent_mol, child_mol, min_atoms_ratio=0.5):
    mcs_result = rdFMCS.FindMCS([parent_mol, child_mol])
    if mcs_result.canceled or not mcs_result.smartsString:
        return None

    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
    parent_match = parent_mol.GetSubstructMatch(mcs_mol)
    child_match = child_mol.GetSubstructMatch(mcs_mol)

    if len(child_match) < min_atoms_ratio * child_mol.GetNumAtoms():
        print(f"MCS too small: {len(child_match)} < {min_atoms_ratio * child_mol.GetNumAtoms()}", flush=True)
        return None

    return parent_match, child_match


def generate_highlighted_overlap_image(parent_smiles, child_smiles, show_atom_mapping = False, min_atoms_ratio = .5):
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    child_mol = Chem.MolFromSmiles(child_smiles)

    if not show_atom_mapping:
        for atom in parent_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        for atom in child_mol.GetAtoms():
            atom.SetAtomMapNum(0)

    result = get_mcs_match(parent_mol, child_mol, min_atoms_ratio=min_atoms_ratio)
    if result is None:
        print(f"[MCS TOO SMALL or NOT FOUND] Parent: {parent_smiles} | Child: {child_smiles}", flush=True)
        return generate_image_from_smiles(parent_smiles)

    rdDepictor.Compute2DCoords(parent_mol)
    parent_match, _ = result
    highlight_atoms = list(parent_match)
    highlight_color = (0.0, 1.0, 0.0)
    highlight_dict = {idx: highlight_color for idx in highlight_atoms}

    drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
    drawer.DrawMolecule(parent_mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_dict)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()
    img = Image.open(BytesIO(img_bytes)).convert("RGBA")
    return img


def generate_highlighted_grid_image(
    parent_smiles_list,
    child_smiles_list,
    parent_reaction_ids,
    page,
    images_per_row,
    rows_per_grid,
    image_size,
    show_atom_mapping=False
):
    if len(parent_smiles_list) != len(child_smiles_list):
        raise ValueError("parent_smiles_list and child_smiles_list must have the same length.")

    images_per_grid = images_per_row * rows_per_grid
    start = page * images_per_grid
    end = start + images_per_grid

    parent_subset = parent_smiles_list[start:end]
    child_subset = child_smiles_list[start:end]
    reaction_subset = parent_reaction_ids[start:end] if parent_reaction_ids else [None] * len(parent_subset)

    images = []

    for i, (parent_smiles, child_smiles) in enumerate(zip(parent_subset, child_subset)):
        parent_mol = Chem.MolFromSmiles(parent_smiles)
        child_mol = Chem.MolFromSmiles(child_smiles)

        if not show_atom_mapping:
            for atom in parent_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in child_mol.GetAtoms():
                atom.SetAtomMapNum(0)

        result = get_mcs_match(parent_mol, child_mol, min_atoms_ratio=0.5)
        if result is None:
            print(f"[NO MATCH or MCS TOO SMALL] Parent: {parent_smiles}\nChild:  {child_smiles}", flush=True)
            images.append(generate_image_from_smiles(child_smiles))
            continue

        parent_match, _ = result
        highlight_atoms = list(parent_match)
        highlight_color = (0.0, 1.0, 0.0)
        highlight_dict = {idx: highlight_color for idx in highlight_atoms}

        drawer = rdMolDraw2D.MolDraw2DCairo(image_size, image_size)
        drawer.DrawMolecule(parent_mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_dict)
        drawer.FinishDrawing()

        img_bytes = drawer.GetDrawingText()
        img = Image.open(BytesIO(img_bytes)).convert("RGBA")
        if reaction_subset[i]:
            img = add_label_inside_image(img, str(reaction_subset[i]))

        images.append(img)

    rows = (len(images) + images_per_row - 1) // images_per_row
    grid_img = Image.new('RGBA', (images_per_row * image_size, rows * image_size), (255, 255, 255, 0))
    for i, img in enumerate(images):
        if img:
            x_pos = (i % images_per_row) * image_size
            y_pos = (i // images_per_row) * image_size
            grid_img.paste(img, (x_pos, y_pos))

    return grid_img


def make_normal_grid_image(
    smiles_list,
    page,
    images_per_row,
    rows_per_grid,
    image_size,
):
    images_per_grid = images_per_row * rows_per_grid
    start = page * images_per_grid
    end = start + images_per_grid
    subset = smiles_list[start:end]

    if not subset:
        return TRANSPARENT_PNG_BASE64

    images = [generate_image_from_smiles(smi) for smi in subset if smi]
    rows = (len(images) + images_per_row - 1) // images_per_row

    grid_img = Image.new('RGBA', (images_per_row * image_size, rows * image_size), (255, 255, 255, 0))
    for i, img in enumerate(images):
        if img:
            x_pos = (i % images_per_row) * image_size
            y_pos = (i // images_per_row) * image_size
            grid_img.paste(img, (x_pos, y_pos))

    return grid_img



def generate_image_from_smiles(smiles, show_atom_mapping=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES for image generation: {smiles}")
        return None
    if not show_atom_mapping:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()
    img = Image.open(BytesIO(img_bytes)).convert("RGBA")
    return img

def image_to_base64(img):
    buffer = BytesIO()
    img.save(buffer, format='PNG')
    return f"data:image/png;base64,{base64.b64encode(buffer.getvalue()).decode()}"

def find_nearest_bin(center_array, val):
    idx = np.abs(center_array - val).argmin()
    return idx

def plot_interactive_hs(
    embedded, target_ids, hc_orders, ref_ids, parent_target_ids, parent_ref_ids, closest_indices_all, parent_reaction_ids, target_inchikeys,
    ref_inchikeys, special_ref_inchikeys=None, special_unique_ref_inchikeys=None, start_port=8050, max_mols_per_point=500, images_per_row = 5, rows_per_grid = 4, image_size = 200 ):
    idx = 0
    data_by_hc_order = {}

    for i, coords in enumerate(embedded):
        hc_order_val = hc_orders[i]

        rows = []
        for j in range(coords.shape[0]):
            dim_1, dim_2 = map(float, coords[j])
            ref_idx = closest_indices_all[idx]
            rows.append({
                "target_id": target_ids[idx],
                "parent_target_id": parent_target_ids[idx],
                "ref_id": ref_ids[ref_idx],
                "parent_ref_id": parent_ref_ids[ref_idx],
                "dim_1": dim_1,
                "dim_2": dim_2,
                "hc_order": hc_order_val,
                "target_inchikey": target_inchikeys[idx],
                "ref_inchikey": ref_inchikeys[ref_idx],
                "parent_reaction_id": parent_reaction_ids[idx],
            })
            idx += 1
        data_by_hc_order.setdefault(hc_order_val, []).extend(rows)

    for i, (hc_order_val, rows) in enumerate(data_by_hc_order.items()):
        df = pd.DataFrame(rows)

        # Compute frequency as the number of targets per (dim_1, dim_2, ref_id)
        freq = (
            df.groupby(["dim_1", "dim_2", "ref_id"])
            .size()
            .reset_index(name="frequency")
        )

        agg_targets = (
            df.groupby(["dim_1", "dim_2", "ref_id"])
            .apply(lambda group: pd.DataFrame({
                "target_id_list": [list(group["target_id"])],
                "parent_target_id_list": [list(group["parent_target_id"])],
                "parent_reaction_id_list": [list(group["parent_reaction_id"])]
            }))
            .reset_index()
            .drop(columns=["level_3"])  # adjust if needed
        )

        df_ref_inchikeys = df[['ref_id', 'ref_inchikey']].drop_duplicates()

        agg_ref = (
            df.groupby(["dim_1", "dim_2", "ref_id"])
            .agg({
                "parent_ref_id": "first"
            })
            .reset_index()
        )
        agg_ref = agg_ref.merge(df_ref_inchikeys, on='ref_id', how='left')

        agg_combined = agg_ref.merge(agg_targets, on=["dim_1", "dim_2", "ref_id"], how="outer")

        agg_combined = agg_combined.merge(
            df[["dim_1", "dim_2", "ref_id", "parent_reaction_id"]],
            on=["dim_1", "dim_2", "ref_id"],
            how="left"
        )

        agg_combined = agg_combined.merge(freq, on=["dim_1", "dim_2", "ref_id"], how="left")

        target_inchikey_group = df.groupby(["dim_1", "dim_2"])["target_inchikey"].unique().reset_index()
        target_inchikey_group.rename(columns={"target_inchikey": "target_inchikey_list"}, inplace=True)
        agg_combined = agg_combined.merge(target_inchikey_group, on=["dim_1", "dim_2"], how="left")

        # Remove duplicate rows caused by merges
        agg_combined = agg_combined.drop_duplicates(subset=["dim_1", "dim_2", "ref_id"])

        # Save selected fields to CSV with hc_order in filename
        save_fields = ["dim_1", "dim_2", "ref_id", "parent_ref_id", "ref_inchikey", "parent_reaction_id", "frequency"]
        csv_filename = f"agg_combined_hc_order_{hc_order_val}.csv"
        agg_combined[save_fields].to_csv(csv_filename, index=False)
        print(f"Saved {csv_filename} with selected fields.")

        # Clip target lists to max_mols_per_point
        agg_combined["target_id_list"] = agg_combined["target_id_list"].apply(lambda lst: lst[:max_mols_per_point] if isinstance(lst, list) else [])
        agg_combined["parent_target_id_list"] = agg_combined["parent_target_id_list"].apply(lambda lst: lst[:max_mols_per_point] if isinstance(lst, list) else [])
        agg_combined["parent_reaction_id_list"] = agg_combined["parent_reaction_id_list"].apply(lambda lst: lst[:max_mols_per_point] if isinstance(lst, list) else [])

        bin_size = 1
        max_dim = 2 ** hc_order_val

        x_min = 0
        x_max = min(agg_combined.dim_1.max(), max_dim)
        y_min = 0
        y_max = min(agg_combined.dim_2.max(), max_dim)

        x_bins = np.arange(x_min, x_max + bin_size, bin_size)
        y_bins = np.arange(y_min, y_max + bin_size, bin_size)

        agg_combined['x_bin'] = np.digitize(agg_combined.dim_1, x_bins) - 1
        agg_combined['x_bin'] = agg_combined['x_bin'].clip(0, len(x_bins) - 2) 
        agg_combined['y_bin'] = np.digitize(agg_combined.dim_2, y_bins) - 1
        agg_combined['y_bin'] = agg_combined['y_bin'].clip(0, len(y_bins) - 2)

        x_centers = x_bins[:-1] + bin_size / 2
        y_centers = y_bins[:-1] + bin_size / 2

        canonical_smiles_to_coord = {}

        # Build bin_data_map with keys = (y_bin, x_bin, ref_id)
        bin_data_map = {}
        for _, row in agg_combined.iterrows():
            smiles = row.get('ref_id')
            if pd.notnull(smiles):
                canonical_smiles_to_coord[smiles] = (row['dim_1'], row['dim_2'])
            key = (row['y_bin'], row['x_bin'], row['ref_id'])
            if key not in bin_data_map:
                bin_data_map[key] = {
                    'ref_id': row['ref_id'],
                    'target_id_list': row['target_id_list'],
                    'parent_ref_id': row['parent_ref_id'],
                    'parent_target_id_list': row['parent_target_id_list'],
                    'parent_reaction_id_list': row['parent_reaction_id_list'],
                    'parent_reaction_id': row['parent_reaction_id']
                }

        n_bins = 10
        # We'll need freq_matrix to create bins, so generate it here
        freq_df = agg_combined.groupby(['y_bin', 'x_bin'])["frequency"].sum().unstack(fill_value=0)
        freq_matrix = freq_df.values
        bins = create_custom_linear_bins(freq_matrix, n_bins)
        n_bins = len(bins) - 1
        colors = sns.color_palette("cubehelix", n_colors=n_bins - 1)
        colors = colors[::-1]  # Invert color scheme: light = low, dark = high
        # First color is white for bin 0 (frequency == 0)
        discrete_colorscale = [(0.0, "rgb(255,255,255)")]
        discrete_colorscale += [
            ((i + 1) / (n_bins - 1), rgb_tuple_to_str(color))
            for i, color in enumerate(colors)
        ]
        # Only generate the frequency plot
        fig, freq_bins, min_raw_freq, max_raw_freq = generate_heatmap_figure(
            agg_combined, "frequency", None,
            x_centers, y_centers, discrete_colorscale, n_bins,
            special_ref_inchikeys=special_ref_inchikeys,
            special_unique_ref_inchikeys=special_unique_ref_inchikeys
        )

        figs = {"frequency": fig}
        legend_data = {
            "frequency": {
                "colors": ["rgb(255,255,255)"] + [rgb_tuple_to_str(c) for c in colors],
                "bins": freq_bins
            }
        }

        app = Dash(__name__)

        def to_human_notation(val):
            if val == 0:
                return "0"
            elif val < 1:
                return f"{val:.0f}"
            elif val < 1_000:
                return f"{int(val)}"
            elif val < 1_000_000:
                return f"{int(val/1_000)}K"
            elif val < 1_000_000_000:
                return f"{int(val/1_000_000)}M"
            else:
                return f"{int(val/1_000_000_000)}B"

        # Show marks for all bins using just the bin edge values
        slider_marks = {i: to_human_notation(legend_data['frequency']['bins'][i]) for i in range(n_bins)}

        # Reverse marks for the grid slider so 0 is at the top
        def make_reversed_marks(max_page):
            return {i: str(max_page - i) for i in range(max_page + 1)}

        app.layout = html.Div([
            dcc.Store(id='hover-coords', data={'x': 0, 'y': 0}),
            dcc.Store(id='mouse-coords', data={'x': 0, 'y': 0}),
            html.Img(
                id='hover-img',
                src=TRANSPARENT_PNG_BASE64,
                style={
                    'position': 'fixed',
                    'top': '0px',
                    'left': '0px',
                    'display': 'none',
                    'zIndex': 1000,
                    'border': '2px solid #333',
                    'backgroundColor': 'white',
                    'height': '120px'
                }
            ),
            # dcc.Store removed, now using clickData
            html.Div([
                html.H2(
                    f"HS Coordinates (hc_order = {hc_order_val})",
                    style={"textAlign": "center", "marginBottom": "10px"}
                ),
                
                html.Div([
                    dcc.RadioItems(
                        id='mode-toggle',
                        options=[
                            {'label': 'Use Ref/Target Scaffolds', 'value': 'normal'},
                            {'label': 'Use Parent Ref/Target Product SMILES', 'value': 'parent'}
                        ],
                        value='parent',
                        labelStyle={'display': 'inline-block', 'marginRight': '15px'},
                    )
                ], style={"textAlign": "center", "marginBottom": "15px"}),

            # Plot type dropdown removed
            ]),

            html.Div([
                html.Label("Search Reference (SMILES or InChIKey):", style={"fontWeight": "bold"}),
                dcc.Input(
                    id='search-input',
                    type='text',
                    debounce=True,
                    placeholder="Enter SMILES or InChIKey",
                    style={'width': '400px'}
                ),
                html.Div(id='search-status', style={"marginTop": "5px", "color": "red"})
            ], style={"textAlign": "center", "marginBottom": "20px"}),
            html.Div([
                html.Label("Search Target (SMILES or InChIKey):", style={"fontWeight": "bold"}),
                dcc.Input(
                    id='target-search-input',
                    type='text',
                    debounce=True,
                    placeholder="Enter SMILES or InChIKey",
                    style={'width': '400px'}
                ),
                html.Div(id='target-search-status', style={"marginTop": "5px", "color": "red"})
            ], style={"textAlign": "center", "marginBottom": "20px"}),



            html.Div([
                html.Label("Frequency Range:", style={"marginBottom": "5px"}),
                dcc.RangeSlider(
                    id='freq-range-slider',
                    min=0,
                    max=n_bins-1,
                    step=1,
                    value=[0, n_bins-1],
                    marks=slider_marks,
                    allowCross=False
                )
            ], style={
                "width": "80%",
                "margin": "0 auto 10px auto",
                "padding": "0",
                "backgroundColor": "transparent",
                "border": "none"
            }),


            html.Div([
                dcc.Graph(
                    id='heatmap',
                    figure=figs['frequency'],
                    style={
                        'height': '700px',
                        'width': '80%'
                    }
                ),
                html.Div(id='color-legend', style={
                    'marginLeft': '5px',
                    'alignSelf': 'center',
                })
            ], style={
                'display': 'flex',
                'flexDirection': 'row',
                'justifyContent': 'flex-center',
                'alignItems': 'center',
            }),

            html.Div([
                html.Div("Reference Space Structure", style={"fontWeight": "bold", "fontSize": 18, "marginBottom": "5px"}),
                html.Img(id='ref-img', src="", style={
                    'height': '150px',
                    'display': 'none',
                    'border': '1px solid #ccc',
                    'backgroundColor': 'white',
                    'margin': '0 auto 10px auto',
                    'maxWidth': '200px',
                    'objectFit': 'contain'
                }),
            ], style={'textAlign': 'center'}),

            html.Div([
                html.Div("Corresponding Target Space Structures", style={
                    "fontWeight": "bold",
                    "fontSize": 18,
                    "marginBottom": "5px",
                    "textAlign": "center"
                }),

                html.Div([
                    html.Img(id='grid-img', src="", style={
                        'height': '600px',
                        'border': '1px solid #ccc',
                        'margin': '10px',
                        'maxWidth': '100%',
                        'objectFit': 'contain'
                    }),

                    html.Div([
                        html.Label("Scroll Pages", style={"fontWeight": "bold", "marginBottom": "10px"}),
                        dcc.Slider(
                            id='grid-page-slider',
                            min=0,
                            max=0,  # will be set dynamically
                            step=1,
                            value=0,
                            marks={0: "0"},
                            vertical=True,
                            updatemode='drag',
                            tooltip={"always_visible": True, "placement": "right"},
                        )
                    ], style={
                        'width': '60px',
                        'height': '600px',
                        'marginLeft': '10px',
                        'display': 'flex',
                        'flexDirection': 'column',
                        'justifyContent': 'center'
                    })
                ], style={
                    'display': 'flex',
                    'flexDirection': 'row',
                    'justifyContent': 'center',
                    'alignItems': 'center',
                    'marginBottom': '30px'
                })
            ])

        ])

        # Clientside callback to update mouse-coords with actual mouse position over the graph
        app.clientside_callback(
            """
            function(n_intervals) {
                // Only attach the event listener once
                if (!window._dash_mousemove_attached) {
                    var graphDiv = document.querySelector('#heatmap');
                    if (graphDiv) {
                        graphDiv.addEventListener('mousemove', function(e) {
                            var coords = {'x': e.clientX, 'y': e.clientY};
                            window.dash_clientside.callback_context.setStoreData('mouse-coords', coords);
                        });
                        window._dash_mousemove_attached = true;
                    }
                }
                return window._dash_last_mouse_coords || {'x': 0, 'y': 0};
            }
            """,
            Output('mouse-coords', 'data'),
            Input('heatmap', 'n_clicks'),  # dummy input to trigger on load
            prevent_initial_call=True
        )

        @app.callback(
            Output('hover-img', 'src'),
            Output('hover-img', 'style'),
            Input('heatmap', 'hoverData')
        )
        def update_hover_img(hoverData):
            if not hoverData or "points" not in hoverData:
                return TRANSPARENT_PNG_BASE64, {
                    'position': 'fixed',
                    'top': '0px',
                    'left': '0px',
                    'display': 'none',
                    'zIndex': 1000,
                    'border': '2px solid #333',
                    'backgroundColor': 'white',
                    'height': '120px'
                }
            # Get data coordinates
            point = hoverData["points"][0]
            x = point["x"]
            y = point["y"]

            # Find the nearest bin and ref_id
            x_bin = find_nearest_bin(x_centers, x)
            y_bin = find_nearest_bin(y_centers, y)
            possible_keys = [key for key in bin_data_map.keys() if key[0] == y_bin and key[1] == x_bin]
            ref_id = None
            if possible_keys:
                ref_id = bin_data_map[possible_keys[0]]['ref_id']

            # Generate the image for the ref_id, or fallback to not rendering
            if ref_id:
                img = generate_image_from_smiles(ref_id)
                img_src = image_to_base64(img) if img else TRANSPARENT_PNG_BASE64

                # Assume graph is 80% of window width and 700px height, and centered horizontally
                import dash
                window_width = 1200  # fallback default
                graph_width = int(window_width * 0.8)
                graph_height = 700
                graph_left = int((window_width - graph_width) / 2)
                graph_top = 300  # estimate vertical offset from top of page

                # Get axis ranges (assume linear, not log)
                x_min = min(x_centers)
                x_max = max(x_centers)
                y_min = min(y_centers)
                y_max = max(y_centers)

                # Map data coordinates to pixel coordinates
                px = int((x - x_min) / (x_max - x_min) * graph_width) + graph_left
                py = int((y - y_min) / (y_max - y_min) * graph_height) + graph_top

                style = {
                    'position': 'fixed',
                    'top': f"{py + 10}px",
                    'left': f"{px + 10}px",
                    'display': 'block',
                    'zIndex': 1000,
                    'border': '2px solid #333',
                    'backgroundColor': 'white',
                    'height': '120px'
                }
                return img_src, style
            else:
                # No ref_id, do not render image
                return TRANSPARENT_PNG_BASE64, {
                    'display': 'none'
                }

        @app.callback(
            Output('heatmap', 'figure'),
            Output('search-status', 'children'),
            Output('target-search-status', 'children'),
            Input('search-input', 'value'),
            Input('target-search-input', 'value'),
            Input('freq-range-slider', 'value'),
        )
        def update_figure_with_search(search_input, target_search_input, freq_range):
            # freq_range is now [bin_start_idx, bin_end_idx]
            bins = legend_data["frequency"]["bins"]
            # Use all indices as valid
            bin_start = min(range(n_bins), key=lambda x: abs(x - freq_range[0]))
            bin_end = min(range(n_bins), key=lambda x: abs(x - freq_range[1]))
            # Convert to frequency values
            min_freq = bins[bin_start]
            max_freq = bins[bin_end+1] if bin_end+1 < len(bins) else bins[-1]
            freq_range_actual = [min_freq, max_freq]

            fig, freq_bins, min_freq, max_freq = generate_heatmap_figure(
                agg_combined, "frequency", freq_range_actual,
                x_centers, y_centers, discrete_colorscale, n_bins,
                special_ref_inchikeys=special_ref_inchikeys,
                special_unique_ref_inchikeys=special_unique_ref_inchikeys
            )

            status = ""
            if search_input:
                search_input = search_input.strip()
                # Try SMILES first
                mol = Chem.MolFromSmiles(search_input)
                if mol is not None:
                    try:
                        inchikey = Chem.inchi.MolToInchiKey(mol)
                        match = agg_combined[agg_combined["ref_inchikey"] == inchikey]
                        if not match.empty:
                            x, y = match.iloc[0]['dim_1'], match.iloc[0]['dim_2']
                            fig.add_trace(go.Scatter(
                                x=[x], y=[y], mode='markers+text',
                                marker=dict(size=15, color='red', symbol='x')
                            ))
                            status = f"Match found for SMILES (InChIKey: {inchikey})"
                        else:
                            status = f"No match found for SMILES (InChIKey: {inchikey})"
                    except Exception as e:
                        status = f"Error parsing SMILES: {e}"
                else:
                    # Check if input looks like an InChIKey (27 chars, 2 hyphens, uppercase)
                    ik = search_input.strip().upper()
                    if len(ik) == 27 and ik.count('-') == 2 and all(c.isalnum() or c == '-' for c in ik):
                        match = agg_combined[agg_combined["ref_inchikey"] == ik]
                        if not match.empty:
                            x, y = match.iloc[0]['dim_1'], match.iloc[0]['dim_2']
                            fig.add_trace(go.Scatter(
                                x=[x], y=[y], mode='markers+text',
                                marker=dict(size=15, color='blue', symbol='circle'),
                                text=["Match (InChIKey)"], textposition="top center"
                            ))
                            status = f"Match found for InChIKey: {ik}"
                        else:
                            status = f"No match found for InChIKey: {ik}"
                    else:
                        status = "Input is not a valid SMILES or InChIKey."

            target_status = ""
            if target_search_input:
                target_search_input = target_search_input.strip()
                # Try SMILES first
                mol = Chem.MolFromSmiles(target_search_input)
                if mol is not None:
                    try:
                        inchikey = Chem.inchi.MolToInchiKey(mol)
                        match = agg_combined[agg_combined["target_inchikey_list"].apply(lambda lst: inchikey in lst if isinstance(lst, (list, np.ndarray)) else False)]
                        if not match.empty:
                            xs = match['dim_1'].tolist()
                            ys = match['dim_2'].tolist()
                            fig.add_trace(go.Scatter(
                                x=xs, y=ys, mode='markers+text',
                                marker=dict(size=15, color='green', symbol='diamond'),
                                text=["Target"]*len(xs), textposition="top center"
                            ))
                            target_status = f"Match found for Target SMILES (InChIKey: {inchikey})"
                        else:
                            target_status = f"No match found for Target SMILES (InChIKey: {inchikey})"
                    except Exception as e:
                        target_status = f"Error parsing SMILES: {e}"
                else:
                    # Check if input looks like an InChIKey (27 chars, 2 hyphens, uppercase)
                    ik = target_search_input.strip().upper()
                    if len(ik) == 27 and ik.count('-') == 2 and all(c.isalnum() or c == '-' for c in ik):
                        match = agg_combined[agg_combined["target_inchikey_list"].apply(lambda lst: ik in lst if isinstance(lst, (list, np.ndarray)) else False)]
                        if not match.empty:
                            xs = match['dim_1'].tolist()
                            ys = match['dim_2'].tolist()
                            fig.add_trace(go.Scatter(
                                x=xs, y=ys, mode='markers+text',
                                marker=dict(size=15, color='orange', symbol='star'),
                                text=["Target (InChIKey)"]*len(xs), textposition="top center"
                            ))
                            target_status = f"Match found for Target InChIKey: {ik}"
                        else:
                            target_status = f"No match found for Target InChIKey: {ik}"
                    else:
                        target_status = "Input is not a valid SMILES or InChIKey."

            return fig, status, target_status




        @app.callback(
            Output('grid-img', 'src'),
            Input('heatmap', 'clickData'),
            Input('mode-toggle', 'value'),
            Input('grid-page-slider', 'value')
        )
        def update_grid_image(clickData, mode, page):
            if not clickData or "points" not in clickData:
                return TRANSPARENT_PNG_BASE64

            point = clickData["points"][0]
            x, y = point["x"], point["y"]

            x_bin = find_nearest_bin(x_centers, x)
            y_bin = find_nearest_bin(y_centers, y)

            possible_keys = [key for key in bin_data_map.keys() if key[0] == y_bin and key[1] == x_bin]
            if not possible_keys:
                return TRANSPARENT_PNG_BASE64
            selected_key = possible_keys[0]

            # Use page directly (no inversion)
            if mode == 'normal':
                smiles_list = bin_data_map[selected_key]['target_id_list']
                grid_img = make_normal_grid_image(
                    smiles_list,
                    page=page,
                    images_per_row=images_per_row,
                    rows_per_grid=rows_per_grid,
                    image_size=image_size
                )
            else:
                grid_img = generate_highlighted_grid_image(
                    parent_smiles_list=bin_data_map[selected_key]['parent_target_id_list'],
                    child_smiles_list=bin_data_map[selected_key]['target_id_list'],
                    parent_reaction_ids=bin_data_map[selected_key]['parent_reaction_id_list'],
                    page=page,
                    images_per_row=images_per_row,
                    rows_per_grid=rows_per_grid,
                    image_size=image_size
                )

            return image_to_base64(grid_img)

        @app.callback(
            [Output('ref-img', 'src'), Output('ref-img', 'style')],
            [Input('heatmap', 'clickData'),
             Input('mode-toggle', 'value')]
        )
        def update_ref_image_and_position(clickData, mode):
            if not clickData or "points" not in clickData:
                return TRANSPARENT_PNG_BASE64, {
                    'height': '150px',
                    'border': '1px solid #ccc',
                    'backgroundColor': 'white',
                    'margin': '20px auto 10px auto',
                    'display': 'block',  # keep display block so space is reserved
                    'opacity': '0'  # optionally make it fully transparent but still occupies space
                }
            point = clickData["points"][0]
            x, y = point["x"], point["y"]
            x_bin = find_nearest_bin(x_centers, x)
            y_bin = find_nearest_bin(y_centers, y)

            possible_keys = [key for key in bin_data_map.keys() if key[0] == y_bin and key[1] == x_bin]
            if not possible_keys:
                return TRANSPARENT_PNG_BASE64, {
                    'height': '150px',
                    'border': '1px solid #ccc',
                    'backgroundColor': 'white',
                    'margin': '20px auto 10px auto',
                    'display': 'block',  # keep display block so space is reserved
                    'opacity': '0'  # optionally make it fully transparent but still occupies space
                }

            selected_key = possible_keys[0]
            ref_id = bin_data_map[selected_key]['ref_id']
            if not ref_id:
                return TRANSPARENT_PNG_BASE64, {
                    'height': '150px',
                    'border': '1px solid #ccc',
                    'backgroundColor': 'white',
                    'margin': '20px auto 10px auto',
                    'display': 'block',  # keep display block so space is reserved
                    'opacity': '0'  # optionally make it fully transparent but still occupies space
                }
            img = generate_image_from_smiles(ref_id)

            src = image_to_base64(img) if img else ""
            style = {
                'height': '150px',
                'border': '1px solid #ccc',
                'backgroundColor': 'white',
                'margin': '20px auto 10px auto',
                'display': 'block',
            }
            return src, style

        
        @app.callback(
            Output('color-legend', 'children'),
            Input('freq-range-slider', 'value'),
        )
        def update_color_legend(freq_range):
            colors = legend_data["frequency"]["colors"]
            bins = legend_data["frequency"]["bins"]

            def bin_label(i):
                if i == 0:
                    return "0"
                elif i == len(colors) - 1:
                    return f"[{to_human_notation(bins[i])}, ∞)"
                else:
                    return f"[{to_human_notation(bins[i])}, {to_human_notation(bins[i+1])})"

            return html.Div([
                html.H4("Color Legend"),
                html.Ul([
                    html.Li(style={'listStyle': 'none', 'marginBottom': '5px'}, children=[
                        html.Span(style={
                            'display': 'inline-block',
                            'width': '20px',
                            'height': '20px',
                            'backgroundColor': colors[i],
                            'marginRight': '10px',
                            'border': '1px solid #000'
                        }),
                        bin_label(i)
                    ]) for i in reversed(range(len(colors)))
                ])
            ], style={'textAlign': 'left'})

        @app.callback(
            Output('grid-page-slider', 'max'),
            Output('grid-page-slider', 'value'),
            Output('grid-page-slider', 'marks'),
            Input('heatmap', 'clickData')
        )
        def update_slider_max(clickData):
            if not clickData or "points" not in clickData:
                return 0, 0, {0: "0"}

            point = clickData["points"][0]
            x, y = point["x"], point["y"]
            x_bin = find_nearest_bin(x_centers, x)
            y_bin = find_nearest_bin(y_centers, y)

            possible_keys = [key for key in bin_data_map if key[0] == y_bin and key[1] == x_bin]
            if not possible_keys:
                return 0, 0, {0: "0"}

            selected_key = possible_keys[0]
            total = len(bin_data_map[selected_key]['target_id_list'])
            per_page = images_per_row * rows_per_grid
            max_page = (total - 1) // per_page

            # Use increasing marks so 0 is at the top and matches the tooltip
            marks = {i: str(i) for i in range(max_page + 1)}

            return max_page, 0, marks

        port = start_port + i * 10
        print(f"Dash app running on port {port} for hc_order={hc_order_val}")
        app.run_server(port=port, host='0.0.0.0')
