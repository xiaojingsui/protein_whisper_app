import streamlit as st
import pandas as pd
import py3Dmol
import requests
import re
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ============================================================
# Load peptide table
# ============================================================
@st.cache_data
def load_table():
    df = pd.read_excel("Table_S1A.xlsx", header=3)

    # Extract UniProt ID like P02566 from "sp|P02566|UNC-54_CAEEL"
    df["uniprot_id"] = (
        df["Protein ID"]
        .astype(str)
        .str.split("|")
        .str[1]
        .str.strip()
    )

    # Clean gene symbol
    df["gene"] = (
        df["Gene symbol"]
        .astype(str)
        .str.replace("CELE_", "", regex=False)
        .str.strip()
    )

    return df


df = load_table()


# ============================================================
# Extract experimental conditions (dynamic)
# ============================================================
def extract_conditions(df):
    conds = []
    for col in df.columns:
        m = re.match(r"AvgLog₂\((.+)\)\.conformation", col)
        if m:
            conds.append(m.group(1))
    return sorted(set(conds))


conditions = extract_conditions(df)


# ============================================================
# AlphaFold downloader (supports CIF & PDB, v6/v4/v3)
# ============================================================
@st.cache_data(show_spinner=True)
def download_structure(uniprot):
    uniprot = uniprot.strip().upper()
    headers = {"User-Agent": "Mozilla/5.0 (Streamlit App)"}

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    r = requests.get(api_url, headers=headers)

    if r.status_code != 200:
        return None, None, None, "API request failed"

    try:
        data = r.json()
    except Exception:
        return None, None, None, "Invalid JSON from API"

    if len(data) == 0:
        return None, None, None, "AlphaFold API returned empty list"

    model_url = None
    file_format = None

    # Prefer PDB if present
    for entry in data:
        if "modelUrl" in entry:
            model_url = entry["modelUrl"]
            file_format = "pdb"
            break

    # Fallback to CIF
    if model_url is None:
        for entry in data:
            if "cifUrl" in entry:
                model_url = entry["cifUrl"]
                file_format = "cif"
                break

    if model_url is None:
        return None, None, None, "No PDB or CIF URL in API response"

    file_data = requests.get(model_url, headers=headers)
    if file_data.status_code != 200 or len(file_data.text) < 500:
        return None, model_url, file_format, "Failed to download structure file"

    return file_data.text, model_url, file_format, None


# ============================================================
# 3D Viewer
# ============================================================
def render_structure(structure_text, segments, file_format, plddt_coloring):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(structure_text, file_format)

    # Base coloring: either pLDDT (B-factor) or neutral white
    if plddt_coloring:
        # AlphaFold-style B-factor gradient
        view.setStyle(
            {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "linear",
                        "min": 50,
                        "max": 100,
                        "colors": ["orange", "yellow", "cyan", "blue"],
                    }
                }
            }
        )
    else:
        view.setStyle({"cartoon": {"color": "white", "opacity": 0.85}})

    # Highlight peptide segments
    for seg in segments:
        view.addStyle(
            {"resi": list(range(seg["start"], seg["end"] + 1))},
            {"cartoon": {"color": seg["color"], "opacity": 1.0}},
        )

    view.zoomTo()
    return view


# ============================================================
# Sidebar UI
# ============================================================
st.sidebar.header("Protein Search")

query = st.sidebar.text_input("Search gene or UniProt ID:", "")

selected_condition = st.sidebar.selectbox("Select condition:", conditions)

fc_cutoff = st.sidebar.number_input(
    "Fold-change cutoff (|AvgLog₂|):",
    min_value=0.0,
    max_value=10.0,
    value=1.0,
    step=0.1,
)

p_cutoff = st.sidebar.number_input(
    "AdjPval cutoff:",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01,
)

color_mode = st.sidebar.selectbox(
    "Peptide color mode:",
    ["Peptide type (magenta/teal)", "Fold-change heatmap"],
)

plddt_coloring = st.sidebar.checkbox(
    "Color backbone by pLDDT (AlphaFold confidence)", value=False
)


# ============================================================
# Main UI
# ============================================================
st.title("🧬 Protein Whisper – Structure Viewer")

if not query:
    st.info("Enter a gene name or UniProt ID to begin.")
    st.stop()

# Match protein
hits = df[
    df["uniprot_id"].str.contains(query, case=False, na=False)
    | df["gene"].str.contains(query, case=False, na=False)
]

if hits.empty:
    st.error(f"No protein found matching '{query}'.")
    st.stop()

protein = hits.iloc[0]
uniprot = protein["uniprot_id"]
gene = protein["gene"]

st.subheader(f"Protein: **{gene}** ({uniprot})")


# ============================================================
# Condition-specific filtering
# ============================================================
avg_col = f"AvgLog₂({selected_condition}).conformation"
pval_col = f"AdjPval({selected_condition}).conformation"

if avg_col not in df.columns or pval_col not in df.columns:
    st.error("Condition does not exist in the dataset.")
    st.stop()

# All peptides for this protein (for volcano)
prot_all = df[df["uniprot_id"] == uniprot].copy()

# Significant peptides for highlighting
peps = prot_all[
    (prot_all[avg_col].abs() >= fc_cutoff) & (prot_all[pval_col] <= p_cutoff)
].copy()

if peps.empty:
    st.warning("No peptides meet the filter criteria.")
else:
    st.success(f"{len(peps)} peptides passed the filters.")


# ============================================================
# Build segments with colors
# ============================================================
segments = []

if not peps.empty:
    if color_mode == "Fold-change heatmap":
        # Diverging color map centered at 0 (blue → white → red)
        fc_vals = peps[avg_col].astype(float)
        fc_max = float(fc_vals.abs().max())
        if fc_max == 0:
            fc_max = 1.0
        norm = mcolors.TwoSlopeNorm(vmin=-fc_max, vcenter=0.0, vmax=fc_max)
        cmap = plt.get_cmap("coolwarm")

    for _, row in peps.iterrows():
        start = row["Start position"]
        end = row["End position"]
        if pd.isna(start) or pd.isna(end):
            continue

        if color_mode == "Peptide type (magenta/teal)":
            color = "magenta" if row["Peptide type"] == "full" else "teal"
        else:
            fc = float(row[avg_col])
            rgba = cmap(norm(fc))
            color = mcolors.to_hex(rgba)

        segments.append(
            {
                "start": int(start),
                "end": int(end),
                "color": color,
            }
        )


# ============================================================
# Download and Display Structure
# ============================================================
structure_text, model_url, file_format, error_msg = download_structure(uniprot)

st.subheader("3D Structure Viewer")

if structure_text is None:
    st.error(f"No AlphaFold structure available for **{uniprot}**.")
    if error_msg:
        st.warning(f"Reason: {error_msg}")
    st.info("Showing peptide table and volcano plot only.")
else:
    st.write(f"Using AlphaFold model: {model_url}")

    viewer = render_structure(structure_text, segments, file_format, plddt_coloring)
    html = viewer._make_html()
    components.html(html, height=650, scrolling=True)


# ============================================================
# Volcano Plot for Selected Condition
# ============================================================
st.subheader("Volcano Plot for Selected Condition")

x_all = prot_all[avg_col].astype(float)
y_all = -np.log10(prot_all[pval_col].astype(float) + 1e-300)

fig, ax = plt.subplots(figsize=(6, 4))

# All peptides (light grey)
ax.scatter(x_all, y_all, s=8, alpha=0.3, color="lightgrey", label="All peptides")

# Highlight significant peptides (those on the structure)
if not peps.empty:
    x_sig = peps[avg_col].astype(float)
    y_sig = -np.log10(peps[pval_col].astype(float) + 1e-300)
    ax.scatter(x_sig, y_sig, s=25, alpha=0.9, color="black", label="Highlighted")

# Threshold lines
ax.axvline(fc_cutoff, color="red", linestyle="--", linewidth=1)
ax.axvline(-fc_cutoff, color="red", linestyle="--", linewidth=1)
ax.axhline(-np.log10(p_cutoff + 1e-300), color="blue", linestyle="--", linewidth=1)

ax.set_xlabel(f"AvgLog₂({selected_condition})")
ax.set_ylabel("-log₁₀(AdjPval)")
ax.legend(loc="best", fontsize=8)
ax.grid(alpha=0.2)

st.pyplot(fig)


# ============================================================
# Peptide Table
# ============================================================
st.subheader("Peptides Used for Highlighting")
st.dataframe(
    peps[
        [
            "Peptide sequence",
            "Start position",
            "End position",
            avg_col,
            pval_col,
        ]
    ]
)
