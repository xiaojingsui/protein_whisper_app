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
# Load peptide-level table (Sheet S1A)
# ============================================================
@st.cache_data
def load_table():
    df = pd.read_excel("Table_S1A.xlsx", header=3)

    df["uniprot_id"] = (
        df["Protein ID"].astype(str).str.split("|").str[1].str.strip()
    )
    df["gene"] = (
        df["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip()
    )
    return df


df = load_table()


# ============================================================
# Load abundance table (Sheet S1B)
# ============================================================
@st.cache_data
def load_abundance_table():
    df2 = pd.read_excel("Table_S1A.xlsx", sheet_name="Protein-level data", header=3)

    df2["uniprot_id"] = (
        df2["Protein ID"].astype(str).str.split("|").str[1].str.strip()
    )
    df2["gene"] = (
        df2["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip()
    )
    return df2


abun_df = load_abundance_table()


# ============================================================
# Extract conformation-only conditions
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
# AlphaFold Downloader
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
    except:
        return None, None, None, "Invalid JSON response"

    if len(data) == 0:
        return None, None, None, "AlphaFold returned empty list"

    model_url, file_format = None, None

    # Prefer PDB
    for entry in data:
        if "modelUrl" in entry:
            model_url = entry["modelUrl"]
            file_format = "pdb"
            break

    # Else CIF
    if model_url is None:
        for entry in data:
            if "cifUrl" in entry:
                model_url = entry["cifUrl"]
                file_format = "cif"
                break

    if model_url is None:
        return None, None, None, "No modelUrl or cifUrl available"

    res = requests.get(model_url, headers=headers)
    if res.status_code != 200 or len(res.text) < 500:
        return None, model_url, file_format, "Failed to download structure file"

    return res.text, model_url, file_format, None


# ============================================================
# 3D Viewer
# ============================================================
def render_structure(structure_text, segments, file_format, plddt_coloring):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(structure_text, file_format)

    if plddt_coloring:
        view.setStyle(
            {"cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "linear",
                    "min": 50,
                    "max": 100,
                    "colors": ["orange", "yellow", "cyan", "blue"]
                }
            }}
        )
    else:
        view.setStyle({"cartoon": {"color": "white", "opacity": 0.85}})

    # Highlight peptides
    for seg in segments:
        view.addStyle(
            {"resi": list(range(seg["start"], seg["end"] + 1))},
            {"cartoon": {"color": seg["color"], "opacity": 1.0}}
        )

    view.zoomTo()
    return view


# ============================================================
# Sidebar UI
# ============================================================
st.sidebar.header("Protein Search")

query = st.sidebar.text_input("Search gene or UniProt ID:", "")

selected_condition = st.sidebar.selectbox("Select conformation condition:", conditions)

fc_cutoff = st.sidebar.number_input(
    "Fold-change cutoff (|AvgLog₂|, conformation only):",
    min_value=0.0,
    max_value=10.0,
    value=1.0,
    step=0.1
)

p_cutoff = st.sidebar.number_input(
    "AdjPval cutoff (conformation only):",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01
)

color_mode = st.sidebar.selectbox(
    "Peptide color mode:",
    ["Peptide type (magenta/teal)", "Fold-change heatmap"]
)

plddt_coloring = st.sidebar.checkbox(
    "Color backbone by pLDDT (AlphaFold confidence)", value=False
)


# ============================================================
# MAIN UI
# ============================================================
st.title("🧬 Protein Whisper – Structure Viewer")

if not query:
    st.info("Enter a gene or UniProt ID.")
    st.stop()

hits = df[
    df["uniprot_id"].str.contains(query, case=False, na=False)
    | df["gene"].str.contains(query, case=False, na=False)
]

if hits.empty:
    st.error(f"No protein found for '{query}'.")
    st.stop()

protein = hits.iloc[0]
uniprot = protein["uniprot_id"]
gene = protein["gene"]

st.subheader(f"Protein: **{gene}** ({uniprot})")


# ============================================================
# Conformation Filtering
# ============================================================
avg_col = f"AvgLog₂({selected_condition}).conformation"
pval_col = f"AdjPval({selected_condition}).conformation"

prot_all = df[df["uniprot_id"] == uniprot]

peps = prot_all[
    (prot_all[avg_col].abs() >= fc_cutoff)
    & (prot_all[pval_col] <= p_cutoff)
]

if peps.empty:
    st.warning("No peptides meet conformation FC/Pval filtering.")
else:
    st.success(f"{len(peps)} peptides highlighted on the structure.")


# ============================================================
# Build peptide highlight segments
# ============================================================
segments = []

if not peps.empty:
    if color_mode == "Fold-change heatmap":
        fc_vals = peps[avg_col].astype(float)
        maxfc = max(1.0, float(fc_vals.abs().max()))
        norm = mcolors.TwoSlopeNorm(vmin=-maxfc, vcenter=0, vmax=maxfc)
        cmap = plt.get_cmap("coolwarm")

    for _, row in peps.iterrows():
        start, end = row["Start position"], row["End position"]
        if pd.isna(start) or pd.isna(end):
            continue

        if color_mode == "Peptide type (magenta/teal)":
            color = "magenta" if row["Peptide type"] == "full" else "teal"
        else:
            fc = float(row[avg_col])
            color = mcolors.to_hex(cmap(norm(fc)))

        segments.append({
            "start": int(start),
            "end": int(end),
            "color": color
        })


# ============================================================
# Structure Viewer
# ============================================================
st.subheader("3D Structure Viewer")

structure_text, model_url, file_format, error_msg = download_structure(uniprot)

if structure_text is None:
    st.error("No AlphaFold model available.")
    st.warning(error_msg)
else:
    st.write(f"Using AlphaFold model: {model_url}")
    viewer = render_structure(structure_text, segments, file_format, plddt_coloring)
    components.html(viewer._make_html(), height=650, scrolling=True)


# ============================================================
# Volcano Plot (Conformation Only)
# ============================================================
st.subheader("Volcano Plot (Conformation)")

x_all = prot_all[avg_col].astype(float)
y_all = -np.log10(prot_all[pval_col].astype(float) + 1e-300)

fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(x_all, y_all, color="lightgrey", alpha=0.5, s=12)

if not peps.empty:
    ax.scatter(
        peps[avg_col].astype(float),
        -np.log10(peps[pval_col].astype(float) + 1e-300),
        color="black",
        s=30
    )

ax.axvline(fc_cutoff, color="red", linestyle="--")
ax.axvline(-fc_cutoff, color="red", linestyle="--")
ax.axhline(-np.log10(p_cutoff), color="blue", linestyle="--")

ax.set_xlabel(f"AvgLog₂({selected_condition})")
ax.set_ylabel("-log10(AdjPval)")
ax.grid(alpha=0.25)

st.pyplot(fig)


# ============================================================
# Abundance Plots (Soluble / Pellet / Total)
# ============================================================
st.subheader("Protein Abundance (Soluble / Pellet / Total)")

prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]

# Metric suffixes
abundance_metrics = {
    "Soluble": "soluble",
    "Pellet": "pellet",
    "Total": "total",
}

fig3, axs = plt.subplots(1, 3, figsize=(9, 3))
fig3.tight_layout(pad=3.0)

cond = selected_condition

for i, (label, suffix) in enumerate(abundance_metrics.items()):
    ax = axs[i]

    avg_col_abun = f"AvgLog₂({cond}).{suffix}"
    pval_col_abun = f"AdjPval({cond}).{suffix}"

    if avg_col_abun not in prot_abun.columns:
        ax.set_title(label)
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
        continue

    yval = float(prot_abun[avg_col_abun].values[0])
    pval = prot_abun[pval_col_abun].values[0] if pval_col_abun in prot_abun.columns else None

    ax.bar([label], [yval], color="#4e79a7")
    ax.set_title(label, fontsize=10)

    if i == 0:
        ax.set_ylabel("AvgLog₂")

    if pval is not None:
        ax.text(
            0,
            yval + (0.05 if yval >= 0 else -0.05),
            f"AdjP = {pval:.3g}",
            ha="center",
            va="bottom" if yval >= 0 else "top",
            fontsize=9
        )

    ax.grid(alpha=0.2)
    ax.set_ylim(
        min(0, yval) - abs(yval) * 0.3 - 0.1,
        max(0, yval) + abs(yval) * 0.3 + 0.1
    )

st.pyplot(fig3)


# ============================================================
# Peptide Table
# ============================================================
st.subheader("Peptides Passing Conformation Filters")
st.dataframe(
    peps[[
        "Peptide sequence",
        "Start position",
        "End position",
        avg_col,
        pval_col
    ]]
)
