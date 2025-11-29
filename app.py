import streamlit as st
import pandas as pd
import py3Dmol
import requests
import re

# ============================================================
# Load peptide table
# ============================================================
@st.cache_data
def load_table():
    df = pd.read_excel("Table_S1A.xlsx", header=3)
    df["uniprot_id"] = df["Protein ID"].str.split("|").str[1]
    df["gene"] = df["Gene symbol"].str.replace("CELE_", "", regex=False)
    return df

df = load_table()

# ============================================================
# Extract available experimental conditions
# ============================================================
def extract_conditions(df):
    conds = []
    for col in df.columns:
        m = re.match(r"AvgLog₂\((.+)\)\.conformation", col)
        if m:
            conds.append(m.group(1))
    return sorted(conds)

conditions = extract_conditions(df)

# ============================================================
# Download AlphaFold PDB on the fly
# ============================================================
@st.cache_data(show_spinner=True)
def download_pdb(uniprot):
    """
    Fetch PDB structure from AlphaFold DB dynamically.
    Works on Streamlit Cloud.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        raise ValueError(f"Could not download PDB for {uniprot}")
    return r.text

# ============================================================
# 3D Viewer
# ============================================================
def render_structure(pdb_data, segments):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "white", "opacity": 0.85}})

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

query = st.sidebar.text_input("Gene name / UniProt ID:", "")

selected_condition = st.sidebar.selectbox(
    "Condition to visualize:",
    conditions
)

fc_cut = st.sidebar.number_input(
    "Fold-change cutoff (|AvgLog₂|):",
    min_value=0.0, max_value=10.0, value=1.0, step=0.1
)

p_cut = st.sidebar.number_input(
    "AdjPval cutoff:",
    min_value=0.0, max_value=1.0, value=0.05, step=0.01
)

# ============================================================
# Main view
# ============================================================
st.title("Protein Whisper – Structure Viewer")

if not query:
    st.info("Search a protein to begin.")
    st.stop()

# ============================================================
# Find protein match
# ============================================================
hits = df[
    (df["uniprot_id"].str.contains(query, case=False, na=False)) |
    (df["gene"].str.contains(query, case=False, na=False))
]

if hits.empty:
    st.error("No matching protein found.")
    st.stop()

protein = hits.iloc[0]
uniprot = protein["uniprot_id"]
gene = protein["gene"]

st.subheader(f"Protein: **{gene}** ({uniprot})")

# ============================================================
# Condition column names
# ============================================================
avg_col = f"AvgLog₂({selected_condition}).conformation"
pval_col = f"AdjPval({selected_condition}).conformation"

if avg_col not in df.columns or pval_col not in df.columns:
    st.error("Selected condition not found in dataset.")
    st.stop()

# ============================================================
# Extract significant peptides for this protein
# ============================================================
sub = df[df["uniprot_id"] == uniprot].copy()

sub = sub[
    (sub[avg_col].abs() >= fc_cut) &
    (sub[pval_col] <= p_cut)
]

if sub.empty:
    st.warning("No peptides meet filter criteria.")
    st.stop()

# Create segment list
segments = []
for _, row in sub.iterrows():
    if pd.isna(row["Start position"]) or pd.isna(row["End position"]):
        continue

    color = "magenta" if row["Peptide type"] == "full" else "teal"

    segments.append({
        "start": int(row["Start position"]),
        "end": int(row["End position"]),
        "color": color
    })

# ============================================================
# Download AlphaFold structure
# ============================================================
with st.spinner("Downloading AlphaFold structure..."):
    pdb_data = download_pdb(uniprot)

# ============================================================
# Show 3D structure
# ============================================================
st.subheader(f"Structure Highlighted by Condition: {selected_condition}")

viewer = render_structure(pdb_data, segments)
viewer.show()

# ============================================================
# Show peptide table
# ============================================================
st.subheader("Highlighted Peptides")
st.dataframe(
    sub[["Peptide sequence", "Start position", "End position", avg_col, pval_col]]
)
