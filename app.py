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

    # Extract UniProt IDs from strings like "sp|P02566|UNC-15_CAEEL"
    df["uniprot_id"] = df["Protein ID"].str.split("|").str[1].str.strip()

    # Clean gene symbol
    df["gene"] = df["Gene symbol"].str.replace("CELE_", "", regex=False).str.strip()

    return df


df = load_table()


# ============================================================
# Extract experimental conditions from column names
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
# Robust AlphaFold downloader (supports PDB & CIF)
# ============================================================
@st.cache_data(show_spinner=True)
def download_structure(uniprot):
    """
    Downloads AlphaFold structure using the JSON API.
    Supports CIF or PDB, all model versions (v6, v4, v3...).
    """

    uniprot = uniprot.strip().upper()
    headers = {"User-Agent": "Mozilla/5.0 (Streamlit App)"}

    # Query AlphaFold API for metadata
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    r = requests.get(api_url, headers=headers)

    if r.status_code != 200:
        return None, None, "API request failed"

    try:
        data = r.json()
    except:
        return None, None, "Invalid JSON returned"

    if len(data) == 0:
        return None, None, "No entries found for this UniProt ID"

    # Detect PDB or CIF
    model_url = None
    file_format = None

    # Prefer PDB if available
    for entry in data:
        if "modelUrl" in entry:
            model_url = entry["modelUrl"]
            file_format = "pdb"
            break

    # Otherwise fallback to CIF
    if model_url is None:
        for entry in data:
            if "cifUrl" in entry:
                model_url = entry["cifUrl"]
                file_format = "cif"
                break

    if model_url is None:
        return None, None, "No PDB or CIF model available"

    # Download the model file
    res = requests.get(model_url, headers=headers)
    if res.status_code != 200 or len(res.text) < 500:
        return None, model_url, "Failed to download model file"

    return res.text, model_url, file_format


# ============================================================
# 3D structure viewer
# ============================================================
def render_structure(structure_text, segments, file_format):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(structure_text, file_format)
    view.setStyle({"cartoon": {"color": "white", "opacity": 0.85}})

    # Highlight peptide segments
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

selected_condition = st.sidebar.selectbox("Select condition:", conditions)

fc_cutoff = st.sidebar.number_input(
    "Fold-change cutoff (|AvgLog₂|):", min_value=0.0, max_value=10.0, value=1.0, step=0.1
)

p_cutoff = st.sidebar.number_input(
    "AdjPval cutoff:", min_value=0.0, max_value=1.0, value=0.05, step=0.01
)

# ============================================================
# Main UI
# ============================================================
st.title("🧬 Protein Whisper – Structure Viewer")

if not query:
    st.info("Enter a gene name or UniProt ID above.")
    st.stop()

# Match protein
hits = df[
    (df["uniprot_id"].str.contains(query, case=False, na=False))
    | (df["gene"].str.contains(query, case=False, na=False))
]

if hits.empty:
    st.error(f"No protein found for: {query}")
    st.stop()

protein = hits.iloc[0]
uniprot = protein["uniprot_id"]
gene = protein["gene"]

st.subheader(f"Protein: **{gene}** ({uniprot})")

# Condition columns
avg_col = f"AvgLog₂({selected_condition}).conformation"
pval_col = f"AdjPval({selected_condition}).conformation"

# Filter peptide rows
peps = df[df["uniprot_id"] == uniprot].copy()
peps = peps[(peps[avg_col].abs() >= fc_cutoff) & (peps[pval_col] <= p_cutoff)]

if peps.empty:
    st.warning("No peptides meet the filters.")
else:
    st.success(f"{len(peps)} peptides passed the filters.")

# Convert peptide rows to segments
segments = []
for _, row in peps.iterrows():
    start, end = row["Start position"], row["End position"]
    if pd.isna(start) or pd.isna(end):
        continue

    color = "magenta" if row["Peptide type"] == "full" else "teal"

    segments.append({
        "start": int(start),
        "end": int(end),
        "color": color,
    })

# ============================================================
# Download structure
# ============================================================
structure, model_url, file_format = download_structure(uniprot)

if structure is None:
    st.error(f"No AlphaFold model available for UniProt ID **{uniprot}**.")
    st.warning(f"Reason: {file_format}")
    st.info("Structure view skipped; peptide table is shown below.")
else:
    st.subheader("3D Structure Viewer")
    st.write(f"Using AlphaFold model: {model_url}")
    viewer = render_structure(structure, segments, file_format)
    viewer.show()

# ============================================================
# Peptide table
# ============================================================
st.subheader("Peptides Used for Highlighting")
st.dataframe(
    peps[["Peptide sequence", "Start position", "End position", avg_col, pval_col]]
)
