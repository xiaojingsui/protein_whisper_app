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

    # Extract UniProt ID from Protein ID like "sp|P02566|UNC-15_CAEEL"
    df["uniprot_id"] = df["Protein ID"].str.split("|").str[1].str.strip()

    # Clean gene symbol
    df["gene"] = df["Gene symbol"].str.replace("CELE_", "", regex=False).str.strip()

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
    return sorted(set(conds))

conditions = extract_conditions(df)


# ============================================================
# Download AlphaFold structure from API (robust to v6/v4/v3/etc)
# ============================================================
@st.cache_data(show_spinner=True)
def download_pdb(uniprot):
    uniprot = uniprot.strip().upper()

    # 1. Query AlphaFold API for metadata
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    r = requests.get(api_url)

    if r.status_code != 200:
        return None, None

    try:
        data = r.json()
        model_url = data[0]["modelUrl"]  # e.g., AF-P02566-F1-model_v6.pdb
    except:
        return None, None

    # 2. Fetch the actual PDB file
    pdb_response = requests.get(model_url)
    if pdb_response.status_code == 200 and "NoSuchKey" not in pdb_response.text:
        return pdb_response.text, model_url

    return None, None


# ============================================================
# 3D Viewer
# ============================================================
def render_structure(pdb_data, segments):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(pdb_data, "pdb")
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

selected_condition = st.sidebar.selectbox(
    "Select condition:",
    conditions,
)

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


# ============================================================
# Main UI
# ============================================================
st.title("🧬 Protein Whisper – Structure Viewer")


# Stop if no query
if not query:
    st.info("Enter a gene name or UniProt ID to start.")
    st.stop()


# ============================================================
# Find matching protein entries
# ============================================================
hits = df[
    (df["uniprot_id"].str.contains(query, case=False, na=False))
    | (df["gene"].str.contains(query, case=False, na=False))
]

if hits.empty:
    st.error(f"No protein found matching: {query}")
    st.stop()

protein = hits.iloc[0]
uniprot = protein["uniprot_id"]
gene = protein["gene"]

st.subheader(f"Protein: **{gene}**  ({uniprot})")


# ============================================================
# Condition-specific fields
# ============================================================
avg_col = f"AvgLog₂({selected_condition}).conformation"
pval_col = f"AdjPval({selected_condition}).conformation"

if avg_col not in df.columns or pval_col not in df.columns:
    st.error("Condition does not exist in the dataset.")
    st.stop()


# ============================================================
# Filter peptides for this protein
# ============================================================
peps = df[df["uniprot_id"] == uniprot].copy()

peps = peps[
    (peps[avg_col].abs() >= fc_cutoff)
    & (peps[pval_col] <= p_cutoff)
]

if peps.empty:
    st.warning("No peptides meet the filter criteria.")
else:
    st.success(f"{len(peps)} peptides passed the filters.")


# Build list of peptide segments
segments = []
for _, row in peps.iterrows():
    start = row["Start position"]
    end = row["End position"]

    if pd.isna(start) or pd.isna(end):
        continue

    pep_type = row["Peptide type"]
    color = "magenta" if pep_type == "full" else "teal"

    segments.append({
        "start": int(start),
        "end": int(end),
        "color": color,
    })


# ============================================================
# Download PDB from AlphaFold (API)
# ============================================================
pdb_data, model_url = download_pdb(uniprot)

if pdb_data is None:
    st.warning(f"No AlphaFold model available for UniProt ID **{uniprot}**.")
    st.info("Peptide table is shown below, but structure visualization is skipped.")
else:
    st.write(f"AlphaFold model URL: {model_url}")
    st.subheader("3D Structure Viewer")
    viewer = render_structure(pdb_data, segments)
    viewer.show()


# ============================================================
# Show peptide table
# ============================================================
st.subheader("Peptides Used for Highlighting")
st.dataframe(
    peps[[
        "Peptide sequence",
        "Start position",
        "End position",
        avg_col,
        pval_col,
    ]]
)
