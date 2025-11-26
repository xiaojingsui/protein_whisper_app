import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

# ----------------------------
# Settings
# ----------------------------
EXCEL_FILE = "Table_S1A.xlsx"
PEP_SHEET = "Peptide-level data"
HEADER_ROW = 3   # Data starts at row 4

st.set_page_config(page_title="Protein Whisper – Peptide Explorer", layout="wide")

@st.cache_data
def load_peptide_data():
    df_pep = pd.read_excel(EXCEL_FILE, sheet_name=PEP_SHEET, header=HEADER_ROW)
    return df_pep

# ----------------------------
# Load data
# ----------------------------
df = load_peptide_data()

st.title("🔬 Protein Whisper – Peptide-level Explorer")
st.write("Interactive exploration of **peptide-level conformation data** from your Table S1A dataset.")

# ----------------------------
# Sidebar search
# ----------------------------
st.sidebar.header("🔍 Search peptides/proteins")

search_gene = st.sidebar.text_input("Gene symbol contains:")
search_peptide = st.sidebar.text_input("Peptide sequence contains:")

df_filtered = df.copy()

if search_gene:
    df_filtered = df_filtered[df_filtered["Gene symbol"].astype(str).str.contains(search_gene, case=False)]

if search_peptide:
    df_filtered = df_filtered[df_filtered["Peptide sequence"].astype(str).str.contains(search_peptide, case=False)]

st.write(f"### Showing {len(df_filtered)} peptides")

st.dataframe(df_filtered, use_container_width=True)

# ----------------------------
# Volcano plot (optional)
# ----------------------------
st.write("---")
st.write("## 🧨 Volcano Plot (Choose condition)")

# Detect available conformation columns
conformation_cols = [c for c in df.columns if "conformation" in c and "AvgLog₂" in c]

if conformation_cols:
    selected_cond = st.selectbox("Select comparison:", conformation_cols)

    # Corresponding p-values
    pval_col = selected_cond.replace("AvgLog₂", "Pval")

    if pval_col in df.columns:
        df_volcano = df[[selected_cond, pval_col, "Gene symbol", "Peptide sequence"]].dropna()

        fig = px.scatter(
            df_volcano,
            x=selected_cond,
            y=-df_volcano[pval_col].apply(lambda p: np.log10(p) * -1),
            hover_data=["Gene symbol", "Peptide sequence"],
            title=f"Volcano plot: {selected_cond}"
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("No p-value column found for this comparison.")
else:
    st.info("No conformation columns found in peptide sheet.")
