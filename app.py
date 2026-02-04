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
# Page Configuration
# ============================================================
st.set_page_config(
    page_title="Protein Whisper",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ============================================================
# 4. GLOBAL CSS & NAVIGATION STYLING
# ============================================================
st.markdown("""
    <style>
    /* --- GLOBAL FONTS & MAIN CONTAINER --- */
    /* Force Arial on everything in the app */
    html, body, [data-testid="stAppViewContainer"], .stApp, p, h1, h2, h3, h4, h5, h6, span, div {
        font-family: Arial, Helvetica, sans-serif !important;
        background-color: #FBFEFF;
    }
    
    /* Remove standard top padding so the navbar sits at the very top */
    .block-container {
        padding-top: 1rem !important;
        padding-bottom: 5rem !important;
    }

    /* Hide standard Streamlit header elements */
    [data-testid="stHeader"] { display: none !important; }
    #MainMenu { visibility: hidden; }
    footer { visibility: hidden; }

    /* --- CUSTOM NAVBAR (STYLING THE RADIO BUTTON) --- */
    div[role="radiogroup"] {
        position: fixed !important;      
        top: 0 !important;                
        left: 0 !important;               
        width: 100vw !important;          
        z-index: 99999 !important;        
        
        background-color: #FFFFFF;        
        
        display: flex !important;
        justify-content: center !important; 
        padding: 10px 0 !important; 
        
        align-items: center !important; 
        border-bottom: 1px solid #E0E0E0;
    }

    div[role="radiogroup"] label > div:first-child {
        display: none !important;
    }

    div[role="radiogroup"] label {
        margin-right: 0px !important;
        background-color: transparent !important;
        border: none !important;
    }

    div[role="radiogroup"] p {
        font-family: Arial, Helvetica, sans-serif !important;  
        font-size: 18px !important;
        font-weight: 600 !important;
        color: #445550 !important; 
        cursor: pointer;
        padding: 8px 20px;
        border-radius: 20px;
        transition: all 0.3s ease;
        background-color: transparent;
    }

    div[role="radiogroup"] p:hover {
        background-color: #F0FBFC;
        color: #006064 !important;
    }
    
    /* Highlight the selected item (Streamlit usually wraps selected text in a generic span, 
       but we can approximate focus/active states or just rely on the UI feedback) */
    div[role="radiogroup"] label[data-testid="stRadioOption"] {
        background-color: transparent !important;
    }
    </style>
""", unsafe_allow_html=True)

# ============================================================
# Data Loading & Caching
# ============================================================
@st.cache_data
def load_table():
    try:
        df = pd.read_excel("Table_S1A.xlsx", header=3)
        df["uniprot_id"] = df["Protein ID"].astype(str).str.split("|").str[1].str.strip()
        df["gene"] = df["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip()
        return df
    except FileNotFoundError:
        return pd.DataFrame()

@st.cache_data
def load_abundance_table():
    try:
        df2 = pd.read_excel("Table_S1A.xlsx", sheet_name="Protein-level data", header=3)
        df2["uniprot_id"] = df2["Protein ID"].astype(str).str.split("|").str[1].str.strip()
        df2["gene"] = df2["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip()
        return df2
    except:
        return pd.DataFrame()

df = load_table()
abun_df = load_abundance_table()

# Stop if data didn't load
if df.empty:
    st.error("Data file 'Table_S1A.xlsx' not found or empty.")
    st.stop()

# ============================================================
# Helper Functions
# ============================================================
def extract_conditions(df):
    conds = []
    for col in df.columns:
        m = re.match(r"AvgLog₂\((.+)\)\.conformation", col)
        if m:
            conds.append(m.group(1))
    return sorted(set(conds))

def sort_conditions(conds):
    priority = [
        "wt day6/wt day1", "wt day9/wt day1", "Q35/wt aging",
        "wt heat-shock 35°C/wt 20°C", "Q35/Q24", "Q40/Q24",
        "myosin-ts 15°C/wt 15°C", "myosin-ts 25°C/wt 25°C",
        "paramyosin-ts 15°C/wt 15°C", "paramyosin-ts 25°C/wt 25°C"
    ]
    ordered = [c for c in priority if c in conds]
    remaining = sorted(set(conds) - set(ordered))
    return ordered + remaining

conditions_raw = extract_conditions(df)
conditions = sort_conditions(conditions_raw)

@st.cache_data(show_spinner=True)
def download_structure(uniprot):
    uniprot = uniprot.strip().upper()
    headers = {"User-Agent": "Mozilla/5.0 (Streamlit App)"}
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    
    try:
        r = requests.get(api_url, headers=headers)
        if r.status_code != 200: return None, None, None, "API request failed"
        data = r.json()
    except:
        return None, None, None, "Invalid JSON response"

    if len(data) == 0: return None, None, None, "AlphaFold returned empty list"

    model_url, file_format = None, None
    for entry in data:
        if "modelUrl" in entry:
            model_url = entry["modelUrl"]
            file_format = "pdb"
            break
    if model_url is None:
        for entry in data:
            if "cifUrl" in entry:
                model_url = entry["cifUrl"]
                file_format = "cif"
                break

    if model_url is None: return None, None, None, "No modelUrl or cifUrl found"

    res = requests.get(model_url, headers=headers)
    if res.status_code != 200 or len(res.text) < 500:
        return None, model_url, file_format, "Failed to download structure"

    return res.text, model_url, file_format, None

def render_structure(structure_text, segments, file_format, plddt_coloring):
    view = py3Dmol.view(width=900, height=650)
    view.addModel(structure_text, file_format)

    if plddt_coloring:
        view.setStyle(
            {"cartoon": {
                "colorscheme": {
                    "prop": "b",
                    "gradient": "linear",
                    "min": 50, "max": 100,
                    "colors": ["orange", "yellow", "cyan", "blue"]
                }
            }}
        )
    else:
        view.setStyle({"cartoon": {"color": "white", "opacity": 0.85}})

    for seg in segments:
        view.addStyle(
            {"resi": list(range(seg["start"], seg["end"] + 1))},
            {"cartoon": {"color": seg["color"], "opacity": 1.0}}
        )

    view.zoomTo()
    return view

# ============================================================
# NAVIGATION CONTROLLER
# ============================================================
# Use a horizontal radio button, hidden labels, and CSS to make it look like a navbar
page = st.radio(
    "Main Navigation", 
    ["Search", "About", "Guides"], 
    horizontal=True,
    label_visibility="collapsed"
)

# SPACER: Because the navbar is fixed (position: fixed), it sits ON TOP of the content.
# We need to push the content down so it doesn't hide behind the navbar.
st.markdown('<div style="height: 80px;"></div>', unsafe_allow_html=True)

# ============================================================
# PAGE: ABOUT
# ============================================================
if page == "About":
    st.title("About Protein Whisper")
    st.markdown("""
    **Protein Whisper** is an interactive visualization tool designed to explore Limited Proteolysis-Mass Spectrometry (LiP-MS) data. 
    
    It allows researchers to map peptide-level structural alterations directly onto 3D protein structures predicted by AlphaFold.

    ### Key Features
    * **Structure Mapping:** Visualizes peptides that undergo significant structural changes.
    * **Dual-Layer Data:** Integrates conformational data with protein solubility and total abundance.
    * **AlphaFold Integration:** Automatically fetches predicted structures from the EBI AlphaFold database.
    
    ### Data Source
    This tool utilizes data from **Table S1A**, containing LiP-MS peptide intensity data.
    """)

# ============================================================
# PAGE: GUIDES
# ============================================================
elif page == "Guides":
    st.title("User Guide")
    
    st.markdown("### 1. How to Search")
    st.info("Click the **Search** tab above to begin.")
    st.markdown("""
    1.  **Gene/Protein:** Enter a *C. elegans* gene symbol (e.g., `unc-54`) or a UniProt ID in the sidebar.
    2.  **Condition:** Select the stress condition you wish to compare.
    3.  **Filters:** Adjust the Fold-Change and P-value cutoffs.
    """)

    st.markdown("### 2. Visualization")
    st.markdown("""
    * **3D Viewer:** Backbone is white; **Red** = fully-tryptic, **Cyan** = semi-tryptic (structural change).
    * **Volcano Plot:** Shows statistical significance.
    * **Abundance Plot:** Displays Soluble vs. Pellet vs. Total abundance.
    """)

# ============================================================
# PAGE: SEARCH
# ============================================================
elif page == "Search":
    # --- Sidebar for inputs (Only visible on Search page) ---
    st.sidebar.header("Search Parameters")
    query = st.sidebar.text_input("Search gene or UniProt ID:", "")
    selected_condition = st.sidebar.selectbox("Stress condition:", conditions)
    
    with st.sidebar.expander("Filter Settings", expanded=True):
        fc_cutoff = st.number_input(
            "Fold-change cutoff (|AvgLog₂|):", value=1.0, min_value=0.0, max_value=10.0, step=0.1)
        p_cutoff = st.number_input(
            "AdjPval cutoff:", value=0.05, min_value=0.0, max_value=1.0, step=0.01)

    with st.sidebar.expander("Visualization Settings"):
        color_mode = st.selectbox(
            "Peptide color mode:",
            ["Peptide type (red/cyan)", "Fold-change heatmap"]
        )
        plddt_coloring = st.checkbox(
            "Color backbone by AlphaFold pLDDT", value=False
        )

    # --- Main Content ---
    if not query:
        st.markdown("<h2 style='text-align: center; color: #444;'>Structure Viewer</h2>", unsafe_allow_html=True)
        st.markdown(
        """
        <style>
        .stAlert { background-color: #e0faff !important; border: 1px solid #19CFE2 !important; border-radius: 8px; }
        </style>
        """, unsafe_allow_html=True)

        st.info("Type a C. elegans gene name or UniProt ID in the **sidebar** to explore a protein. Try 'unc-54'.")

        st.markdown(
            """
            <div style="text-align: center; margin-top: 50px;">
                <span style="font-size: 50px; color: #19CFE2;">⬅</span>
                <br>
                <span style="font-size: 24px; color: #aaa; font-weight: bold;">Start your search in the sidebar</span>
            </div>
            """, unsafe_allow_html=True)
    else:
        # Search Logic
        hits = df[
            df["uniprot_id"].str.contains(query, case=False, na=False)
            | df["gene"].str.contains(query, case=False, na=False)
        ]

        if hits.empty:
            st.error(f"No protein found for '{query}'.")
        else:
            protein = hits.iloc[0]
            uniprot = protein["uniprot_id"]
            gene = protein["gene"]
            
            st.markdown(f"### Protein: **{gene}** ({uniprot})")
            
            # --- Filtering ---
            avg_col = f"AvgLog₂({selected_condition}).conformation"
            pval_col = f"AdjPval({selected_condition}).conformation"
            prot_all = df[df["uniprot_id"] == uniprot]
            peps = prot_all[
                (prot_all[avg_col].abs() >= fc_cutoff) & (prot_all[pval_col] <= p_cutoff)
            ]

            if peps.empty:
                st.warning("No peptides meet filters.")
            else:
                st.success(f"{len(peps)} significant peptides found.")

            # --- Segments ---
            segments = []
            full_color = "#E1341E"
            half_color = "#1ECBE1"

            if not peps.empty:
                if color_mode == "Fold-change heatmap":
                    fc_vals = peps[avg_col].astype(float)
                    maxfc = max(1.0, float(fc_vals.abs().max()))
                    cmap = plt.get_cmap("coolwarm")
                    norm = mcolors.TwoSlopeNorm(vmin=-maxfc, vcenter=0, vmax=maxfc)

                for _, row in peps.iterrows():
                    start, end = row["Start position"], row["End position"]
                    if pd.isna(start) or pd.isna(end): continue
                    
                    if color_mode == "Peptide type (red/cyan)":
                        color = full_color if row["Peptide type"] == "full" else half_color
                    else:
                        fc = float(row[avg_col])
                        color = mcolors.to_hex(cmap(norm(fc)))
                    
                    segments.append({"start": int(start), "end": int(end), "color": color})

            # --- Structure Viewer ---
            st.markdown("---")
            structure_text, model_url, file_format, error_msg = download_structure(uniprot)

            if structure_text is None:
                st.error("No AlphaFold model available.")
            else:
                viewer = render_structure(structure_text, segments, file_format, plddt_coloring)
                components.html(viewer._make_html(), height=650, scrolling=True)

                if color_mode == "Peptide type (red/cyan)":
                    st.markdown(f"""
                        <div style="text-align:center;">
                            <span style="color:{full_color};">● Fully-tryptic</span> &nbsp;&nbsp;
                            <span style="color:{half_color};">● Semi-tryptic</span>
                        </div>
                        """, unsafe_allow_html=True)

                if color_mode == "Fold-change heatmap":
                    fig_cb, ax_cb = plt.subplots(figsize=(1.5, 0.1))
                    cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb, orientation="horizontal")
                    cb1.set_label("Log2FC", fontsize=6)
                    cb1.ax.tick_params(labelsize=6)  
                    st.pyplot(fig_cb, use_container_width=False)

            # --- Plots (Volcano + Abundance) ---
            st.markdown("---")
            st.subheader("Conformation & Abundance")
            
            x_all = prot_all[avg_col].astype(float)
            y_all = -np.log10(prot_all[pval_col].astype(float) + 1e-300)

            fig, (ax_volc, ax_abun) = plt.subplots(1, 2, figsize=(10, 4))
            fig.subplots_adjust(wspace=0.6)

            # Volcano
            ax_volc.scatter(x_all, y_all, color="lightgrey", alpha=0.5, s=12)
            if not peps.empty:
                x_sig = peps[avg_col].astype(float)
                y_sig = -np.log10(peps[pval_col].astype(float) + 1e-300)
                ax_volc.scatter(x_sig, y_sig, facecolors="none", edgecolors="black", linewidths=1.5, s=40)
            ax_volc.set_title("Volcano Plot")
            ax_volc.set_xlabel("Log2FC")
            ax_volc.set_ylabel("-log10(P)")
            
            # Abundance Bar
            prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]
            labels = ["Soluble", "Pellet", "Total"]
            suffix_map = {"Soluble": "soluble", "Pellet": "pellet", "Total": "total"}
            y_vals, p_vals = [], []
            for label in labels:
                suffix = suffix_map[label]
                ac, pc = f"AvgLog₂({selected_condition}).{suffix}", f"AdjPval({selected_condition}).{suffix}"
                y_vals.append(float(prot_abun[ac].values[0]) if ac in prot_abun.columns else np.nan)
                p_vals.append(prot_abun[pc].values[0] if pc in prot_abun.columns else None)
            
            x_pos = np.arange(3)
            ax_abun.bar(x_pos, np.nan_to_num(y_vals), fill=False, edgecolor="black")
            ax_abun.set_xticks(x_pos)
            ax_abun.set_xticklabels(labels)
            ax_abun.set_title("Abundance")
            
            st.pyplot(fig)
            
            st.subheader("Peptide Data")
            st.dataframe(peps[["Peptide sequence", avg_col, pval_col]])