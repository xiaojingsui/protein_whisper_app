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
    initial_sidebar_state="expanded"
)

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
        st.error("File 'Table_S1A.xlsx' not found. Please ensure the data file is in the directory.")
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
    st.stop()

# ============================================================
# Helper Functions (Conditions, 3D, AlphaFold)
# ============================================================
def extract_conditions(df):
    conds = []
    for col in df.columns:
        m = re.match(r"AvgLog₂\((.+)\)\.conformation", col)
        if m:
            conds.append(m.group(1))
    return sorted(set(conds))

def sort_conditions(conds):
    aging = ["wt day6/wt day1", "wt day9/wt day1", "Q35/wt aging"]
    heatshock = ["wt heat-shock 35°C/wt 20°C"]
    polyq = ["Q35/Q24", "Q40/Q24"]
    myosin = ["myosin-ts 15°C/wt 15°C", "myosin-ts 25°C/wt 25°C"]
    paramyosin = ["paramyosin-ts 15°C/wt 15°C", "paramyosin-ts 25°C/wt 25°C"]
    
    priority = aging + heatshock + polyq + myosin + paramyosin
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
# Top Navigation
# ============================================================
st.title("Protein Whisper")

# Using horizontal radio buttons to act as a top navbar
# label_visibility="collapsed" hides the word "Navigation" to look cleaner
page = st.radio(
    "Navigation", 
    ["Search", "About", "Guides"], 
    horizontal=True, 
    label_visibility="collapsed"
)

st.markdown("---")

# ============================================================
# PAGE: ABOUT
# ============================================================
if page == "About":
    st.header("About Protein Whisper")
    st.markdown("""
    **Protein Whisper** is an interactive visualization tool designed to explore Limited Proteolysis-Mass Spectrometry (LiP-MS) data. 
    
    It allows researchers to map peptide-level structural alterations directly onto 3D protein structures predicted by AlphaFold, providing a bridge between quantitative proteomics and structural biology.

    ### Key Features
    * **Structure Mapping:** Visualizes peptides that undergo significant structural changes.
    * **Dual-Layer Data:** Integrates conformational data with protein solubility and total abundance.
    * **AlphaFold Integration:** Automatically fetches predicted structures from the EBI AlphaFold database.
    
    ### Data Source
    This tool utilizes data from **Table S1A**, containing LiP-MS peptide intensity data and protein-level abundance measurements across various stress conditions in *C. elegans*.
    """)

# ============================================================
# PAGE: GUIDES
# ============================================================
elif page == "Guides":
    st.header("User Guide")
    
    st.markdown("### 1. How to Search")
    st.info("Click the **Search** tab above to begin.")
    st.markdown("""
    1.  **Gene/Protein:** Enter a *C. elegans* gene symbol (e.g., `unc-54`) or a UniProt ID in the sidebar.
    2.  **Condition:** Select the stress condition you wish to compare against the wild-type/control.
    3.  **Filters:** Adjust the Fold-Change and P-value cutoffs to refine which peptides are highlighted.
    """)

    st.markdown("### 2. Understanding the Visualization")
    st.markdown("""
    * **3D Viewer:** * The protein backbone is shown in white (or colored by pLDDT confidence).
        * **Red segments** indicate fully-tryptic peptides (cleaved at both ends).
        * **Cyan segments** indicate semi-tryptic peptides (cleaved at one end), often indicative of a structural cleavage event.
    * **Volcano Plot:** Shows the statistical significance of peptide changes. Hollow circles represent significant peptides highlighted in the 3D viewer.
    * **Abundance Plot:** Displays the total protein level changes (Soluble vs. Pellet vs. Total) to help distinguish specific structural changes from general degradation or aggregation.
    """)

# ============================================================
# PAGE: SEARCH (Main Functionality)
# ============================================================
elif page == "Search":
    # --- Sidebar Controls (Specific to Search) ---
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

    # --- Main Search Logic ---
    st.header("Structure Viewer")

    if not query:
        # Intro Screen
        st.markdown(
        """
        <style>
        .stAlert { background-color: #e0faff !important; border: 1px solid #19CFE2 !important; border-radius: 8px; }
        .stAlert p { color: #0b7d8a !important; font-weight: 500; }
        </style>
        """, unsafe_allow_html=True)

        st.info("Type a C. elegans gene name or UniProt ID in the sidebar to explore a protein. Try 'unc-54'.")

        st.markdown(
            """
            <style>
            @keyframes blink { 0% { opacity: 1; } 50% { opacity: 0; } 100% { opacity: 1; } }
            .hint-container { display: flex; align-items: center; margin-top: 10px; }
            .blink-arrow { font-size: 38px; color: skyblue; font-weight: bold; animation: blink 1s infinite; margin-right: 8px; }
            .start-text { font-size: 30px; font-weight: 600; color: skyblue; }
            </style>
            <div class="hint-container">
                <div class="blink-arrow">⬅</div>
                <div class="start-text">Start your search here</div>
            </div>
            """, unsafe_allow_html=True)
    else:
        # Perform Search
        hits = df[
            df["uniprot_id"].str.contains(query, case=False, na=False)
            | df["gene"].str.contains(query, case=False, na=False)
        ]

        if hits.empty:
            st.error(f"No protein found for '{query}'.")
        else:
            protein = hits.iloc[0]
            uniprot = protein["uniprot_id"]
            protein_name = protein["Description"].split("OS=")[0].strip()
            gene = protein["gene"]

            st.subheader(f"Protein: {gene} ({uniprot}) — {protein_name}")

            # --- Filtering ---
            avg_col = f"AvgLog₂({selected_condition}).conformation"
            pval_col = f"AdjPval({selected_condition}).conformation"
            prot_all = df[df["uniprot_id"] == uniprot]
            peps = prot_all[
                (prot_all[avg_col].abs() >= fc_cutoff)
                & (prot_all[pval_col] <= p_cutoff)
            ]

            if peps.empty:
                st.warning("No peptides meet conformation filters.")
            else:
                st.success(f"{len(peps)} peptides highlighted.")

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

            # --- Structure ---
            st.markdown("---")
            st.subheader("3D Structure Viewer")
            structure_text, model_url, file_format, error_msg = download_structure(uniprot)

            if structure_text is None:
                st.error("No AlphaFold model available.")
                st.warning(error_msg)
            else:
                viewer = render_structure(structure_text, segments, file_format, plddt_coloring)
                components.html(viewer._make_html(), height=650, scrolling=True)

                if color_mode == "Peptide type (red/cyan)":
                    st.markdown(f"""
                        <div style="text-align:center; font-size:16px;">
                            <span style="color:{full_color};">● Fully-tryptic peptide</span> &nbsp;&nbsp;&nbsp;&nbsp;
                            <span style="color:{half_color};">● Semi-tryptic peptide</span>
                        </div>
                        """, unsafe_allow_html=True)

                if color_mode == "Fold-change heatmap":
                    st.markdown("**Fold-change color scale:**")
                    fig_cb, ax_cb = plt.subplots(figsize=(1.5, 0.1))
                    cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb, orientation="horizontal")
                    cb1.set_label(f"AvgLog₂({selected_condition})", fontsize=6)
                    cb1.ax.tick_params(labelsize=6)  
                    st.pyplot(fig_cb, use_container_width=False)

            # --- Plots ---
            st.markdown("---")
            st.subheader("Protein Conformation & Solubility & Total Abundance")
            
            x_all = prot_all[avg_col].astype(float)
            y_all = -np.log10(prot_all[pval_col].astype(float) + 1e-300)

            fig, (ax_volc, ax_abun) = plt.subplots(1, 2, figsize=(10, 4))
            fig.subplots_adjust(wspace=0.6)
            fig.tight_layout(pad=3.0)

            # Volcano
            ax_volc.set_title("Protein Conformation (Peptides)", fontsize=14, fontweight="bold")
            ax_volc.scatter(x_all, y_all, color="lightgrey", alpha=0.5, s=12)
            if not peps.empty:
                x_sig = peps[avg_col].astype(float)
                y_sig = -np.log10(peps[pval_col].astype(float) + 1e-300)
                ax_volc.scatter(x_sig, y_sig, facecolors="none", edgecolors="black", linewidths=1.5, s=40)

            ax_volc.axvline(fc_cutoff, color="dimgray", linestyle="--", linewidth=1)
            ax_volc.axvline(-fc_cutoff, color="dimgray", linestyle="--", linewidth=1)
            ax_volc.axhline(-np.log10(p_cutoff), color="dimgray", linestyle="--", linewidth=1)
            ax_volc.set_xlabel(f"AvgLog₂({selected_condition})")
            ax_volc.set_ylabel("-log₁₀(AdjPval)")
            ax_volc.grid(alpha=0.25)

            # Abundance
            prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]
            labels = ["Soluble", "Pellet", "Total"]
            suffix_map = {"Soluble": "soluble", "Pellet": "pellet", "Total": "total"}
            y_vals, p_vals = [], []

            for label in labels:
                suffix = suffix_map[label]
                avg_col_abun = f"AvgLog₂({selected_condition}).{suffix}"
                pval_col_abun = f"AdjPval({selected_condition}).{suffix}"
                if avg_col_abun in prot_abun.columns:
                    y_vals.append(float(prot_abun[avg_col_abun].values[0]))
                    p_vals.append(prot_abun[pval_col_abun].values[0])
                else:
                    y_vals.append(np.nan)
                    p_vals.append(None)

            x_pos = np.arange(len(labels))
            ax_abun.set_title("Protein Solubility & Abundance", fontsize=14, fontweight="bold")
            ax_abun.bar(x_pos, np.nan_to_num(y_vals, nan=0.0), width=0.66, fill=False, edgecolor="black", linewidth=1.5)
            ax_abun.set_xticks(x_pos)
            ax_abun.set_xticklabels(labels)
            ax_abun.set_ylabel(f"AvgLog₂({selected_condition})")

            finite_vals = [v for v in y_vals if np.isfinite(v)]
            if finite_vals:
                ymin, ymax = min(finite_vals), max(finite_vals)
                margin = 0.3 * (abs(ymax) + abs(ymin) + 0.1)
                ax_abun.set_ylim(min(0, ymin) - margin, max(0, ymax) + margin)
            ax_abun.axhline(0, color="dimgray", linestyle="--", linewidth=1)
            
            for xpos, yv, pv in zip(x_pos, y_vals, p_vals):
                if np.isfinite(yv) and pv is not None:
                    ax_abun.text(xpos, yv + (0.05 if yv >= 0 else -0.05), f"AdjP = {pv:.3g}", ha="center", va="bottom" if yv >= 0 else "top")

            ax_abun.grid(alpha=0.2)
            st.pyplot(fig)

            # --- Table ---
            st.subheader("Peptides Passing Conformation Filters")
            st.dataframe(peps[["Peptide sequence", "Start position", "End position", avg_col, pval_col]])
            