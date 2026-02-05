import streamlit as st
import pandas as pd
import requests
import re
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import json
import os
from PIL import Image

# ============================================================
# MATPLOTLIB FONT CONFIGURATION (Force Arial for static plots)
# ============================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ============================================================
# Page Configuration
# ============================================================
st.set_page_config(
    page_title="Protein Whisper",
    layout="wide",
    initial_sidebar_state="expanded" 
)

# ============================================================
# GLOBAL CSS & NAVIGATION STYLING
# ============================================================
st.markdown("""
    <style>
    /* --- 1. GLOBAL FONTS (FIXED) --- */
    html, body, p, h1, h2, h3, h4, h5, h6, li, a, label, button, input, select, textarea,
    [data-testid="stTextInput"] input {  /* <--- ADD THIS LINE */
        font-family: Arial, Helvetica, sans-serif !important;
    }
    
    /* --- 2. HEADER & SIDEBAR TOGGLE FIX --- */
    [data-testid="stHeader"] {
        background-color: transparent !important;
        height: 4rem !important;
        pointer-events: none !important;
        z-index: 1000001 !important;
    }
    .stTextInput > div > div > input,
    [data-testid="stNumberInput"] input {
        font-family: Arial, Helvetica, sans-serif !important;
        font-size: 14px !important;  /* <--- Adds consistency with labels/dropdowns */
    }
    
    [data-testid="stHeader"] button[title="View fullscreen"], 
    [data-testid="stHeader"] button[kind="header"] {
        pointer-events: auto !important;
        color: #445550 !important;
    }
    
    [data-testid="stDecoration"], [data-testid="stHeader"] > div:last-child {
        display: none !important;
    }

    /* --- 3. LAYOUT PADDING --- */
    .block-container {
        padding-top: 0rem !important; 
        padding-bottom: 5rem !important;
    }

    /* --- 4. CUSTOM NAVBAR --- */
    div[role="radiogroup"] {
        position: fixed !important;        
        top: 0 !important;                  
        left: 0 !important;                 
        width: 100vw !important;            
        z-index: 1000000 !important;    
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
    
    div[role="radiogroup"] label[data-testid="stRadioOption"] {
        background-color: transparent !important;
    }

    /* --- 5. LOGO STYLING --- */
    .logo-text {
        position: fixed;
        top: 12px;
        left: 25px;
        font-size: 22px;
        font-weight: bold;
        color: #445550;
        z-index: 1000002;
        font-family: Arial, Helvetica, sans-serif;
        pointer-events: none;
    }

    /* --- 6. TEXT-ONLY LINK BUTTONS --- */
    div[data-testid="stSidebar"] .stButton > button {
        background-color: transparent !important;
        border: none !important;
        box-shadow: none !important;
        color: #006064 !important;
        text-decoration: underline !important;
        padding: 0px !important;
        font-family: Arial, Helvetica, sans-serif !important;
    }

    div[data-testid="stSidebar"] .stButton > button:hover {
        color: #19CFE2 !important;       
        text-decoration: none !important;
    }
    
    div[data-testid="stSidebar"] .stButton > button:focus,
    div[data-testid="stSidebar"] .stButton > button:active {
        border: none !important;
        box-shadow: none !important;
        outline: none !important;
        color: #006064 !important;
        background-color: transparent !important;
    }
    </style>
""", unsafe_allow_html=True)

# ============================================================
# INJECT LOGO
# ============================================================
st.markdown('<div class="logo-text">ProteinWhisper</div>', unsafe_allow_html=True)

# ============================================================
# Data Loading & Caching
# ============================================================
@st.cache_data
def load_table():
    try:
        df = pd.read_excel("Table_S1A.xlsx", header=3)
        df["uniprot_id"] = df["Protein ID"].astype(str).str.split("|").str[1].str.strip()
        
        # Force Uppercase for Gene Symbol
        df["gene"] = df["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip().str.upper()
        
        # --- CLEANING DESCRIPTION LOGIC ---
        if "Protein names" in df.columns:
            df["desc"] = df["Protein names"].astype(str).str.split(";").str[0].str.split(" OS=").str[0].str.strip()
        elif "Description" in df.columns:
            df["desc"] = df["Description"].astype(str).str.split(" OS=").str[0].str.strip()
        else:
            df["desc"] = ""
            
        return df
    except FileNotFoundError:
        return pd.DataFrame()

@st.cache_data
def load_abundance_table():
    try:
        df2 = pd.read_excel("Table_S1A.xlsx", sheet_name="Protein-level data", header=3)
        df2["uniprot_id"] = df2["Protein ID"].astype(str).str.split("|").str[1].str.strip()
        df2["gene"] = df2["Gene symbol"].astype(str).str.replace("CELE_", "", regex=False).str.strip().str.upper()
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

# ============================================================
# NEW: INTERACTIVE COMBINED VIEWER (Structure + Volcano)
# ============================================================
def render_interactive_dashboard(pdb_content, file_format, peptides_df, fc_cutoff, p_cutoff, selected_condition):
    """
    Creates a single HTML/JS component linking 3Dmol.js (structure) and Plotly.js (volcano).
    """
    
    # 1. Prepare Peptide Data for JavaScript
    js_peptides = []
    
    # Normalize PDB content for JS injection (escape backticks)
    pdb_clean = json.dumps(pdb_content) 
    
    # Determine plot threshold line values for JS
    neg_log_p_thresh = -np.log10(p_cutoff)

    for i, row in peptides_df.iterrows():
        js_peptides.append({
            "index": i,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "x": row["log2fc"],
            "y": row["neg_log_p"],
            "color": row["color"], # Significant gets color, others grey
            "seq": row["sequence"],
            "is_sig": row["is_significant"]
        })

    js_peptides_json = json.dumps(js_peptides)

    # 2. The HTML/JS Template
    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            /* Force Arial in the HTML Component */
            body, .tooltip, .js-plotly-plot {{
                font-family: Arial, Helvetica, sans-serif !important;
            }}
            .container {{
                display: flex;
                flex-direction: row;
                width: 100%;
                height: 600px;
                gap: 20px;
            }}
            #mol-container {{
                width: 60%;
                height: 100%;
                border: 1px solid #eee;
                position: relative;
            }}
            #plot-container {{
                width: 40%;
                height: 100%;
                border: 1px solid #eee;
            }}
            .tooltip {{
                position: absolute;
                background: rgba(0, 0, 0, 0.7);
                color: white;
                padding: 5px;
                border-radius: 4px;
                pointer-events: none;
                display: none;
                z-index: 100;
                font-size: 12px;
                font-family: Arial, Helvetica, sans-serif; /* Explicit Font */
            }}
        </style>
    </head>
    <body>

    <div class="container">
        <div id="mol-container"></div>
        <div id="plot-container"></div>
    </div>
    <div id="tooltip" class="tooltip"></div>

    <script>
        // --- DATA INJECTION ---
        const pdbData = {pdb_clean};
        const fileFormat = "{file_format}";
        const peptides = {js_peptides_json};
        const fcThresh = {fc_cutoff};
        const pThresh = {neg_log_p_thresh};
        const xLabel = "Log<sub>2</sub> ({selected_condition})";

        // --- 1. SETUP 3DMOL VIEWER ---
        const viewer = $3Dmol.createViewer("mol-container", {{
            defaultcolors: $3Dmol.rasmolElementColors
        }});
        
        viewer.addModel(pdbData, fileFormat);
        
        // Base Style: Transparent White Cartoon
        viewer.setStyle({{}}, {{cartoon: {{color: 'white', opacity: 0.6}}}});
        
        // Apply Peptide Colors ONLY for Significant Peptides
        peptides.forEach(p => {{
            if (p.is_sig) {{
                viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}});
            }}
        }});

        viewer.zoomTo();
        viewer.render();

        // Helper to handle range lists for 3Dmol
        function splitRange(start, end) {{
            let arr = [];
            for (let i = start; i <= end; i++) arr.push(i);
            return arr;
        }}

        // --- 2. SETUP PLOTLY VOLCANO ---
        const xVals = peptides.map(p => p.x);
        const yVals = peptides.map(p => p.y);
        const colors = peptides.map(p => p.color);
        const hoverText = peptides.map(p => `Seq: ${{p.seq}}<br>Start: ${{p.start}}<br>End: ${{p.end}}`);

        const trace = {{
            x: xVals,
            y: yVals,
            mode: 'markers',
            type: 'scatter',
            text: hoverText,
            marker: {{
                size: 8,
                color: colors,
                line: {{color: 'black', width: 0.5}},
                opacity: 0.8
            }},
            hoverinfo: 'text+x+y'
        }};

        const layout = {{
            title: 'PK Accessibility',
            font: {{ family: 'Arial, Helvetica, sans-serif' }}, // Force Plotly Font
            hovermode: 'closest',
            margin: {{t: 40, l: 50, r: 20, b: 40}},
            xaxis: {{
                title: xLabel,
                zeroline: false  // Removes the default x=0 line
            }},
            yaxis: {{
                title: '-Log<sub>10</sub> (FDR)',
                rangemode: 'tozero' // Forces Y-axis to start at 0
            }},
            shapes: [
                // Vertical Line +FC
                {{
                    type: 'line', x0: fcThresh, x1: fcThresh, y0: 0, y1: 1, yref: 'paper',
                    line: {{color: 'grey', width: 1.5, dash: 'dash'}}
                }},
                // Vertical Line -FC
                {{
                    type: 'line', x0: -fcThresh, x1: -fcThresh, y0: 0, y1: 1, yref: 'paper',
                    line: {{color: 'grey', width: 1.5, dash: 'dash'}}
                }},
                // Horizontal Line P-val
                {{
                    type: 'line', x0: 0, x1: 1, xref: 'paper', y0: pThresh, y1: pThresh,
                    line: {{color: 'grey', width: 1.5, dash: 'dash'}}
                }}
            ]
        }};

        Plotly.newPlot('plot-container', [trace], layout);

        // --- 3. INTERACTION: PLOT -> STRUCTURE ---
        const plotDiv = document.getElementById('plot-container');

        plotDiv.on('plotly_hover', function(data){{
            const pt = data.points[0];
            const idx = pt.pointIndex; // Index in our peptides array
            const pep = peptides[idx];

            // Highlight in 3D
            // 1. Add a temporary highlight style (Thicker, brighter)
            viewer.addStyle({{resi: splitRange(pep.start, pep.end)}}, 
                {{cartoon: {{color: '#FFFF00', thickness: 1.0, opacity: 1.0}}}}
            );
            viewer.render();
        }});

        plotDiv.on('plotly_unhover', function(data){{
            // Reset Styles
            viewer.setStyle({{}}, {{cartoon: {{color: 'white', opacity: 0.6}}}});
            peptides.forEach(p => {{
                if(p.is_sig) {{
                    viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}});
                }}
            }});
            viewer.render();
        }});

        // --- 4. INTERACTION: STRUCTURE -> PLOT ---
        let lastHoveredIdx = -1;

        viewer.setHoverable({{}}, true, function(atom, viewer, event, container) {{
            if(!atom) return;
            
            // Find which peptide contains this residue
            let foundIdx = -1;
            for(let i=0; i<peptides.length; i++) {{
                if(atom.resi >= peptides[i].start && atom.resi <= peptides[i].end) {{
                    foundIdx = i;
                    break;
                }}
            }}

            if (foundIdx !== -1 && foundIdx !== lastHoveredIdx) {{
                lastHoveredIdx = foundIdx;
                
                // Highlight Plot Point
                const newColors = [...colors];
                const newSizes = Array(peptides.length).fill(8);
                
                newColors[foundIdx] = '#FFFF00'; // Yellow highlight
                newSizes[foundIdx] = 18;         // Bigger size

                Plotly.restyle('plot-container', {{
                    'marker.color': [newColors],
                    'marker.size': [newSizes],
                    'marker.line.width': [Array(peptides.length).fill(0.5).map((v, i) => i === foundIdx ? 2 : 0.5)]
                }});
                
                // Show custom tooltip on structure
                const tooltip = document.getElementById('tooltip');
                tooltip.style.display = 'block';
                tooltip.style.left = event.x + 'px';
                tooltip.style.top = event.y + 'px';
                tooltip.innerHTML = 'Residue: ' + atom.resi + '<br>Seq: ' + peptides[foundIdx].seq;
            }}
        }}, function(atom, viewer, event, container) {{
            // On Leave
            if (lastHoveredIdx !== -1) {{
                // Reset Plot
                Plotly.restyle('plot-container', {{
                    'marker.color': [colors],
                    'marker.size': [Array(peptides.length).fill(8)],
                    'marker.line.width': [Array(peptides.length).fill(0.5)]
                }});
                lastHoveredIdx = -1;
                document.getElementById('tooltip').style.display = 'none';
            }}
        }});

    </script>
    </body>
    </html>
    """
    return html_code

# ============================================================
# NAVIGATION CONTROLLER
# ============================================================
page = st.radio(
    "Main Navigation", 
    ["Search", "About", "Guides"], 
    horizontal=True,
    label_visibility="collapsed"
)

# ============================================================
# DYNAMIC SIDEBAR VISIBILITY LOGIC
# ============================================================
# If the user is on 'About' or 'Guides', we inject CSS to completely hide the sidebar.
if page in ["About", "Guides"]:
    st.markdown("""
        <style>
        [data-testid="stSidebar"] {
            display: none !important;
        }
        [data-testid="stSidebarCollapsedControl"] {
            display: none !important;
        }
        </style>
    """, unsafe_allow_html=True)

# SPACER for fixed navbar
st.markdown('<div style="height: 80px;"></div>', unsafe_allow_html=True)

# ============================================================
# PAGE CONTENT
# ============================================================

if page == "About":
    st.markdown("<h1>About Protein Whisper</h1>", unsafe_allow_html=True)
    st.markdown("""
    **Protein Whisper** is an interactive visualization tool designed to explore TMT-based Limited Proteolysis-Mass Spectrometry (TMT-LiP-MS) data. It allows researchers to map peptide-level structural alterations directly onto AlphaFold2-predicted 3D protein structures.
    """)
    
    # --- Workflow Image (Half Size using columns) ---
    if os.path.exists("workflow.png"):
        c_img, c_spacer = st.columns([1, 1])
        with c_img:
            st.image("workflow.png", caption="Protein Whisper Workflow", use_container_width=True)
    
    st.markdown("---")
    
    # --- Table Information (Detailed) ---
    st.markdown("### Dataset Overview")
    st.markdown("""
    The data summarizes integrated proteomics measurements across aging, myosin-ts, paramyosin-ts, polyQ, and heat-shock conditions at peptide levels. 
    Each row corresponds to an identified peptide mapped to its parent protein, with experimental log₂-fold changes derived from TMT-based quantitative proteomics.
    """)
    
    c1, c2 = st.columns(2)
    
    with c1:
        st.markdown("**Experimental Conditions:**")
        st.markdown("""
        1. **Myosin-ts:** 15 °C and 25 °C
        2. **Paramyosin-ts:** 15 °C and 25 °C
        3. **PolyQ:** Q24/35/40
        4. **Heat-shock:** 35 °C vs 20 °C
        5. **WT and Q35 Aging:** Day 1, day 6, and day 9 time points
        """)
        
        st.markdown("**Fractions:**")
        st.markdown("""
        - **Conformation (PK):** Protease accessibility (structural) measurements.
        - **Solubility:** Soluble and Pellet protein abundance measurements.
        - **Total Abundance:** Total protein abundance measurements.
        """)

    with c2:
        st.markdown("**Column Descriptions:**")
        st.markdown("""
        - **Log₂(…)**: Replicate-level changes in abundance or PK accessibility relative to the indicated control.
        - **AvgLog₂(…)**: Average log₂-fold change across biological replicates.
        - **Pval(…) and AdjPval(…)**: Statistical significance and FDR-adjusted values, respectively.
        """)
    
    # NOTE: The Color Code Legend has been removed as requested.


elif page == "Guides":
    st.title("User Guide")

    st.markdown("### 1. Navigation")
    st.info("Click the **Search** tab in the top navigation bar to begin your analysis.")

    st.markdown("---")

    st.markdown("### 2. Configuring Your Search")
    st.markdown("""
    Once on the Search page, use the **Sidebar** (left panel) to control the analysis:
    
    1.  **Search Query:** Enter a *C. elegans* gene symbol (e.g., `UNC-54`) or a UniProt ID.
    2.  **Stress Condition:** Select the specific experimental comparison (e.g., *wt heat-shock* or *aging*).
    3.  **Filter Settings:**
        * **Fold-change cutoff:** Set the threshold for magnitude of change (default: 1.0).
        * **AdjPval cutoff:** Set the threshold for statistical significance (default: 0.05).
    4.  **Visualization:** Toggle between "Heatmap" (intensity) or "Peptide type" coloring modes.
    """)

    st.markdown("---")

    st.markdown("### 3. Interactive Dashboard")
    st.markdown("""
    The main view connects structural locations with mass spectrometry data:

    * **The Structure (Left):** An AlphaFold2 predicted model of the protein.
    * **The Volcano Plot (Right):** * **Dots:** Represent individual peptides.
        * **Grey dots:** Non-significant peptides (below your filter settings).
        * **Colored dots:** Significant peptides.
    
    **How to Interact:**
    * **Hover over the Plot:** Hovering over a dot highlights that specific peptide sequence in **Yellow** on the 3D structure.
    * **Hover over the Structure:** Hovering over the protein residues highlights the corresponding data point on the volcano plot (the dot becomes larger).
    """)

    st.markdown("---")

    st.markdown("### 4. Abundance Data")
    st.markdown("""
    Scroll to the bottom of the Search page to view static bar charts:
    * **Solubility Changes:** Compares abundance in Soluble vs. Pellet fractions.
    * **Total Abundance:** Shows changes in the overall protein levels.
    * *P-values represent the statistical significance of these specific changes.*
    """)

elif page == "Search":
    # --- Sidebar for inputs (ONLY Visible here) ---
    with st.sidebar:
        st.header("Search Parameters")
        
        # Session State Logic - Default to UNC-54
        if 'search_term' not in st.session_state:
            st.session_state.search_term = "UNC-54"

        # 1. Search Box
        query = st.text_input("Search gene or UniProt ID:", key="search_term")
        st.write("") 

        # 3. Rest of controls
        selected_condition = st.selectbox("Stress condition:", conditions)
        
        with st.expander("Filter Settings", expanded=True):
            fc_cutoff = st.number_input(r"Fold-change cutoff (|Log$_2$|):", value=1.0, min_value=0.0, step=0.1)
            p_cutoff = st.number_input("AdjPval cutoff:", value=0.05, min_value=0.0, max_value=1.0, step=0.01)

        with st.expander("Visualization Settings"):
            color_mode = st.selectbox("Peptide color mode:", ["Fold-change heatmap", "Peptide type (red/cyan)"])

    # --- Main Content ---
    if not query:
        st.markdown("<h2 style='text-align: center; color: #19CFE2;'>Structure Viewer</h2>", unsafe_allow_html=True)
        st.info("Type a C. elegans gene name or UniProt ID in the **sidebar** to explore a protein.")
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
            
            # Retrieve Description safely
            desc = protein["desc"] if "desc" in protein else ""
            
            # --- Header Display ---
            if desc:
                st.markdown(f"### Protein: **{gene}** ({desc}) ([{uniprot}](https://www.uniprot.org/uniprot/{uniprot}))")
            else:
                st.markdown(f"### Protein: **{gene}** ([{uniprot}](https://www.uniprot.org/uniprot/{uniprot}))")

            st.markdown("#### Protein Conformational Changes")
            
            # --- Filtering ---
            avg_col = f"AvgLog₂({selected_condition}).conformation"
            pval_col = f"AdjPval({selected_condition}).conformation"
            prot_all = df[df["uniprot_id"] == uniprot]
            
            # Identify Significant vs Non-Significant
            peps = prot_all.copy()
            peps["is_sig"] = (peps[avg_col].abs() >= fc_cutoff) & (peps[pval_col] <= p_cutoff)
            
            sig_count = peps["is_sig"].sum()
            
            if peps.empty:
                  st.warning("No peptides found for this protein.")
            else:
                  if sig_count > 0:
                      st.success(f"{sig_count} significant peptides found out of {len(peps)} total.")
                  else:
                      st.info(f"No significant peptides found (Total: {len(peps)}). Adjust filters to see changes.")

            # --- Prepare Data for Visualization ---
            viz_data = []
            full_color = "#E1341E"
            half_color = "#1ECBE1"
            grey_color = "#D3D3D3" # Light grey for non-sig

            if not peps.empty:
                # Normalization for heatmap (Calculated on significant range usually, or global)
                fc_vals = peps[avg_col].astype(float)
                maxfc = float(fc_vals.abs().max())
                if maxfc == 0: maxfc = 1.0
                
                cmap = plt.get_cmap("coolwarm")
                norm = mcolors.TwoSlopeNorm(vmin=-maxfc, vcenter=0, vmax=maxfc)

                for _, row in peps.iterrows():
                    start, end = row["Start position"], row["End position"]
                    if pd.isna(start) or pd.isna(end): continue
                    
                    fc = float(row[avg_col])
                    pval = float(row[pval_col])
                    seq = row["Peptide sequence"]
                    is_sig = row["is_sig"]

                    # Determine Color
                    if not is_sig:
                        color = grey_color
                    else:
                        if color_mode == "Peptide type (red/cyan)":
                            color = full_color if row["Peptide type"] == "full" else half_color
                        else:
                            color = mcolors.to_hex(cmap(norm(fc)))
                    
                    viz_data.append({
                        "start": start, 
                        "end": end, 
                        "color": color, 
                        "log2fc": fc,
                        "neg_log_p": -np.log10(pval + 1e-300),
                        "sequence": seq,
                        "is_significant": bool(is_sig)
                    })
            
            viz_df = pd.DataFrame(viz_data)

            # --- Download Structure ---
            st.markdown("---")
            structure_text, model_url, file_format, error_msg = download_structure(uniprot)

            if structure_text is None:
                st.error("No AlphaFold model available.")
            elif viz_df.empty:
                st.warning("Structure available, but no peptides to map.")
            else:
                # --- RENDER INTERACTIVE DASHBOARD ---
                # Pass thresholds and condition name to function
                html_view = render_interactive_dashboard(
                    structure_text, file_format, viz_df, fc_cutoff, p_cutoff, selected_condition
                )
                components.html(html_view, height=650)
                
                # Legend for Colors (Only show if Heatmap is active and there are sig peptides)
                if color_mode == "Fold-change heatmap" and sig_count > 0:
                    c_left, c_rest = st.columns([1, 2])

                    with c_left:
                          fig_cb, ax_cb = plt.subplots(figsize=(1, 0.1))
                          cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb, orientation="horizontal")
                          cb1.set_label(r"Log$_2$FC", fontsize=6, fontname='Arial')
                          cb1.ax.tick_params(labelsize=6)
                          for label in cb1.ax.get_xticklabels():
                              label.set_fontname('Arial')
                          st.pyplot(fig_cb, use_container_width=False)
                
                st.info("💡 **Interactive:** Hover over a dot on the volcano plot to see it on the structure. Grey dots are non-significant.")

            # --- Abundance Plot ---
            st.markdown("---")
            
            # 1. FORCE ARIAL GLOBALLY
            plt.rcParams.update({
                'font.family': 'sans-serif',
                'font.sans-serif': ['Arial'],
                'axes.unicode_minus': False
            })
            
            prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]

            if not prot_abun.empty:
                # 2. Helper to fetch Data
                def get_data(suffix):
                    fc_col = f"AvgLog₂({selected_condition}).{suffix}"
                    fc = 0.0
                    if fc_col in prot_abun.columns and pd.notna(prot_abun[fc_col].values[0]):
                        fc = float(prot_abun[fc_col].values[0])
                    
                    pval_col = f"AdjPval({selected_condition}).{suffix}"
                    pval = 1.0 
                    if pval_col in prot_abun.columns and pd.notna(prot_abun[pval_col].values[0]):
                        pval = float(prot_abun[pval_col].values[0])
                    
                    return fc, pval

                fc_sol, p_sol = get_data("soluble")
                fc_pel, p_pel = get_data("pellet")
                fc_tot, p_tot = get_data("total")

                # 3. Setup Columns
                ab_col1, ab_col2 = st.columns(2)

                # 4. Define Dynamic Y-Axis Label
                dynamic_ylabel = f"Log\u2082 ({selected_condition})"

                # 5. Plotting Helper
                def plot_with_pval(ax, x, y, pvals, labels, ylabel_text=None, bar_width=0.6):
                    # Draw Bars
                    bars = ax.bar(x, y, fill=False, edgecolor="black", width=bar_width)
                    ax.axhline(0, color='grey', linewidth=0.8)
                    
                    # X-Axis Styling (Force Arial)
                    ax.set_xticks(x)
                    ax.set_xticklabels(labels, fontsize=8, fontname='Arial')
                    
                    # Y-Axis Label (Force Arial)
                    if ylabel_text:
                        ax.set_ylabel(ylabel_text, fontsize=8, fontname='Arial')
                    
                    # General Tick Styling (Force Arial)
                    ax.tick_params(axis='both', labelsize=8)
                    for label in ax.get_yticklabels() + ax.get_xticklabels():
                        label.set_fontname('Arial')

                    # --- DYNAMIC TEXT OFFSET ---
                    if len(y) > 0:
                        y_max = max(max(y), 0)
                        y_min = min(min(y), 0)
                        y_range = abs(y_max - y_min)
                        if y_range == 0: y_range = 1.0 
                        offset = y_range * 0.05
                    else:
                        offset = 0.1

                    # Add P-value Text (Force Arial)
                    for bar, p in zip(bars, pvals):
                        height = bar.get_height()
                        p_str = f"p={p:.1e}" if p < 0.001 else f"p={p:.3f}"
                        
                        if height >= 0:
                            y_pos = height + offset
                            va = 'bottom'
                        else:
                            y_pos = height - offset
                            va = 'top'
                        
                        ax.text(bar.get_x() + bar.get_width()/2, y_pos, p_str,
                                ha='center', va=va, fontsize=7, fontname='Arial', color='black')
                    
                    # Ensure Auto-scaling
                    ax.autoscale(enable=True, axis='y', tight=False)
                    ax.margins(y=0.25)

                # --- PLOT 1: Solubility Changes ---
                with ab_col1:
                    st.markdown("#### Solubility Changes") 
                    
                    fig1, ax1 = plt.subplots(figsize=(1.8, 1.8))
                    
                    plot_with_pval(ax1, [0, 1], [fc_sol, fc_pel], [p_sol, p_pel], 
                                   ["Soluble", "Pellet"], ylabel_text=dynamic_ylabel, bar_width=0.6)
                    
                    st.pyplot(fig1, use_container_width=False)

                # --- PLOT 2: Total Abundance Changes ---
                with ab_col2:
                    # 1. Title
                    st.markdown("#### Total Abundance Changes")
                    
                    # 2. Create Plot
                    fig2, ax2 = plt.subplots(figsize=(1, 1.8))
                    
                    # 3. Plot Data
                    plot_with_pval(ax2, [0], [fc_tot], [p_tot], 
                                   ["Total"], ylabel_text=dynamic_ylabel, bar_width=0.5)

                    # 4. CRITICAL FIX: Set x-limits to add whitespace
                    ax2.set_xlim(-1, 1)

                    st.pyplot(fig2, use_container_width=False)

            else:
                st.write("No abundance data available.")

# ============================================================
# FOOTER / CONTACT INFO (Added to Right Bottom Corner)
# ============================================================
st.write("")
st.write("")
st.write("")

# Layout columns to position the box on the far right
# 3:1 Ratio creates a "bottom right corner" effect
col_space, col_contact = st.columns([1.5, 3]) 

with col_contact:
    st.markdown("""
        <div style="
            border-left: 5px solid #006064; 
            background-color: transparent; 
            padding: 8px 12px; 
            font-family: Arial, sans-serif;
        ">
            <strong style="color: #006064; display: block; margin-bottom: 4px;">APP Support</strong>
            <div style="font-size: 15px; color: #333;">
                Xiaojing Sui: <a href="mailto:xiaojing.sui@northwestern.edu" style="color: #006064; text-weight: bold;">xiaojing.sui@northwestern.edu</a>
            </div>
        </div>
    """, unsafe_allow_html=True)