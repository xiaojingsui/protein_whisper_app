import streamlit as st
import pandas as pd
import requests
import re
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import json

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
    /* --- 1. GLOBAL FONTS --- */
    html, body, p, h1, h2, h3, h4, h5, h6 {
        font-family: Arial, Helvetica, sans-serif !important;
    }
    
    /* --- 2. HEADER & SIDEBAR TOGGLE FIX --- */
    [data-testid="stHeader"] {
        background-color: transparent !important;
        height: 4rem !important;
        pointer-events: none !important;
        z-index: 1000001 !important;
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
        padding-top: 6rem !important; 
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

# ============================================================
# NEW: INTERACTIVE COMBINED VIEWER (Structure + Volcano)
# ============================================================
def render_interactive_dashboard(pdb_content, file_format, peptides_df):
    """
    Creates a single HTML/JS component linking 3Dmol.js (structure) and Plotly.js (volcano).
    """
    
    # 1. Prepare Peptide Data for JavaScript
    # We create a list of dicts: {index, start, end, x, y, color, seq}
    js_peptides = []
    
    # Normalize PDB content for JS injection (escape backticks)
    pdb_clean = json.dumps(pdb_content) 
    
    # Create colors for the plot
    # We'll use a simple logic here: Grey for non-sig (handled by filters), 
    # but we pass the specific color calculated in Python.
    for i, row in peptides_df.iterrows():
        js_peptides.append({
            "index": i,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "x": row["log2fc"],
            "y": row["neg_log_p"],
            "color": row["color"],
            "seq": row["sequence"]
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
                font-family: sans-serif;
                font-size: 12px;
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

        // --- 1. SETUP 3DMOL VIEWER ---
        const viewer = $3Dmol.createViewer("mol-container", {{
            defaultcolors: $3Dmol.rasmolElementColors
        }});
        
        viewer.addModel(pdbData, fileFormat);
        
        // Base Style: Transparent White Cartoon
        viewer.setStyle({{}}, {{cartoon: {{color: 'white', opacity: 0.6}}}});
        
        // Apply Peptide Colors
        peptides.forEach(p => {{
            viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}});
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
                size: 10,
                color: colors,
                line: {{color: 'black', width: 1}}
            }},
            hoverinfo: 'text+x+y'
        }};

        const layout = {{
            title: 'Volcano Plot',
            hovermode: 'closest',
            margin: {{t: 40, l: 50, r: 20, b: 40}},
            xaxis: {{title: 'Log2FC'}},
            yaxis: {{title: '-Log10(P)'}}
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
                viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}});
            }});
            viewer.render();
        }});

        // --- 4. INTERACTION: STRUCTURE -> PLOT ---
        // We use setHoverable to detect mouse over residues
        
        // This is tricky because atoms are many, points are few.
        // Logic: When hovering an atom, find if it belongs to any peptide range.
        
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
                // We use Plotly.restyle to change the color/size of the specific point
                
                // Construct new color/size arrays
                const newColors = [...colors];
                const newSizes = Array(peptides.length).fill(10);
                
                newColors[foundIdx] = '#FFFF00'; // Yellow highlight
                newSizes[foundIdx] = 20;         // Bigger size

                Plotly.restyle('plot-container', {{
                    'marker.color': [newColors],
                    'marker.size': [newSizes],
                    'marker.line.width': [Array(peptides.length).fill(1).map((v, i) => i === foundIdx ? 3 : 1)]
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
                    'marker.size': [Array(peptides.length).fill(10)],
                    'marker.line.width': [Array(peptides.length).fill(1)]
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

# SPACER for fixed navbar
st.markdown('<div style="height: 80px;"></div>', unsafe_allow_html=True)

# ============================================================
# PAGE CONTENT
# ============================================================

if page == "About":
    with st.sidebar:
        st.info("Navigate to the **Search** tab to explore proteins.")
    
    st.title("About Protein Whisper")
    st.markdown("""
    **Protein Whisper** is an interactive visualization tool designed to explore Limited Proteolysis