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
# MATPLOTLIB FONT CONFIGURATION (Strict Arial Enforcement)
# ============================================================
# 1. Force the base font family
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# 2. CRITICAL: Make Math text (like the subscript 2) use the same Arial font
# Without this, "Log2FC" will use a different "Math" font for the number 2.
plt.rcParams['mathtext.default'] = 'regular' 

# ============================================================
# Page Configuration
# ============================================================
st.set_page_config(
    page_title="Protein Whisper",
    layout="wide",
    initial_sidebar_state="expanded" 
)

# ============================================================
# GLOBAL CSS (FIXED TO SAVE ICONS)
# ============================================================
st.markdown("""
    <style>
    /* --- 1. TARGETED FONT APPLICATION --- */
    /* The previous error was caused by targeting 'div' and 'span' globally.
       Streamlit uses specific font-family settings on spans to render icons.
       We must only target TEXT containers to avoid breaking the "arrow_drop_down" icons.
    */
    html, body, p, h1, h2, h3, h4, h5, h6, li, a, label, button, input, select, textarea {
        font-family: Arial, Helvetica, sans-serif !important;
    }
    
    /* Target the specific Streamlit markdown container, but not the icon containers */
    .stMarkdown, .stText {
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

# ... [KEEP ALL THE LOADING FUNCTIONS AS THEY WERE] ...
# (load_table, load_abundance_table, extract_conditions, sort_conditions, download_structure)

# ============================================================
# NEW: INTERACTIVE COMBINED VIEWER (Structure + Volcano)
# ============================================================
def render_interactive_dashboard(pdb_content, file_format, peptides_df, fc_cutoff, p_cutoff, selected_condition):
    # ... [Keep this function exactly the same] ...
    # Note: The HTML inside this function already has: 
    # font-family: Arial, Helvetica, sans-serif !important;
    # So it is safe.
    
    # 1. Prepare Peptide Data for JavaScript
    js_peptides = []
    pdb_clean = json.dumps(pdb_content) 
    neg_log_p_thresh = -np.log10(p_cutoff)

    for i, row in peptides_df.iterrows():
        js_peptides.append({
            "index": i,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "x": row["log2fc"],
            "y": row["neg_log_p"],
            "color": row["color"],
            "seq": row["sequence"],
            "is_sig": row["is_significant"]
        })

    js_peptides_json = json.dumps(js_peptides)

    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body, .tooltip, .js-plotly-plot {{
                font-family: Arial, Helvetica, sans-serif !important;
            }}
            .container {{ display: flex; flex-direction: row; width: 100%; height: 600px; gap: 20px; }}
            #mol-container {{ width: 60%; height: 100%; border: 1px solid #eee; position: relative; }}
            #plot-container {{ width: 40%; height: 100%; border: 1px solid #eee; }}
            .tooltip {{
                position: absolute; background: rgba(0, 0, 0, 0.7); color: white; padding: 5px;
                border-radius: 4px; pointer-events: none; display: none; z-index: 100;
                font-size: 12px; font-family: Arial, Helvetica, sans-serif;
            }}
        </style>
    </head>
    <body>
    <div class="container"><div id="mol-container"></div><div id="plot-container"></div></div>
    <div id="tooltip" class="tooltip"></div>
    <script>
        const pdbData = {pdb_clean};
        const fileFormat = "{file_format}";
        const peptides = {js_peptides_json};
        const fcThresh = {fc_cutoff};
        const pThresh = {neg_log_p_thresh};
        const xLabel = "Log<sub>2</sub> ({selected_condition})";

        const viewer = $3Dmol.createViewer("mol-container", {{ defaultcolors: $3Dmol.rasmolElementColors }});
        viewer.addModel(pdbData, fileFormat);
        viewer.setStyle({{}}, {{cartoon: {{color: 'white', opacity: 0.6}}}});
        peptides.forEach(p => {{ if (p.is_sig) {{ viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}}); }} }});
        viewer.zoomTo(); viewer.render();

        function splitRange(start, end) {{ let arr = []; for (let i = start; i <= end; i++) arr.push(i); return arr; }}

        const xVals = peptides.map(p => p.x);
        const yVals = peptides.map(p => p.y);
        const colors = peptides.map(p => p.color);
        const hoverText = peptides.map(p => `Seq: ${{p.seq}}<br>Start: ${{p.start}}<br>End: ${{p.end}}`);

        const trace = {{
            x: xVals, y: yVals, mode: 'markers', type: 'scatter', text: hoverText,
            marker: {{ size: 8, color: colors, line: {{color: 'black', width: 0.5}}, opacity: 0.8 }},
            hoverinfo: 'text+x+y'
        }};

        const layout = {{
            title: 'PK Accessibility',
            font: {{ family: 'Arial, Helvetica, sans-serif' }}, 
            hovermode: 'closest', margin: {{t: 40, l: 50, r: 20, b: 40}},
            xaxis: {{ title: xLabel, zeroline: false }},
            yaxis: {{ title: '-Log<sub>10</sub> (FDR)', rangemode: 'tozero' }},
            shapes: [
                {{ type: 'line', x0: fcThresh, x1: fcThresh, y0: 0, y1: 1, yref: 'paper', line: {{color: 'grey', width: 1.5, dash: 'dash'}} }},
                {{ type: 'line', x0: -fcThresh, x1: -fcThresh, y0: 0, y1: 1, yref: 'paper', line: {{color: 'grey', width: 1.5, dash: 'dash'}} }},
                {{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: pThresh, y1: pThresh, line: {{color: 'grey', width: 1.5, dash: 'dash'}} }}
            ]
        }};

        Plotly.newPlot('plot-container', [trace], layout);

        const plotDiv = document.getElementById('plot-container');
        plotDiv.on('plotly_hover', function(data){{
            const idx = data.points[0].pointIndex; const pep = peptides[idx];
            viewer.addStyle({{resi: splitRange(pep.start, pep.end)}}, {{cartoon: {{color: '#FFFF00', thickness: 1.0, opacity: 1.0}}}});
            viewer.render();
        }});
        plotDiv.on('plotly_unhover', function(data){{
            viewer.setStyle({{}}, {{cartoon: {{color: 'white', opacity: 0.6}}}});
            peptides.forEach(p => {{ if(p.is_sig) {{ viewer.addStyle({{resi: splitRange(p.start, p.end)}}, {{cartoon: {{color: p.color, opacity: 1.0}}}}); }} }});
            viewer.render();
        }});

        let lastHoveredIdx = -1;
        viewer.setHoverable({{}}, true, function(atom, viewer, event, container) {{
            if(!atom) return;
            let foundIdx = -1;
            for(let i=0; i<peptides.length; i++) {{ if(atom.resi >= peptides[i].start && atom.resi <= peptides[i].end) {{ foundIdx = i; break; }} }}
            if (foundIdx !== -1 && foundIdx !== lastHoveredIdx) {{
                lastHoveredIdx = foundIdx;
                const newColors = [...colors]; const newSizes = Array(peptides.length).fill(8);
                newColors[foundIdx] = '#FFFF00'; newSizes[foundIdx] = 18;
                Plotly.restyle('plot-container', {{ 'marker.color': [newColors], 'marker.size': [newSizes], 'marker.line.width': [Array(peptides.length).fill(0.5).map((v, i) => i === foundIdx ? 2 : 0.5)] }});
                const tooltip = document.getElementById('tooltip');
                tooltip.style.display = 'block'; tooltip.style.left = event.x + 'px'; tooltip.style.top = event.y + 'px';
                tooltip.innerHTML = 'Residue: ' + atom.resi + '<br>Seq: ' + peptides[foundIdx].seq;
            }}
        }}, function(atom, viewer, event, container) {{
            if (lastHoveredIdx !== -1) {{
                Plotly.restyle('plot-container', {{ 'marker.color': [colors], 'marker.size': [Array(peptides.length).fill(8)], 'marker.line.width': [Array(peptides.length).fill(0.5)] }});
                lastHoveredIdx = -1; document.getElementById('tooltip').style.display = 'none';
            }}
        }});
    </script>
    </body>
    </html>
    """
    return html_code

# ... [KEEP NAVIGATION AND SEARCH LOGIC UNCHANGED] ...
# (Page logic, search logic, data processing)

# ... [SCROLL DOWN TO WHERE THE PLOTS ARE RENDERED] ...

            # [Inside the Search Loop, after rendering the Structure]

                # Legend for Colors (Only show if Heatmap is active and there are sig peptides)
                if color_mode == "Fold-change heatmap" and sig_count > 0:
                    c_left, c_rest = st.columns([1, 2])

                    with c_left:
                          fig_cb, ax_cb = plt.subplots(figsize=(1, 0.1))
                          cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb, orientation="horizontal")
                          
                          # --- REVISED: Explicitly set Arial for the labels ---
                          cb1.set_label(r"Log$_2$FC", fontsize=8, fontname='Arial') 
                          cb1.ax.tick_params(labelsize=8)
                          
                          # Force the ticks to use Arial explicitly (double insurance)
                          for t in cb1.ax.get_xticklabels():
                              t.set_fontname('Arial')
                              
                          st.pyplot(fig_cb, use_container_width=False)
                
                st.info("💡 **Interactive:** Hover over a dot on the volcano plot to see it on the structure. Grey dots are non-significant.")

            # --- Abundance Plot ---
            st.subheader("Protein Abundance")
            prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]
            
            if not prot_abun.empty:
                # ... [Keep logic] ...
                
                fig_abun, ax_abun = plt.subplots(figsize=(6, 3))
                x_pos = np.arange(3)
                ax_abun.bar(x_pos, y_vals, fill=False, edgecolor="black", width=0.6)
                ax_abun.axhline(0, color='grey', linewidth=0.8)
                ax_abun.set_xticks(x_pos)
                ax_abun.set_xticklabels(labels, fontname='Arial') # Force Arial
                
                # Force Arial on Y Label and Ticks
                ax_abun.set_ylabel(r"Log$_2$ Fold Change", fontname='Arial')
                for tick in ax_abun.get_yticklabels():
                    tick.set_fontname('Arial')
                    
                st.pyplot(fig_abun)
            else:
                st.write("No abundance data available.")

            # ... [Rest of the code] ...