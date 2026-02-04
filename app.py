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

    /* --- 5. TEXT-ONLY LINK BUTTONS --- */
    div[data-testid="stSidebar"] .stButton > button {
        background-color: transparent !important;
        border: none !important;
        box-shadow: none !important;
        color: #006064 !important; 
        text-decoration: underline !important;
        padding: 0px !important;
        width: auto !important;
        height: auto !important;
        min-height: 0px !important;
        margin-top: 5px !important; 
        line-height: 1.2 !important;
        font-weight: normal !important;
        display: inline-block !important;
    }

    div[data-testid="stSidebar"] .stButton > button * {
        font-size: 14px !important;     
    }

    div[data-testid="stSidebar"] .stButton > button:hover {
        color: #19CFE2 !important;     
        text-decoration: none !important; 
        background-color: transparent !important;
    }
    
    div[data-testid="stSidebar"] .stButton > button:focus,
    div[data-testid="stSidebar"] .stButton > button:active {
        border: none !important;
        box-shadow: none !important;
        outline: none !important;
        color: #006064 !important;
        background-color: transparent !important;
    }

    .example-label {
        margin-top: 5px !important; 
        font-size: 14px !important;        
        color: #666 !important;
        font-weight: 600 !important;
        margin-right: 5px !important;
        white-space: nowrap !important;
    }

    /* --- 6. BRAND LOGO (NEW) --- */
    .nav-logo {
        position: fixed !important;
        top: 14px !important;
        left: 30px !important;
        z-index: 1000002 !important; /* Sit on top of navbar */
        font-family: Arial, Helvetica, sans-serif !important;
        font-size: 24px !important;
        font-weight: bold !important;
        color: #006064 !important; /* Teal Brand Color */
        pointer-events: none; /* Allows clicks to pass through if necessary */
    }
    </style>
""", unsafe_allow_html=True)