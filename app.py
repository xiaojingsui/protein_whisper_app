# ... (Keep the sidebar code as is) ...

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
            # 2. Capitalize Gene Name
            gene = protein["gene"].upper() 
            
            # 3. Link UniProt ID to database
            uniprot_link = f"[{uniprot}](https://www.uniprot.org/uniprot/{uniprot})"
            st.markdown(f"### Protein: **{gene}** ({uniprot_link})")
            
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
                 # 1. SWAP: Interactive Message moved to TOP
                 st.info("💡 **Interactive:** Hover over a dot on the volcano plot to see it on the structure. Grey dots are non-significant.")

            # --- Prepare Data for Visualization ---
            viz_data = []
            full_color = "#E1341E"
            half_color = "#1ECBE1"
            grey_color = "#D3D3D3" # Light grey for non-sig

            if not peps.empty:
                # Normalization for heatmap 
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
                html_view = render_interactive_dashboard(
                    structure_text, file_format, viz_df, fc_cutoff, p_cutoff, selected_condition
                )
                components.html(html_view, height=650)
                
                # Legend for Colors
                if color_mode == "Fold-change heatmap" and sig_count > 0:
                    c_left, c_rest = st.columns([1, 2])

                    with c_left:
                          fig_cb, ax_cb = plt.subplots(figsize=(1, 0.1))
                          cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb, orientation="horizontal")
                          cb1.set_label(r"Log$_2$FC", fontsize=8)
                          cb1.ax.tick_params(labelsize=8)
                          st.pyplot(fig_cb, use_container_width=False)
                
                # 1. SWAP: Significant Count Message moved to BOTTOM
                if sig_count > 0:
                      st.success(f"{sig_count} significant peptides found out of {len(peps)} total.")
                else:
                      st.info(f"No significant peptides found (Total: {len(peps)}). Adjust filters to see changes.")

            # --- Abundance Plot ---
            st.subheader("Protein Abundance")
            prot_abun = abun_df[abun_df["uniprot_id"] == uniprot]
            
            if not prot_abun.empty:
                labels = ["Soluble", "Pellet", "Total"]
                suffix_map = {"Soluble": "soluble", "Pellet": "pellet", "Total": "total"}
                y_vals = []
                for label in labels:
                    suffix = suffix_map[label]
                    ac = f"AvgLog₂({selected_condition}).{suffix}"
                    y_vals.append(float(prot_abun[ac].values[0]) if ac in prot_abun.columns else 0)
                
                fig_abun, ax_abun = plt.subplots(figsize=(6, 3))
                x_pos = np.arange(3)
                ax_abun.bar(x_pos, y_vals, fill=False, edgecolor="black", width=0.6)
                ax_abun.axhline(0, color='grey', linewidth=0.8)
                ax_abun.set_xticks(x_pos)
                ax_abun.set_xticklabels(labels)
                ax_abun.set_ylabel(r"Log$_2$ Fold Change")
                st.pyplot(fig_abun)
            else:
                st.write("No abundance data available.")

            st.markdown("---")
            st.subheader("Peptide Data")
            if not peps.empty:
                st.dataframe(peps[["Peptide sequence", "Start position", "End position", avg_col, pval_col, "is_sig"]])