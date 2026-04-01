import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

# Set page config for a cleaner look
st.set_page_config(page_title="Chiral Centers Finder", page_icon="🧪")

st.title("Chiral Centers Finder")
st.write("This app determines the chiral centers of a molecule from its SMILES string and displays its 2D structure.")

# Input for SMILES string
# Default to Artemether's SMILES from the original file
default_smiles = "CO[C@H]1C[C@@H]2[C@@H]3CCC(C)[C@H](O2)[C@@]4(OOC(C)C)[C@H]3CC[C@@H]14"
smiles = st.text_input("Enter SMILES String:", value=default_smiles)

if smiles:
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        st.error("Invalid SMILES string! Please enter a valid one.")
    else:
        # Layout columns for side-by-side display
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Molecule Structure")
            
            # Using tabs for 2D and 3D
            tab1, tab2 = st.tabs(["2D View", "3D View"])
            
            with tab1:
                try:
                    # Compute 2D coordinates for a better layout
                    from rdkit.Chem import AllChem
                    AllChem.Compute2DCoords(mol)
                    
                    # Generate an SVG of the molecule which is more reliable for web browsers
                    img_svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(300, 300), useSVG=True)
                    
                    # Display the SVG image using st.image or st.markdown
                    st.markdown(img_svg, unsafe_allow_html=True)
                except Exception as e:
                    st.warning(f"Could not generate 2D image: {e}")
                    
            with tab2:
                try:
                    from rdkit.Chem import AllChem
                    import py3Dmol
                    from stmol import showmol
                    
                    # Create a copy so we don't mess up 2D coordinates for the chiral calculation later
                    mol_3d = Chem.Mol(mol)
                    # Add hydrogens and generate 3D conformer
                    mol_3d = Chem.AddHs(mol_3d)
                    AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(mol_3d)
                    
                    # Convert to MolBlock for py3Dmol
                    molblock = Chem.MolToMolBlock(mol_3d)
                    
                    # Create viewer
                    viewer = py3Dmol.view(width=400, height=400)
                    viewer.addModel(molblock, 'mol')
                    viewer.setStyle({'stick': {}})
                    viewer.zoomTo()
                    
                    # Display in streamlit
                    showmol(viewer, height=400, width=400)
                except Exception as e:
                    st.warning(f"Could not generate 3D image: {e}")
        
        with col2:
            st.subheader("Chiral Centers Results")
            
            # Prepare molecule for stereochemistry
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Assign stereochemistry
            Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
            
            # Find chiral centers
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
            
            if not chiral_centers:
                st.info("No assigned chiral centers found in this molecule.")
            else:
                st.write(f"**Total number of chiral centers:** {len(chiral_centers)}")
                
                # Format chiral centers data into a list of dictionaries for the table
                centers_data = [{"Atom Index": idx, "Configuration": config} for idx, config in chiral_centers]
                
                # Display simply as a dataframe/table
                st.dataframe(centers_data, use_container_width=True)
