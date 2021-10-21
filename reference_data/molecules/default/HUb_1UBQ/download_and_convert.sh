# Clone the NWX_TA repo
git clone git@github.com:NWChemEx-Project/NWX_TA.git nwx_ta

# Convert the files
obabel ./nwx_ta/HUb_1UBQ/struct_h_added/*.pdb \
       -oxyz \
       -Oubiquitin_.xyz \
       -m

# Remove structure not based from [Ruger:2015]
rm ubiquitin_step1_pdbreader.xyz

# Remove downloaded repo
rm -rf nwx_ta