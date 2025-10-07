#!/bin/bash

# Check if arguments were provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 wat_site1 wat_site2 wat_site3 ..."
    exit 1
fi

/gist-post-processing/gistpp -i gist-Esw-dens.dx -op multconst -opt const 0.5 -o half-gist-Esw-dens.dx
/gist-post-processing/gistpp -i gist-Eww-dens.dx -i2 half-gist-Esw-dens.dx -op add -o gist-Etot-dens.dx
/gist-post-processing/gistpp -i gist-Etot-dens.dx -op multconst -opt const 0.125 -o gist-Etot.dx

/gist-post-processing/gistpp -i gist-dTStrans-dens.dx -i2 gist-dTSorient-dens.dx -op add -o gist-dTStot-dens.dx
/gist-post-processing/gistpp -i gist-dTStot-dens.dx -op multconst -opt const 0.125 -o gist-dTStot.dx
/gist-post-processing/gistpp -i gist-Etot.dx -i2 gist-dTStot.dx -op sub -o gist-dA.dx


# Loop through all provided wat_site arguments
for wat_site in "$@"; do
    echo "Processing wat_site: $wat_site"
    
    # Run the commands with the current wat_site variable
    /gist-post-processing/gistpp -i gist-gO.dx -i2 wat_site_${wat_site}.pdb -op defbp -opt const 2.5 -o bp_wat_${wat_site}.dx

    /gist-post-processing/gistpp -i bp_wat_${wat_site}.dx -i2 gist-Etot.dx -op mult -o wat_Etot_${wat_site}.dx
    /gist-post-processing/gistpp -i bp_wat_${wat_site}.dx -i2 gist-dTStot.dx -op mult -o wat_dTStot_${wat_site}.dx
    /gist-post-processing/gistpp -i bp_wat_${wat_site}.dx -i2 gist-dA.dx -op mult -o wat_dA_${wat_site}.dx

    /gist-post-processing/gistpp -i wat_Etot_${wat_site}.dx -op sum
    /gist-post-processing/gistpp -i wat_dTStot_${wat_site}.dx -op sum
    /gist-post-processing/gistpp -i wat_dA_${wat_site}.dx -op sum

    echo "Completed processing wat_site: $wat_site"
    echo "-------------------------------------"
done



