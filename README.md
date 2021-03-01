# cellphonedb_shiny
Shiny app for visualisation of CellPhoneDB data

CellPhoneDB documentation can be found [here](https://github.com/Teichlab/cellphonedb)

Requires CellPhoneDB output files

Available at: https://saezlab.shinyapps.io/cellphonedb_shiny

## Instructions

Simple to use shiny application for cell-cell interaction analysis on scRNAseq data

1. Run CellPhoneDB using the instructions at the in the [CellphoneDB documentation](https://github.com/Teichlab/cellphonedb).
2. Upload the significant_meants.txt file and pvalues.txt file from the output of the above run by clicking browse or drag and drop.
3. Upload the metadata file you supplied to the CellPhoneDB command by clicking browse or drag and drop.
4. When all three files above show "Upload Complete" in the load bar, hit the "Generate report" button.
5. If you want to download any plots, select them from the dropdown menu and hit "Download Plot".
6. To view a specific feature in the Sankey plot, enter the feature name in the search box and hit the "Search" button.
