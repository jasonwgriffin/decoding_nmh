# **Code Directory for the following:**

## **Title:** Decoding the temporal dynamics of face-specific neural representations in autism
**Authors:** Jason W. Griffin, Alan H. Gerber, Kaylee Litson, Susan Faja, Shafali Jeste, Natalia Kleinhans, Geraldine Dawson, Adam Naples, April R. Levin, Sara Jane Webb, Frederick Shic, Catherine Sugar, James Dziura, and James C. McPartland for the Autism Biomarkers Consortium for Clinical Trials

**Journal:** Nature Mental Health

 ## Decription of Code Directory

 1. MVPA/ - This is matlab code related to running the MVPA analyses. The code reads in data, runs the mvpa analyses, and saves the output for the 4 conditions: image-selectivity, face-selectivity, identity-selectivity, and orientation-selectivity.
 2. R/ - This is the directory for reading in the MVPA output into R and conducting statistical analyses and data visualization.
 3. R/decoding_analysis.Rmd - This is the top level markdown file that runs analyses and scripts for the src folder
 
 4. R/src/Figure1.R - This is source code to generate Figure 1.
 5. R/src/Figure2.R - This is source code to generate Figure 2.
 6. R/src/Figure3.R - This is source code to generate Figure 3.
 7. R/src/Figure4.R - This is source code to generate Figure 4.
 8. R/src/ED_Figure2.R - This is source code to generate Extended Data Figure 2.
 9. R/src/Figure4.R - This is source code to generate Figure 4.
 10. R/src/Figure4.R - This is source code to generate Figure 4.
 11. R/src/make_table1.R - This is source code to generate Table 1.
 12. R/src/make_table2.R - This is source code to generate Table 2.
 13. R/src/process_output.R - This is source code for reading in the MVPA output
 14. R/src/CNR_analysis.R - This is source code for conducting sensitivity analyses, including the derivation of the Contrast-to-noise (CNR) ratio scores. This also produces Extended Data Table 1.


## How to use

The matlab files in the MVPA directory will run the decoding/MVPA analyses. The primary input data for this is acquired from NDAR (Collection ID #2288). This code will derive and save output that will be the input data for the top level R Markdown file `decoding_analysis.Rmd`. `decoding_analysis.Rmd.` will generate all of the figures and tables of the paper, which contain the primary results. 

## Questions about the code?
You can contact jasongriffin@uh.edu for any questions
