# The skin microbiome of patients with atopic dermatitis normalizes gradually during treatment

#### Veda D. Khadka*, Felix Michael Key*, Carolina Romo-González, Adrián Martínez-Gayosso, Blanca L Campos-Cabrera, Armando Gerónimo-Gallegos, Tucker C Lynn, Carola Durán-Mckinster, Rafael Coria Jimenez, Tami D Lieberman**, and Maria Teresa García-Romero**
 '*' First authors <br>
'**' Last authors

Original Research, Front. Cell. Infect. Microbiol. - Microbiome in Health and Disease

Background: Atopic dermatitis (AD) is characterized by an altered skin microbiome dominantly colonised by S. aureus. Standard treatment includes emollients, anti-inflammatory medications and antiseptics.
Objectives: To characterize changes in the skin microbiome during treatment for AD.

Methods: The skin microbiomes of children with moderate-to-severe AD and healthy children were investigated in a longitudinal prospective study. Patients with AD were randomized to receive either standard treatment with emollients and topical corticosteroids or standard treatment with the addition of dilute bleach baths (DBB) and sampled at four visits over a three-month period. At each visit, severity of AD was measured, swabs were taken from four body sites and the composition of the microbiome at those sites was assessed using 16S rRNA amplification.

Results: We included 14 healthy controls and 28 patients. We found high relative abundances of S. aureus in patients, which correlated with AD severity and reduced apparent alpha diversity. As disease severity improved with treatment, the abundance of S. aureus decreased, gradually becoming more similar to the microbiomes of healthy controls. After treatment, patients who received DBB had a significantly lower abundance of S. aureus than those who received only standard treatment.

Conclusions: There are clear differences in the skin microbiome of healthy controls and AD patients that diminish with treatment. After three months, the addition of DBB to standard treatment had significantly decreased the S. aureus burden, supporting its use as a therapeutic option. Further study in double-blinded trials is needed

 <a href = "https://www.ncbi.nlm.nih.gov/bioproject/759575"> Sequence data can be found here </a>

## Included: 

<ul>
  <li> <strong> Qiime2_Analysis </strong> </li>
  <ul>
  <li> Qiime2 bash script named "qiime2_dataprocessing.sh". Run this script on raw sequencing data to get the ball rolling. </li>
  <li> Metadata table required by <a href = "https://qiime2.org/"> Qiime2 </a> named "metatable_sra_16s_public.tsv"  </li>
  <li> Classifier trained on <a href = "https://www.arb-silva.de/" > Sliva </a> ASVs, required by Qiime2 script named "classifier.qza.gz" </li>
  </ul> <br>
 <li> <strong> R_Analysis </strong> </li> 
 R workflow for analysing Qiime2 output (level-7) <br>
  <ul>
   <li> Functions for running R scripts, named "16StaphAD_functions.R" </li>
   <li> Code for preprocessing Qiime2 output, named "Preprocessing_QC.Rmd". Run this script on the Qiime2 output generated above. </li>
   <li> Code for generating figures from the paper, named "Figure_Generation.Rmd" </li> <br>
 <li> <strong> Data files </strong> </li>
Qiime2 output necessary for running R preprocessing and figure generation script:
  <ul>
   <li> "level-7.csv" </li>
   <li> Other .csv files in this folder are dataframes generated through the "Preprocessing_QC.Rmd" script. You may use these dataframes in lieu of creating them yourself in the Figure_Generation.Rmd script. </li>
   </ul> 
  </ul> 
  
</ul>
