# pcreff
Set the Efficiency of a RT-PCR to use in the relative quantification of transcripts.<br/>
Reverse transcription(RT) followed by PCR is a powerful tool for the
detection and quantification of mRNA. It is the most sensitive method for
the detection and quantification of gene expression levels, in particular
for low abundance mRNA.
The relative quantification is based on the expression ratio of a target
gene versus a reference gene. Some mathematical models have already been
developed to calculate the relative expression ratios, with or without
efficiency correction. Normally the PCR efficiency is set at 2 (the max 
possible value) for the reference and target gene, but a difference in 
PCR efficiency of 0.03 between the target and reference gene, the falsely
calculated difference in expression ratio is 46% in case of Et<Er and 209%
in the case of Et>Er. The difference will increase dramatically by higher
efficiency differences: i.e. DE=0.05 (27% and 338%) and DE=0.1 (7.2% and
1083%)
This function computes the efficiency of PCR reaction and is based on MYREGR 
function. If it is not present on the computer, pcreff will try to download
it from FEX

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008) PCREfficiency: set the Efficiency of a RT-PCR to use
in the relative quantification of transcripts.
http://www.mathworks.com/matlabcentral/fileexchange/20887
