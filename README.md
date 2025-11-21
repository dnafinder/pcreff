[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/pcreff)

🌐 Overview
pcreff.m computes the PCR efficiency from RT-PCR calibration data, using a linear regression between log10-transformed input amounts and crossing point (Ct) values. It relies on the external function myregr.m to perform a rich least-squares regression and then derives the amplification efficiency and its uncertainty from the fitted slope. The function also reports a simple quality index indicating whether the calibration data are suitable for reliable efficiency estimation.

⭐ Features
- Uses log10(NG) versus CP regression to estimate PCR efficiency.
- Delegates regression to myregr.m, preserving its detailed statistical outputs.
- Computes the efficiency as E = 10^(-1/slope) - 1 (0–1, corresponding to 0–100%).
- Provides a standard error and 95% confidence interval for the efficiency.
- Computes a quality index based on regression statistics.
- Classifies the calibration as “good” or “not a good” calibrator based on the quality index.
- Prints a compact summary table with efficiency, error, confidence interval, quality and comment.

🛠️ Installation
1. Clone or download the repository:
   https://github.com/dnafinder/pcreff
2. Add the folder containing pcreff.m to your MATLAB path.
3. Install and add myregr.m and its related files from:
   https://github.com/dnafinder/myregression
4. Optionally, open and run pcreff.m directly in MATLAB Online using the badge at the top of this file.

▶️ Usage
Basic example:
    ng = [0.12 0.5 3 15 30];
    cp = [29.51 26.65 23.89 20.77 19.12; ...
          29.59 26.54 23.90 20.78 20.23; ...
          29.68 26.68 23.77 20.86 20.07];
    E = pcreff(ng, cp);

With non-verbose mode:
    E = pcreff(ng, cp, 0);

When VERBOSE = 1 (default), the underlying myregr.m function will display its own regression summary and plots, in addition to the efficiency summary printed by pcreff.m.

🔣 Inputs
- NG
  Type: numeric vector (row or column).
  Description: input amounts used for the RT-PCR calibration (e.g., ng of mRNA or cDNA). All values must be positive and finite.

- CP
  Type: numeric array (2D).
  Description: crossing point (Ct) data. CP must have size M×N, where N = numel(NG). Each column corresponds to a given NG value; rows can represent replicate measurements of the same concentration. The full CP matrix is passed to myregr.m, which handles repeated measures.

- VERBOSE (optional)
  Type: scalar (0 or 1).
  Description: controls whether myregr.m prints its detailed summary and plots.
  Default: 1 (verbose mode enabled).

📤 Outputs
The function returns a struct E with the following fields when called with an output argument:

- E.value
  Estimated PCR efficiency in fractional units (0–1), corresponding to 0–100% amplification per cycle.

- E.err
  Standard error of the efficiency estimate, propagated from the standard error of the regression slope provided by myregr.m.

In addition, pcreff.m prints a table with:
- Efficiency (in percent)
- Standard error
- 95% confidence interval for the efficiency
- Quality index
- Textual comment on calibrator suitability

📘 Interpretation
- Efficiency close to 1 (or 100%) indicates near-ideal doubling of product per cycle.
- Efficiency significantly below 1 indicates suboptimal amplification.
- The 95% confidence interval helps assess the precision and uncertainty of the efficiency estimate.
- The quality index is derived from regression statistics (STAT fields returned by myregr.m). A quality value below 0.1 is interpreted as a good calibrator; a value greater than or equal to 0.1 suggests that the calibration data may not be adequate for reliable efficiency estimation.

📝 Notes
- pcreff.m does not change the behavior of myregr.m. It simply calls myregr.m on log10-transformed NG versus CP and uses the returned slope and statistics.
- The function assumes a simple log-linear relation between log10(NG) and CP over the experimental range.
- Care should be taken to ensure that the selected NG values and CP measurements are within the linear dynamic range of the assay.
- The refactoring in 2025 preserves the original mathematical core while improving robustness, input checking, and documentation.

📚 Citation
If you use this function in scientific publications or reports, please cite:

Cardillo G. (2008) PCREfficiency: set the Efficiency of a RT-PCR to use in the relative quantification of transcripts.  
GitHub repository: https://github.com/dnafinder/pcreff

👤 Author
Author : Giuseppe Cardillo  
Email  : giuseppe.cardillo.75@gmail.com  
GitHub : https://github.com/dnafinder  

⚖️ License
This code is distributed under the MIT License. You are free to use, modify, and redistribute it, provided that proper credit is given to the original author and the source repository.
