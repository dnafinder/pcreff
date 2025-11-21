function E = pcreff(ng, cp, varargin)
%PCREFF Set the efficiency of an RT-PCR for relative quantification.
%
%   E = PCREFF(NG, CP) estimates the PCR efficiency from calibration data
%   of input amount NG (e.g., ng of mRNA or cDNA) versus crossing points
%   CP (Ct values). The regression is performed by the external function
%   MYREGR and the efficiency is derived from the fitted slope.
%
%   E = PCREFF(NG, CP, VERBOSE) controls the verbosity of MYREGR:
%       VERBOSE = 0  -> no detailed output from MYREGR
%       VERBOSE = 1  -> MYREGR prints its summary (default)
%
%   NG must be a numeric vector (row or column) containing the input
%   amounts used in the calibration. CP must be a numeric array with the
%   same number of columns as the number of elements in NG. When CP is a
%   matrix, each column corresponds to a given NG value and rows represent
%   replicate crossing points for that concentration. CP is passed as-is
%   to MYREGR, which will handle repeated measures according to its own
%   logic (typically using the mean).
%
%   The efficiency E is computed from the slope of the regression line
%   between log10(NG) and CP as:
%
%       E.value = 10^(-1 / slope.value) - 1
%
%   so that 0 <= E.value <= 1 corresponds to 0–100% amplification
%   efficiency per cycle.
%
%   ------------------------------------------------------------------
%   Syntax:
%       E = pcreff(NG, CP)
%       E = pcreff(NG, CP, VERBOSE)
%
%   Inputs:
%       NG      - Numeric vector (row or column) of input amounts
%                 (e.g., ng of mRNA or cDNA).
%
%       CP      - Numeric array of crossing-point (Ct) data. Must have
%                 size MxN, where N = numel(NG). Each column corresponds
%                 to a given NG, and rows may represent replicates.
%
%       VERBOSE - (optional) scalar flag:
%                   0 -> no detailed output from MYREGR
%                   1 -> MYREGR prints its own summary and plots
%                 Default: 1
%
%   Output:
%       E       - Structure with the following fields:
%                   E.value  : estimated PCR efficiency (0–1)
%                   E.err    : standard error of the efficiency
%
%                 The function also prints a summary table with:
%                   - efficiency (in %)
%                   - standard error
%                   - 95% confidence interval
%                   - quality index
%                   - comment about calibrator suitability
%
%   ------------------------------------------------------------------
%   Notes:
%   ------------------------------------------------------------------
%   - This function relies on MYREGR for linear regression and for
%     obtaining slope and related statistics. MYREGR must be available
%     on the MATLAB path. The current implementation expects the
%     repository:
%
%         https://github.com/dnafinder/myregression
%
%   - The quality index is computed from MYREGR outputs as:
%
%         quality = ((STAT.cv * STAT.rse / slope.value)^2) / STAT.sse
%
%     A value >= 0.1 is interpreted as "not a good calibrator".
%
%   ------------------------------------------------------------------
%   Metadata:
%       Author : Giuseppe Cardillo
%       Email  : giuseppe.cardillo.75@gmail.com
%       GitHub : https://github.com/dnafinder
%       Created: 2008-01-01
%       Updated: 2025-11-21
%       Version: 2.0.0
%
%   Citation:
%       Cardillo G. (2008) PCREfficiency: set the Efficiency of a RT-PCR
%       to use in the relative quantification of transcripts.
%       GitHub: https://github.com/dnafinder/pcreff
%
%   License:
%       This code is distributed under the MIT License.
%   ------------------------------------------------------------------

%% Input validation
p = inputParser;
p.FunctionName = 'pcreff';

addRequired(p, 'ng', @(x) validateattributes( ...
    x, {'numeric'}, {'vector','real','finite','nonnan','nonempty','>',0}));

addRequired(p, 'cp', @(x) validateattributes( ...
    x, {'numeric'}, {'2d','real','finite','nonnan','nonempty'}));

addOptional(p, 'verbose', 1, @(x) ...
    isnumeric(x) && isscalar(x) && ismember(x, [0 1]));

parse(p, ng, cp, varargin{:});
ng      = p.Results.ng;
cp      = p.Results.cp;
verbose = p.Results.verbose;
clear p

% Ensure NG is a row vector for consistency
ng = ng(:).';  % 1 x N

% Check consistency between NG and CP
n = numel(ng);
if size(cp, 2) ~= n
    error('pcreff:SizeMismatch', ...
        'Number of columns in CP (%d) must match numel(NG) (%d).', ...
        size(cp, 2), n);
end

%% Check dependency: MYREGR
if exist('myregr', 'file') ~= 2
    error('pcreff:MissingDependency', ...
        ['The function MYREGR is required but was not found on the path.\n' ...
         'Please install it from: https://github.com/dnafinder/myregression']);
end

%% Call MYREGR on log10-transformed NG
[slope, ~, stat] = myregr(log10(ng), cp, verbose);

%% Compute efficiency and its standard error
E.value = 10^(-1 / slope.value) - 1;  % (0 <= E.value <= 1)
E.err   = slope.se * (1 + E.value) * log(10) / (slope.value^2);

% 95% confidence interval for efficiency (using n-2 degrees of freedom)
ci = E.value + [-1 1] .* tinv(0.975, n - 2) .* E.err;

% Quality index based on MYREGR statistics
quality = ((stat.cv * stat.rse / slope.value) ^ 2) / stat.sse;

if quality >= 0.1
    txt = {'This is not a good calibrator'};
else
    txt = {'This is a good calibrator'};
end

%% Display summary table
disp([blanks(41) 'PCR Efficiency' blanks(41)])
disp(repmat('-', 1, 96));
disp(table(E.value * 100, E.err, ci, quality, txt, ...
    'VariableNames', {'Efficiency','SE','Conf_Interval','Quality','Comment'}));

end
