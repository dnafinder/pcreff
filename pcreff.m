function E=pcreff(ng,cp,varargin)
%PCREFF Set the Efficiency of a RT-PCR to use in the relative
%quantification of transcripts.
%
%Reverse transcription(RT) followed by PCR is a powerful tool for the
%detection and quantification of mRNA. It is the most sensitive method for
%the detection and quantification of gene expression levels, in particular
%for low abundance mRNA.
%The relative quantification is based on the expression ratio of a target
%gene versus a reference gene. Some mathematical models have already been
%developed to calculate the relative expression ratios, with or without
%efficiency correction. Normally the PCR efficiency is set at 2 (the max 
%possible value) for the reference and target gene, but a difference in 
%PCR efficiency of 0.03 between the target and reference gene, the falsely
%calculated difference in expression ratio is 46% in case of Et<Er and 209%
%in the case of Et>Er. The difference will increase dramatically by higher
%efficiency differences: i.e. DE=0.05 (27% and 338%) and DE=0.1 (7.2% and
%1083%)
%This function computes the efficiency of PCR reaction and is based on MYREGR 
%function. If it is not present on the computer, pcreff will try to download
%it from FEX
%
% Syntax: 	pcreff(ng,cp,verbose)
%      
%     Inputs:
%           NG - Array of the ng of mRNA or cDNA used 
%           CP - Crossing points data. These data can be inserted as an
%           array or as a matrix (2 or more replicates for each ng data).
%           In the last case the mean will be calculated
%           verbose - Flag to display all information (default=1)
%     Outputs:
%           - Summary of MYREGR function
%           - PCR efficiency and error.
% Example:
%       ng=[0.1200 0.5000 3.0000 15.0000 30.0000];
%       cp=[29.5100 26.6500 23.8900 20.7700 19.1200; ...
% 29.5900 26.5400 23.9000 20.7800 20.2300; ...
% 29.6800 26.6800 23.7700 20.8600 20.0700];
%
%   Calling on Matlab the function: 
%             pcreff(ng,cp)
%
%   Answer is:
%
% (...) All the outputs of MYREGR function + calibration plot
%                                          PCR Efficency                                         
% ------------------------------------------------------------------------------------------------
%     Value     Percent       SE        Conf_Interval      Quality               Comment          
%     ______    _______    ________    ________________    ________    ___________________________
% 
%     1.7649    76.494     0.014658    1.7362    1.7937    0.018221    'This is a good calibrator'
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) PCREfficiency: set the Efficiency of a RT-PCR to use
% in the relative quantification of transcripts.
% http://www.mathworks.com/matlabcentral/fileexchange/20887

%Input errors handling
p = inputParser;
addRequired(p,'ng',@(x) validateattributes(x,{'numeric'},{'row','real','finite','nonnan','nonempty'}));
addRequired(p,'cp',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'verbose',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,ng,cp,varargin{:});
verbose=p.Results.verbose;
clear p

if exist('myregr.m','file')==0
    filename=unzip('https://it.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/15473/versions/10/download/zip','prova');
    Index = contains(filename,'myregr.m');
    current=cd;
    copyfile(filename{Index},current)
    rmdir('prova','s')
    clear filename Index current 
end

[slope,~,stat]=myregr(log10(ng),cp,verbose);
A.value=10^(-1/slope.value); %PCR efficiency (1<=A<=2)
A.err=1/slope.value^2*10^(-1/slope.value)*slope.se; %Error propagation
ci=A.value+[-1.96 +1.96].*A.err;
quality=((stat.cv*stat.rse/slope.value)^2)/stat.sse;
if quality>=0.1
    txt='This is not a good calibrator';
else
    txt='This is a good calibrator';
end
disp([blanks(41) 'PCR Efficency' blanks(41)])
disp(repmat('-',1,96));
disp(cell2table({A.value,(A.value-1)*100,A.err,ci,quality,txt},'VariableNames',{'Value','Percent','SE','Conf_Interval','Quality','Comment'}))
if nargout
    E=A;
end
end
