%OutputResults.m

%DESCRIPTION:
%This script takes the results from Lemoine's 2016 DICE model in MATLAB
%and outputs them into an Excel file within an output folder. The file name
%and location are specified in the script. This script can be edited to
%modify the list of outputs and reformat the layout of the table.

%INPUTS:
%Note that this code is a subsidary script that should be run within the 
%main_dice2016r script. This function relies on the following inputs
%computed within the DICE model:
%   C               = Consumption over time
%   Params          = A struct containing parameters relevant to the model
%   T =             = Atmospheric Temperature Increase over time
%   savingsrate     = The Savings Rate as determined by optimization
%   Ygross          = Gross Output
%   Ynet            = Net Output
%   k               = Capital
%   inv             = Investment
%   abatecost       = Abatement
%   abaterate       = Emissions Control Rate
%   emsind          = Industrial Emissions
%   iR              = Interest Rate
%   otherems        = Land-Based Emissions
%   emtax_pertCO2   = Carbon Tax
%   Carbon_ppm      = Atmospheric Concentration of Carbon
%                   = Social Cost of Carbon

%OUTPUTS:
%This function does not pass any values on to other scripts because it is
%the last script to run in the entire model. Instead, it outputs the values
%above in a Excel file named "BackEndOutput.xlsx" in a subfolder called
%"output" unless otherwise specified.

%AUTHOR:
%This Script was created by Benjamin Healy and Professor Gilbert Metcalf
%to run in conjuntion with Professor Derek Lemoine's 2016 DICE Model in MATLAB

%LAST EDITED: 2/24/21 by DEREK LEMOINE

%% Gathering Important Outputs

%Consumption per capita
cPc = C*1e12 ./ Params.pop*1e6;
cPc = cPc/1e3;   %thousand dollars per capita

%Climate damages fraction output
damFrac = Params.damcoeff * T.^2;             

%Damages
damages = (Ygross * Params.damcoeff) .* T.^2;

%Investment
inv = savingsrate.*( Ynet - abatecost );  

%Convert exogenous emissions from GtC to GtCO2
otheremsCO2 = otherems * Params.co2_per_c;

%Convert emissions by period to annual averages over the period
emsind_annual = emsind / Params.timestep;
emsind_annualCo2 = emsind_annual * Params.co2_per_c;

%Revert sign-convention for welfare and convert to annual summation
Welfare = -Params.timestep * Welfare;

%Save the Optimization Selection
switch Params.opt
    case 1
        abateDes = 'Optimal';
    case 0
        abateDes = ['Fixed at '  num2str(fixedAbate)];
    case 2
        abateDes = ['Linearly increasing to meet net zero carbon goal by '  num2str(targDate)];
    case 3
        abateDes = ['Avoiding Global Temperature Increases Over '  num2str(maxT) ' deg C'];
end 

%Vector of period increments
periodVec = 1:t;                               


dt = Params.timestep;
numyears = Params.numyears;

%Vector of year increments
yearVec = 2015  : dt : (2015 + numyears - dt); 

%Interest Rate
iR = zeros(t,1);
for i = 1:(t-1)
    iR(i) = (1 + Params.discrate) .* (((C(i+1)/Params.pop(i+1)) / ...
        (C(i)/Params.pop(i))))^ (Params.rra/Params.timestep) - 1;
end

%% Building an Output Table and Exporting to Excel

%Create an empty array to clear previous data
clearedCells = cell(35,510);

%Building data titles for each row
table2 = {'Period'; 'Year'; 'Industrial Emissions (GTCO2 per year)'; ...
    'Atmospheric Concentration of Carbon (ppm)'; 'Atmospheric Temperature Increase (deg C, rel to 1900)';...
    'Gross Output (trillion 2010$)'; 'Net Output (trillion 2010$)'; 'Capital (trillion 2010$)'; 'Labor (Population, millions)'; ...
    'Investment (trillion 2010$)'; 'Savings Rate'; 'TFP'; 'Consumption (trillion 2010$)'; ...
    'Consumption per Capita (thousand 2010$)'; 'Climate Damages Fraction Output'; ...
    'Damages (trillion 2010$)'; 'Emissions Control Rate'; 'Carbon Price (2010$ per tCO2)'; ...
    'Land Emissions (GtCO2 per year'; 'Abatement Cost (trillion 2010$)'; 'Emissions Intensity (GtCO2 per trillion 2010$)'; ...
    'Social Cost of Carbon (2010$ per tCO2)'; 'Annual Interest Rate'};

%Building the table of outputs
table1 = [table2, {periodVec; yearVec; emsind_annualCo2; Carbon_ppm; T; Ygross; Ynet; K; ...
    Params.pop; inv; savingsrate; Params.tfp; C; cPc; damFrac; damages; ...
    abaterate; emtax_pertCO2; otheremsCO2; abatecost; Params.sigma*Params.co2_per_c; SCC; iR}];

%Building a table of user-defined parameters for the model run
modelDetails = {'Date', date; 'Climate Model', Params.climatemodel; ...
    'Carbon Model', Params.carbonmodel; 'Damage Model', ...
    Params.damagemodel; 'Welfare', Welfare; 'Abatement Option', abateDes};

%Declaring the output filename
cd(Params.savedir);
filename = ['Output_' Params.filenaming '.xlsx'];

%Clearing previous data
writecell(clearedCells,filename,'Sheet',1,'Range','A1');

%Filing in Parameter Information
writecell(modelDetails,filename,'Sheet',1,'Range','A1');

%Exporting to the Excel file
writecell(table1,filename,'Sheet',1,'Range','A8');

%Closing and Saving Command Window Outputs to .txt File
fileRunOutput = [maindir, filesep, 'FileRunOutput.txt'];

%Check that File exists.
if isfile(fileRunOutput)
    diary off
    try %in case it fails because e.g. multiple processes
        movefile(fileRunOutput, Params.savedir);
    end
end

%Move to output folder
cd([maindir filesep 'output'])

%Save outputs for next guess in output directory
dtOut = Params.timestep;
nyOut = Params.numyears;
save([Params.carbonmodel '_guess.mat'], 'output_for_guess', ...
    'dtOut', 'nyOut');

%Return directory to Main Folder
cd(maindir);