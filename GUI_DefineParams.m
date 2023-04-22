%GUI_DefineParams_ToLemoine_Healy

%DESCRIPTION:
%This script provide a graphical user interface (GUI) which allows users to
%easily define certain parameters for the 2016 DICE Model. Drop-down menus 
%in the GUI allow users to define the models used for climate, carbon, and 
%damage estimates. Users may also decide whether to optimize over the study
%period using dropdown menu. Finally, users can define the time step and
%time horizon (the end of the study period) using interactive text boxes.

%NOTE: This script is built to run with an adapted version of Professor 
%Derek Lemoine's 2016 DICE Model in Matlab


%INPUTS: none


%OUTPUTS:
%   Params.timestep       = the time step for which the model will solve 
%                           and output results
%   Params.numyears       = the time horizon at which the model will stop   
%   Params.climatemodel   = the type of climate model to be used in solving
%                      1. if 'dice', the model will  use DICE-2016R 
%                         climate dynamics
%                      2. if 'update', the model will use Geoffroy et al 
%                         dynamics (see Dietz et al. 2020)
%                      3. if 'other', the user will provide a script, 
%                         'other_climate.m' to load in a different guess
%
%   Params.carbonmodel = the type of carbon model to be used in solving
%                      1. if 'dice', the model will  use DICE-2016R 
%                         climate dynamics
%                      2. if 'gol', the model will use Golosov et al. 
%                         (2014), via Lemoine (2020) and with starting 
%                         stocks from Dietz et al. (2020)
%                      3. if 'joos' or 'fair', the model will use Joos et 
%                         al. (2013), via Dietz et al. (2020)
%                      4. if 'other', the user will provide a script, 
%                         'other_carbon.m' to load in a different guess
%
%   Params.damagemodel = the type of damage model to be used in solving
%                      1. if 'dice', the model will  use DICE-2016R 
%                         climate dynamics
%                      2. if 'expert', the model will use Pindyck (2019) 
%                         damages as implemented by Lemoine (2021)
%                      3. if 'other', the user will provide a script, 
%                         'other_damage.m' to load in a different guess
%
%   Params.opt         = A boolean operator that toggles optimization on
%                        (= 1) or off (= 0) in the model's solution
%
%AUTHOR:
%This Script was created by Benjamin Healy and Professor Gilbert Metcalf
%to run in conjuntion with Professor Derek Lemoine's 2016 DICE Model in MATLAB

%Last edited: 2/21/21 by Ben Healy


%% Settting Up the Figure

%Creating the Figure for the GUI
diceGUI = figure('Position', [500, 100, 400, 500], 'Name', ...
    'Set User-Defined Inputs', 'numbertitle', 'off');

% Removing Unnecessary Features
set(diceGUI, 'MenuBar', 'none');
set(diceGUI, 'ToolBar', 'none');

%Creating feature spacing within the figure
boxPos = linspace(425, 50, 7);

%Definining Title Formatting for Input Boxes
titleL = 5/500;
titleW = 490/500;
titleH = 60/600;
titleColor = '#c1d9e3';
fSize = 12;

%% Defining Input Options and Defaults

%Defining Options for Each Pop-Up Menu
climateModels = {'dice', 'update', 'other'};
carbonModels = {'dice', 'gol', 'joos', 'fair', 'other'};
damageModels = {'dice', 'expert', 'other'};
optimizingChoices = {'Yes', 'No, I want to fix abatement', ...
    'No, I want to set a zero carbon target date', ...
    'Yes, and constrain maximum temperature increase'};

%Defining Default Values for Input Boxes
dtDefault = 10;
timeHorizonDefault = 500;

%%  Building GUI Features

%Climate Model Selection
annotation('textbox', [titleL, (boxPos(1) + 5)/500, titleW, titleH],...
    'String', 'Select a Climate Model to Use', 'HorizontalAlignment', 'center', ...
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', fSize)
climMod = uicontrol('Style', 'popupmenu', 'String', climateModels, ...
    'Position', [20, boxPos(1), 375, 30]);



%Carbon Model Selection
annotation('textbox', [titleL, (boxPos(2) + 5)/500, titleW, titleH],...
    'String', 'Select a Carbon Model to Use', 'HorizontalAlignment', 'center', ... 
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', fSize)
carbMod = uicontrol('Style', 'popupmenu', 'String', carbonModels, ...
    'Position', [20, boxPos(2), 375, 30]);



%Damage Model Selection
annotation('textbox', [titleL, (boxPos(3) + 5)/500, titleW, titleH],...
    'String', 'Select a Damage Model to Use', 'HorizontalAlignment', 'center', ...
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', fSize)
damMod = uicontrol('Style', 'popupmenu', 'String', damageModels, ...
    'Position', [20, boxPos(3), 375, 30]);


  
%Optimization Selection
annotation('textbox', [titleL, (boxPos(4) + 5)/500, titleW, titleH],...
    'String', 'Do you Want to Choose Abatement to Maximize Welfare?', ...
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', ...
    fSize, 'HorizontalAlignment', 'center')
optSel = uicontrol('Style', 'popupmenu', 'String', optimizingChoices, ...
    'Position', [20, boxPos(4), 375, 30]);



%Timestep Selection
annotation('textbox', [titleL, (boxPos(5) - 5)/500, titleW, titleH],...
    'String', 'Select a Timestep for Calculations Across the Study Period (years)',... 
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', fSize, ...
     'HorizontalAlignment', 'center')
dtSel = uicontrol('Style', 'Edit', 'String', dtDefault, 'Position', ...
    [20, boxPos(5), 360, 20]);



%Time Horizon Selection
annotation('textbox', [titleL, (boxPos(6) - 5)/500, titleW, titleH],...
    'String', 'Select the Time Horizon at Which the Model Should Finish (years)',... 
    'BackgroundColor', titleColor, 'FontWeight', 'bold', 'FontSize', fSize, ...
     'HorizontalAlignment', 'center')
horSel = uicontrol('Style', 'Edit', 'String', timeHorizonDefault, ...
    'Position', [20, boxPos(6), 360, 20]);



%Enter Button
entSel = uicontrol('Style', 'pushbutton', 'String', 'Enter', 'Position', ...
    [150, boxPos(7), 100, 20], 'Callback', 'uiresume(diceGUI)');
uiwait(diceGUI)

%% Save Selected Inputs

%Save Climate Model Selection
switch climMod.Value
    case 1
        Params.climatemodel = climateModels{1};
    case 2
        Params.climatemodel = climateModels{2};
    case 3
        Params.climatemodel = climateModels{3};
end

%Save Carbon Model Selection
switch carbMod.Value
    case 1
        Params.carbonmodel = carbonModels{1};
    case 2
        Params.carbonmodel = carbonModels{2};
    case 3
        Params.carbonmodel = carbonModels{3};
    case 4
        Params.carbonmodel = carbonModels{4};
    case 5
        Params.carbonmodel = carbonModels{5};
end


%Save the Damage Model Selection
switch damMod.Value
    case 1
        Params.damagemodel = damageModels{1};
    case 2
        Params.damagemodel = damageModels{2};
    case 3
        Params.damagemodel = damageModels{3};
end

%Save the Optimization Selection
switch optSel.Value
    case 1
        Params.opt = 1;
    case 2
        Params.opt = 0;
    case 3
        Params.opt = 2;
    case 4
        Params.opt = 3;
end

%Save the Timestep Selection
Params.timestep = str2double(dtSel.String);

%Save the Time horizon Selection
Params.numyears = str2double(horSel.String);

%Check that time horizon and timestep choices are logical
if (Params.timestep <= 0) || (Params.timestep > Params.numyears) || (Params.numyears <= 0)
    close(diceGUI)
    error('Please review your timestep and time horizon for errors')
end

%Display Termination Message and Close GUI
disp('User-Defined Inputs Were Successfully Saved')
close(diceGUI)




%% Creating a follow-up GUI for other damage model prediction by user
if strcmp(Params.damagemodel, 'other')
    
    %Defining title, input box prompts, and default inputs
    dlgtitle1 = 'Set User-Defined Guesses for "Other" Damage Model Selection';
    prompt1 = {'Enter the Damage Coefficient'};
    definput1 = {'0.00236'};
    
    %Creating input dialogue GUI and pausing code for user interaction
    answerDam = inputdlg(prompt1, dlgtitle1,[1, 100], definput1);
    waitfor(answerDam);
    
    %extracting results
    Params.damcoeff = str2double(answerDam{:});

    %Display Termination Message
    disp('User-Defined Inputs Were Successfully Saved for Damage Model')
end




%% Creating a follow-up GUI for other climate model prediction by user

if strcmp(Params.climatemodel, 'other')
    
    %Defining title, input box prompts, and default inputs
    dlgtitle2 = 'Set User-Defined Guesses for "Other" Climate Model Selection';
    prompt2 = {'Enter value for Equilibrium Climate Sensitivity (W/m^2)', ...
        'Enter value for Forcing per Degree Warming (W/m^2/K)', ...
        'Enter value for Warming Delay Parameter',... 
        'Enter value for Transfer Warming Coefficient from ocean to surface',...
        'Enter value for Transfer Warming Coefficient from surface to ocean'};
    definput2 = {'3.6813', '1.1875', '0.1005', '0.088', '0.025'};
    
    %Creating input dialogue GUI and pausing code for user interaction
    answerClim = inputdlg(prompt2, dlgtitle2,[1, 100],definput2);
    waitfor(answerClim);
    
    %extracting results
    Params.f2x = str2double(answerClim{1});
    Params.lambda = str2double(answerClim{2});
    Params.phi1 = str2double(answerClim{3});
    Params.phi3 = str2double(answerClim{4});
    Params.phi4 = str2double(answerClim{5});
    
    %Display Termination Message
    disp('User-Defined Inputs Were Successfully Saved for Climate Model')
end
 



%% Creating a follow-up GUI for other carbon model prediction by user

if strcmp(Params.carbonmodel, 'other')
    %Building the GUI
    carbModGUI = figure('Position', [500, 300, 500, 400], 'Name', ...
    'Set User-Defined "Other" Guesses for Carbon Dynamics', 'numbertitle', 'off');

    % Removing Unnecessary Features
    set(carbModGUI, 'MenuBar', 'none');
    set(carbModGUI, 'ToolBar', 'none');
    
    %Defining default answers
    mMatrix = {0.9722, 0.0455, -0.0000; 0.0279, 0.9529,...
        0.0003; -0.0001, 0.0015, 0.9997};
    emSinks =  {1, 0, 0};
    m0 = {851, 460, 1740};
    
    %Defining prompts for each input parameter
    prompt3 = {'Enter a 3 by 3 transition matrix for carbon sinks', ...
        'Enter a 1 by 3 matrix of initial emissions deposition', ...
        'Enter a 1 by 3 matrix of initial stocks (Gt C)'};
 

    %Creating input boxes for the transition matrix of carbon sinks
    for n = 1:3
        for m = 1:3
            transM(n,m) = uicontrol('Style', 'Edit', 'String', mMatrix{n,m}, ...
                'Position', [100 + (100 * (m-1)), 350 - (25 * n), 100, 20]);
        end
    end
    
    annotation('textbox', [95/500, 270/400, 310/500, 95/400],...
        'String', prompt3{1}, 'BackgroundColor', titleColor, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
        'FontSize', fSize)
    
    
    %Creating  input boxes for the initial  emissions deposition matrix
    for m2 = 1:3
        iEm(m2) = uicontrol('Style', 'Edit', 'String', emSinks{m2}, ...
            'Position', [100 + (100 * (m2-1)), 200, 100, 20]);
    end
    
    annotation('textbox', [95/500, 195/400, 310/500, 45/400],...
        'String', prompt3{2}, 'BackgroundColor', titleColor, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
        'FontSize', fSize)
    
    
    %Creating input boxes for the total stock matrix
    for m3 = 1:3
        tSt(m3) = uicontrol('Style', 'Edit', 'String', m0{m3}, ...
            'Position', [100 + (100 * (m3-1)), 125, 100, 20]);
    end

    annotation('textbox', [95/500, 120/400, 310/500, 45/400],...
        'String', prompt3{3}, 'BackgroundColor', titleColor, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center',...
        'FontSize', fSize)    
    
    
    %Enter Button
    entSel2 = uicontrol('Style', 'pushbutton', 'String', 'Enter', 'Position', ...
        [200, 50, 100, 20], 'Callback', 'uiresume(carbModGUI)');
    uiwait(carbModGUI)

    %Extract results for transition matrix of carbon sinks
    Params.Mmatrix = zeros(3);
    for n4 = 1:3
        for m4 = 1:3
            Params.Mmatrix(n4,m4) = str2double(transM(n4,m4).String);
        end
    end
    
    %extract results for initial emissions deposition
    Params.emsinks = [str2double(iEm(1).String), str2double(iEm(2).String), str2double(iEm(3).String)];
    
    %extract results for initial stocks of carbon
    Params.M0 = [str2double(tSt(1).String), str2double(tSt(2).String), str2double(tSt(3).String)];
    
    %Display Termination Message and Close GUI
    disp('User-Defined Inputs Were Successfully Saved for Carbon Model')
    close(carbModGUI)
  
end

%% Creating a follow-up GUI for choosing fixed abatement

if Params.opt == 0
        
    %Defining title, input box prompts, and default inputs
    dlgtitle2 = 'Fixing Abatement Over the Model Run';
    prompt2 = {'Choose an Abatement Rate (format inputs as decimals)'};
    definput2 = {'0.03'};
    
    %Creating input dialogue GUI and pausing code for user interaction
    fixedAbate = inputdlg(prompt2, dlgtitle2,[1, 100], definput2);
    waitfor(fixedAbate);
    
    %extracting results
    fixedAbate = str2double(fixedAbate);

    %Display Termination Message
    disp('User-Defined Inputs Were Successfully Saved for Fixed Abatement')
    
end


%% Creating a follow-up GUI for choosing a zero carbon target date

if Params.opt == 2
        
    %Defining title, input box prompts, and default inputs
    dlgtitle3 = 'Linearly Increasing Abatement Rate to Reach Net Zero Carbon Emissions';
    prompt3 = {'Set a Target Date for Net Zero Carbon Emissions'};
    definput3 = {'2050'};
    
    targDate = {'inf'};
    while (ceil((str2double(targDate{:}) - 2015) / Params.timestep) + 1) ...
            > (Params.numyears / Params.timestep)
        %Creating input dialogue GUI and pausing code for user interaction
        targDate = inputdlg(prompt3, dlgtitle3,[1, 110], definput3);
        waitfor(targDate);
        if (ceil((str2double(targDate{:}) - 2015) / Params.timestep) + 1) ...
                > (Params.numyears / Params.timestep)
            dlgtitle3 = 'ERROR: Target Date Must Occur Before Time Horizon';
            prompt3 = {['Set a Target Date for Net Zero Carbon Emissions (before '...
                num2str(2015 + (Params.numyears - Params.timestep)) ')']};
        end
    end
    
    
    %extracting results
    targDate = str2double(targDate);

    %Display Termination Message
    disp('User-Defined Inputs Were Successfully Saved for a Zero Carbon Target Date')

end

%% Creating a follow-up GUI for choosing a zero carbon target date

if Params.opt == 3
        
    %Defining title, input box prompts, and default inputs
    dlgtitle3 = 'Optimizing Abatement to Prevent a Defined Temperature Increase';
    prompt3 = {'Choose a Maximum Desired Increase in Global Temperatures (Degrees Celsius)'};
    definput3 = {'12'};
    
    %Creating input dialogue GUI and pausing code for user interaction
    maxT = inputdlg(prompt3, dlgtitle3,[1, 100], definput3);
    waitfor(maxT);
    
    %extracting results
    maxT = str2double(maxT{:});

    %Display Termination Message
    disp('User-Defined Inputs Were Successfully Saved for Maximum Temperature Increase')
    
end