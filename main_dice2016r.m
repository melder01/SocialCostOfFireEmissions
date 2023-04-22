% Replicates DICE-2016R

%Original DICE-2016R code by Derek Lemoine (10/4/2020)
%Edited 2/21/21 by Ben Healy for Gilbert Metcalf
%Last edited 3/20/23 by Molly Elder

%%

clear all;
close all;

global iteration_number

%Exporting Command Window Outputs to .txt file
diary 'FileRunOutput.txt'


%% User Options %%%%%%%%%%%%%%%

customdir = ''; % prefix for save directory; not used if running on HPC

%These parameters are now defined in the GUI script
run GUI_DefineParams

%{
Params.timestep = 5; % timestep, in years; default: 5
Params.numyears = 500; % number of years in horizon; default: 500

Params.climatemodel = 'dice'; % 'dice' or 'update'; determines whether use DICE-2016R climate dynamics or Geoffroy et al dynamics (see Dietz et al. 2020)
Params.carbonmodel = 'fair'; % 'dice','gol','joos','fair'; determines whether use DICE-2016R carbon dynamics or one of the other models (see Dietz et al. 2020)
Params.damagemodel = 'dice'; % which damage distribution: 'dice' is DICE-2016R, and 'expert' is Pindyck (2019) damages as implemented by Lemoine (2021)
%}

Params.fixsavings = -99; % <0 (default): endogenize savings rate; >=0: savings rate fixed at this value
Params.dicenegems = 1; % =1 (default): use DICE's constraint allowing negative emissions; =0: do not allow negative emissions
Params.dicecumulems = 1; % =1 (default): use DICE's constraint on cumulative emissions; =0: no constraint
Params.dicelrsavings = 1; % =1 (default): fix long-run savings rate as in DICE; =0: don't


% Computational options
Params.transitionsasconstraints = 1; % =0: solve by guessing policy and simulating trajectories; =1 (default): solve by also guessing states and treating transition equations as constraints
Params.consumptionascontrol = 0; % =1: consumption is a control and output is a constraint; =0 (default): consumption is deduced from output; only relevant if Params.transitionsasconstraints=1
Params.scaleconstraints = 1; % =1 (default): scale the constraints to be fractional deviations; =0: don't
Params.scalevars = 1; % =1 (default): scale the variables by the initial guess; =0: don't
Params.dohpc = 0; % =0: Are running on personal computer; =1: Are running remotely on high performance computing (HPC)
Params.useknitro = 0; % =1: use Knitro; =0: use fmincon
Params.knitroalgorithm = 3; % algorithm to use with Knitro, can use 1-4; default: 3
Params.fminconalgorithm = 'interior-point'; % only relevant if Params.useknitro=0 or Params.dohpc=1; default: 'interior-point'
Params.parallelize = 0; % Only relevant when Params.dohpc==1 and Params.transitionsasconstraints~=1.  =1: Parallelize gradients in fmincon.  =0: Don't
Params.screenreport = 0; % =1: report output summary to screen as optimize; =0: don't

% Obtain directory
[maindir,~,~] = fileparts(mfilename('fullpath'));
addpath(maindir);



%% Parameterization %%%%%%%%%%%%%%%

% Run parameterization script
run sub_parameters

% More computational parameters
if Params.dohpc==1 % HPC doesn't have knitro license so set up to do fmincon
    Params.useknitro = 0;
end
if Params.transitionsasconstraints==1
    if Params.useknitro
        Params.knitro_options_file = ['dice2016r_knitro_alg' num2str(Params.knitroalgorithm) '.opt'];
        Params.optimize_options_list = [];
    else
        Params.optimize_options_list = optimoptions(@fmincon,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',1e4,'MaxIterations',5e3,'Display','iter-detailed','Algorithm',Params.fminconalgorithm,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
        %Params.optimize_options_list = optimset('GradObj','on','GradConstr','on','DerivativeCheck','on');    
    end    
else
    Params.knitro_options_file = [];
    if Params.dohpc==1 && Params.parallelize==1 % Parallelize gradient calculations
        parpool(5);
        Params.optimize_options_list = optimoptions(@fmincon,'Display','off','Algorithm',Params.fminconalgorithm,'UseParallel',true);
    else
        Params.optimize_options_list = optimoptions(@fmincon,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',1e6,'MaxIterations',5e3,'Display','iter-detailed','Algorithm',Params.fminconalgorithm);
    end
end




%% Error check %%%%

if Params.fixsavings>=1
    error('Savings rate too large');
elseif Params.fixsavings>=0
    disp(['Fixing savings rate at ' num2str(Params.fixsavings)]);
end
if Params.dicelrsavings~=1 && Params.dicelrsavings~=0
    error('Params.dicelrsavings not binary.');
end

switch Params.climatemodel
    case 'dice'
        disp('DICE-2016R climate dynamics');
    case 'update'
        disp('Updated climate dynamics');
    case 'other'
        disp('Loading user-defined guesses for climate dynamics')
    otherwise
        error('Unrecognized climate model')
end

switch Params.carbonmodel
    case 'dice'
        disp('DICE-2016R carbon dynamics');
    case 'gol'
        disp('Golosov et al carbon dynamics');
    case 'joos'
        disp('Joos et al carbon dynamics');
    case 'fair'
        disp('FAIR carbon dynamics, including endogenous feedback parameter');
    case 'other'
        disp('Loading user-defined guesses for carbon dynamics')
    otherwise
        error('Unrecognized carbon model')
end

switch Params.damagemodel
    case 'dice'
        disp('Damages centered around DICE-2016R');
    case 'expert'
        disp('Damages centered around Pindyck (2019) expert survey, via Lemoine (2021)');
    case 'other'
        disp('Loading user-defined guesses for damages')
    otherwise
        error('Unrecognized damage type');
end


%% Directory Structure %%%%%%%%%%%%%%%

if Params.dohpc ~= 1 % set up directory for a run on personal computer
    
    Params.savedir = [maindir filesep 'output' filesep customdir ];    
    
    Params.savedir = [Params.savedir 'climate-' Params.climatemodel '_carbon-' Params.carbonmodel '_dam-' Params.damagemodel];
    Params.filenaming = ['cli-' Params.climatemodel '_car-' Params.carbonmodel '_dam-' Params.damagemodel];

    
    if Params.fixsavings>=0
        Params.savedir = [Params.savedir '_fixedsavings-' num2str(Params.fixsavings)];
        Params.filenaming = [Params.filenaming '_fixedsavings-' num2str(Params.fixsavings)];
    end        
    
    if Params.dicenegems~=1
        Params.savedir = [Params.savedir '_nonegems'];
        Params.filenaming = [Params.filenaming '_nonegems'];
    end
    
    if Params.dicecumulems~=1
        Params.savedir = [Params.savedir '_nocumulems'];
        Params.filenaming = [Params.filenaming '_nocumulems'];
    end
    
    if Params.dicelrsavings~=1
        Params.savedir = [Params.savedir '_freelrsavings'];
        Params.filenaming = [Params.filenaming '_freelrsavings'];
    end
    
    if Params.opt == 1
        Params.savedir = [Params.savedir '_abate-opt'];
        Params.filenaming = [Params.filenaming '_abate-opt'];
    end
    
    if Params.opt == 0
        Params.savedir = [Params.savedir '_abate-fix' num2str(fixedAbate)];
        Params.filenaming = [Params.filenaming '_abate-fix' num2str(fixedAbate)];
    end
    
    if Params.opt == 2
        Params.savedir = [Params.savedir '_abate-target' num2str(targDate)];
        Params.filenaming = [Params.filenaming '_abate-target' num2str(targDate)];
    end
    
    if Params.opt == 3
        Params.savedir = [Params.savedir '_abate-maxT' num2str(maxT)];
        Params.filenaming = [Params.filenaming '_abate-maxT' num2str(maxT)];
    end
        
    if Params.transitionsasconstraints==1
        Params.savedir = [Params.savedir '_constraint-sol'];
        Params.filenaming = [Params.filenaming '_cons-sol'];

    else
        Params.savedir = [Params.savedir '_trajectory-sol'];
        Params.filenaming = [Params.filenaming '_traj-sol'];
    end
 %{   
    if Params.useknitro==1
        Params.savedir = [Params.savedir '_knitro'];
    else
        Params.savedir = [Params.savedir '_fmincon'];
    end
 %}    
    Params.savedir = [Params.savedir filesep];
    
    if ~exist(Params.savedir, 'dir')
        mkdir(Params.savedir);
    end
    cd(Params.savedir);
%{   
    if Params.useknitro==1 && ~isempty(Params.knitro_options_file)
        copyfile([maindir filesep Params.knitro_options_file],Params.savedir);
    end
%}   
    
else % don't have permission to mkdir from Matlab on HPC, so use same directory for saving output
    
    Params.savedir = maindir;
    
end

disp(['Saving to ' Params.savedir]);

if Params.dohpc==1
    delete([Params.savedir 'iteration_report.txt']);
end



%% Optimization %%%%%%%%%%%%%%%

% load guess and variable bounds (ub and lb)
run sub_loadguesses;

% initialize iteration tracking
iteration_number = 0;

% objective and nonlinear constraints
objective = @(x) utilityobjective(x, 1, Params.horizon, Fun, Params);
nonlcon = @(x) nonlcon_utilmax(x, 1, Params.horizon, Fun, Params);
if Params.fixsavings<0
    disp('Beginning to optimize abatement and savings rates.');
else
    disp('Beginning to optimize abatement rate.');
end

% optimize
dosolve=2; % when = 2, will try at most two times
elapsedtime = 0;
tic;
while dosolve>0
    if Params.useknitro == 1
        [out_controls,fval,exitflag,outputstruct,lambdastruct,gradient,hessian] = knitromatlab(objective,guess(:),[],[],[],[],lb(:),ub(:),nonlcon,[],Params.optimize_options_list,Params.knitro_options_file);
        if exitflag<0
            disp(['Warning: Optimization may have failed, with exitflag ' num2str(exitflag)]);
            dosolve = 0;
        else
            disp('Solved successfully')
            dosolve = 0;
        end
    else
        [out_controls,fval,exitflag,outputstruct,lambdastruct,gradient,hessian] = fmincon(objective,guess(:),[],[],[],[],lb(:),ub(:),nonlcon,Params.optimize_options_list);
        if exitflag<=0
            disp(['Warning: Optimization may have failed, with exitflag ' num2str(exitflag)]);
            dosolve = 0; %dosolve - 1;
            %guess = out_controls;
            % try sqp algorithm
            %Params.optimize_options_list = optimoptions(@fmincon,'Algorithm','sqp','StepTolerance',1e-10,'ConstraintTolerance',1e-12,'MaxFunctionEvaluations',1e4,'MaxIterations',5e3,'Display','final-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
        else
            disp('Solved successfully')
            dosolve = 0;
        end
    end
    elapsedtime = elapsedtime + toc;
    disp(['Time spent optimizing: ' num2str(elapsedtime/60) ' minutes.']);
end

% store constraints
[out_cneq,out_ceq] = nonlcon(out_controls);
disp(['Max abs equality constraint is ' num2str(max(abs(out_ceq)))])

% turn back into matrix
if Params.transitionsasconstraints == 1
    out_controls = reshape(out_controls,Params.horizon,[]);
end

% make sure state variable bounds didn't bind
test_ub = ub;
test_lb = lb;
test_ub(:,[Params.col_abaterate Params.col_savingsrate]) = Inf; % don't care about controls binding, so get rid of that
test_lb(:,[Params.col_abaterate Params.col_savingsrate]) = -Inf; % don't care about controls binding, so get rid of that
ub_binding = sum( out_controls>=test_ub | abs(out_controls-test_ub)<=1e-4 , 1);
lb_binding = sum( out_controls<=test_lb | abs(out_controls-test_lb)<=1e-4 , 1);
if sum(ub_binding) > 0
    disp(['Upper bounds on states bind ' mat2str(ub_binding) ' times']);
end
if sum(lb_binding) > 0
    disp(['Lower bounds on states bind ' mat2str(lb_binding) ' times']);
end
clear test_ub test_lb;

% undo normalization
out_controls = out_controls.*Params.normalization;


%% Take Results %%%%%%%%%%%%%%%

% store the policy variables
if Params.transitionsasconstraints~=1
    out_policy = out_controls(:);
    
    out_policy2 = out_policy;
    if strcmp(Params.carbonmodel,'fair')
        alpha = out_policy2(1:Params.horizon,1);
        out_policy2(1:Params.horizon,:) = [];
    end
    abaterate = out_policy2(1:Params.horizon,1); % fraction
    if Params.fixsavings<0
        savingsrate = out_policy2(Params.horizon+1:end,1);
        if length(savingsrate)<Params.horizon
            savingsrate(end+1:Params.horizon,1) = Params.lrsavingsrate;
        end
    else
        savingsrate = Params.fixsavings*ones(Params.horizon,1);
    end
    clear out_policy2;
    
else % turn into appropriate vector    
    if Params.fixsavings>=0
        if strcmp(Params.carbonmodel,'fair')
            out_policy = [ out_controls(:,Params.col_alpha); out_controls(:,Params.col_abaterate) ];
        else
            out_policy = out_controls(:,Params.col_abaterate);
        end
    else
        if strcmp(Params.carbonmodel,'fair')
            out_policy = [ out_controls(:,Params.col_alpha); out_controls(:,Params.col_abaterate); out_controls(:,Params.col_savingsrate) ];
        else
            out_policy = [ out_controls(:,Params.col_abaterate); out_controls(:,Params.col_savingsrate) ];
        end
        if Params.dicelrsavings==1
            out_policy(end-round(5*10/Params.timestep)+1:end,:) = [];
        end                    
    end
    
    abaterate = out_controls(:,Params.col_abaterate);
    if Params.fixsavings<0
        savingsrate = out_controls(:,Params.col_savingsrate);
        if length(savingsrate)<Params.horizon
            savingsrate(end+1:Params.horizon,1) = Params.lrsavingsrate;
        end
    else
        savingsrate = Params.fixsavings*ones(Params.horizon,1);
    end
    if strcmp(Params.carbonmodel,'fair')
        alpha = out_controls(:,Params.col_alpha);
    end    
end

[ C, Ynet, K, T, M, emsind, Tocean, Ygross, abatecost] = trajectory( out_policy, 1, Params.horizon, Fun, Params );

year = 2015 + [0:Params.timestep:Params.numyears-Params.timestep]';

Welfare = sum(Params.discfactor.*Params.pop.*Fun.utility(C,Params.pop));

abatement = Fun.abatement(Params.sigma,abaterate,Ygross); % Gt C

% cumulative industrial emissions
cumulemsind = cumsum(emsind(1:Params.horizon));
cumulemsind = Params.cumulems0 + [0; cumulemsind(1:end-1,:)];

% cumulative fossil use, in Gt C
cumulfossil = cumsum(max(0,emsind));
cumulfossil = Params.cumulems0 + [0; cumulfossil(1:end-1,:)];

% implied tax on emissions, from marginal abatement cost
emtax_pertCO2 = Fun.mac(Params.sigma,abaterate,Ygross,[1:Params.horizon]'); % 2010$/tCO2
% note that have to adjust last argument if not wanting to run from first period

% atmospheric CO2, in ppm
switch Params.carbonmodel
    case 'dice' % atmospheric CO2 is the first reservoir
        Carbon_ppm = M(:,1)/Params.gtc_per_ppm;
    otherwise % atmospheric CO2 is sum of all reservoirs
        Carbon_ppm = sum(M,2)/Params.gtc_per_ppm;
end

% approximate implied alpha for fair model. Otherwise output emissions
if strcmp(Params.carbonmodel,'fair')
    otherems = [100;Fun.otherems([1:Params.horizon-1]')]; % Params.horizon x 1 vector of cumulative emissions from deforestation (Gt C); is the leading 100 Gt C meant to adjust for all emissions prior to 2015?
    alphaimplied = Fun.alpha_analytic( Fun.IRF1(T,M,cumulemsind + cumsum(otherems)) );
else
    otherems = [Params.otherems0;Fun.otherems([1:Params.horizon-1]')]; % Params.horizon x 1 vector of cumulative emissions from deforestation (Gt C); 
    
end

% format for future guess
output_for_guess(:,Params.col_abaterate) = abaterate;
output_for_guess(:,Params.col_savingsrate) = savingsrate;
if Params.consumptionascontrol
    output_for_guess(:,Params.col_C) = C;
end
output_for_guess(:,Params.col_K) = K;
output_for_guess(:,Params.col_T) = T;
output_for_guess(:,Params.col_Tocean) = Tocean;
output_for_guess(:,Params.col_M) = M;
if strcmp(Params.carbonmodel,'fair')
   output_for_guess(:,Params.col_alpha) = alpha;
end

save(['workspace.mat']);


%% Calculate Social Cost of Carbon %%%%%%%%%%%%%%%

%Pre-define SCC array
SCC = zeros(ceil(Params.horizon/2),1) ;

%Pre-allocate
W0 = zeros(ceil(Params.horizon/2),1);
W1e = zeros(ceil(Params.horizon/2),1);


%Call trajectory function to calculate new consumption vector and welfare
[ C0, ~, ~, ~, ~, ~, ~, ~, ~] = trajectory(out_policy, 1, Params.horizon, Fun, Params );

for tstep = 1:(ceil(Params.horizon/2))
   
    %FIX INDEXING - CONS SHOULD NOT BE A COMPLEX NUMBER
    %Calculate marginal welfare over half the horizon
    W0(tstep) = sum(Params.discfactor(1:(tstep +ceil(Params.horizon/2)-1)).*Params.pop(1:(tstep +ceil(Params.horizon/2)-1)).*Fun.utility(C0(1:(ceil(Params.horizon/2)+tstep-1)),Params.pop(1:(tstep +ceil(Params.horizon/2)-1)))); 

end

%Augmenting Carbon emissions and calculating SCC
for tstep = 1:(ceil(Params.horizon/2))

    %Increment emissions in period tstep
    Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + (t==tstep);
    
    %FIX INDEXING - CONS SHOULD NOT BE A COMPLEX NUMBER
    %Call trajectory function to calculate new consumption vector and welfare
    [ C_margems, ~, ~, ~, ~, ~, ~, ~, ~] = trajectory(out_policy, 1, Params.horizon, Fun, Params );
    
    %Calculate marginal welfare over half the horizon
    W1e(tstep) = sum(Params.discfactor(1:(tstep + ceil(Params.horizon/2)-1)).*Params.pop(1:(tstep + ceil(Params.horizon/2)-1)).*Fun.utility(C_margems(1:(tstep + ceil(Params.horizon/2)-1)),Params.pop(1:(tstep + ceil(Params.horizon/2)-1)))); 
    
    %Calculate partial derivative of welfare with respect to consumption
    dWc = Params.discfactor(tstep)*Params.pop(tstep)*Fun.dutility_dC(C(tstep),Params.pop(tstep));
    
    %Calculate SCC
    SCC(tstep,1) = -((W1e(tstep) - W0(tstep))/dWc)*1000/Params.co2_per_c ;  

end

%Restore otherems function
Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5);


%% Calculate Social Cost of Fire %%%%%%%%%%%%%%%

%Pre-define arrays
% W0 = zeros(ceil(Params.horizon/2),1);
% SCfire1 = zeros(ceil(Params.horizon/2),1) ;
% W1fire1 = zeros(ceil(Params.horizon/2),1);
C_diff_ram = zeros(ceil(Params.horizon/2),1);
C_diff_25 = zeros(ceil(Params.horizon/2),1);
C_diff_3 = zeros(ceil(Params.horizon/2),1);
C_diff_5 = zeros(ceil(Params.horizon/2),1);
em_diff = zeros(ceil(Params.horizon/2),1);
emSCC_ram = zeros(ceil(Params.horizon/2),1);
emSCC_25 = zeros(ceil(Params.horizon/2),1);
emSCC_3 = zeros(ceil(Params.horizon/2),1);
emSCC_5 = zeros(ceil(Params.horizon/2),1);


%Call trajectory function to calculate new consumption vector and welfare
[ C0, ~, ~, ~, ~, ~, ~, ~, ~] = trajectory(out_policy, 1, Params.horizon, Fun, Params );

run OutputResults
iR1 = iR + 1 ;

%Calculate marginal welfare over half the horizon
W0 = sum(Params.discfactor(1:ceil(Params.horizon/2)).*Params.pop(1:ceil(Params.horizon/2)).*Fun.utility(C0(1:ceil(Params.horizon/2)),Params.pop(1:ceil(Params.horizon/2)))); 

% %mush the regrowth around for robustness
% %start by putting all regrowth in the first period...
% SCfire_movedregrowth = zeros(59, 1)
% for i=8:59
% Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + ((t==1)*(20141.17*10/1000000000))...
%     + (t==1)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + ((t>1)&(t<=7))*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + (t==i)*(-2.4206e-04);
% %Call trajectory function to calculate new consumption vector and welfare
% [ C_fire1, Ynet_fire1, K_fire1, T_fire1, M_fire1, emsind_fire1, Tocean_fire1, Ygross_fire1, abatecost_fire1] = trajectory(out_policy, 1, Params.horizon, Fun, Params );
% %ReCalculate marginal welfare over half the horizon
% W1fire1 = sum(Params.discfactor(1:ceil(Params.horizon/2)).*Params.pop(1:ceil(Params.horizon/2)).*Fun.utility(C_fire1(1:ceil(Params.horizon/2)),Params.pop(1:ceil(Params.horizon/2)))); 
% %Calculate partial derivative of welfare with respect to consumption
% dWc = Params.discfactor(1)*Params.pop(1)*Fun.dutility_dC(C0(1),Params.pop(1));
% %Calculate the discounted welfare difference
% SCfire1 = -((W1fire1 - W0)/dWc)*1000000000000*(258.8/218.1) ;  %convert from trillions of dollars to dollars
% %store result to vector
% SCfire_movedregrowth(i) = SCfire1
% end

%primary model -- law et al. high-severity
%change emissions to incorporate fire
Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + ((t==1)*(20141.17*10/1000000000))...
    + (t==1)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
    + ((t>1)&(t<=58))*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
    + (t==59)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2)))); 
    %7 5-year periods for postivie emissions only
    %58 periods for full vector, plus most of 59


% %primary model -- law et al. high-severity
% %change emissions to incorporate fire
% %250 years of original fire, followed by another ignition
% Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + ((t==1)*(20141.17*10/1000000000))...
%     + (t==1)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + ((t>1)&(t<=50))*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + (t==51)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)-1))/95.63)/1.098).^2))))...
%     + ((t==50)*(20141.17*10/1000000000))...
%     + (t==50)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + ((t>50)&(t<=58))*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + (t==59)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2))))...; 
%     %7 5-year periods for postivie emissions only
%     %58 periods for full vector, plus most of 59

% %medium severity -- factor of 0.6
% Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + ((t==1)*(0.6)*(20141.17*10/1000000000))...
%     + (t==1)*(0.6)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + ((t>1)&(t<=58))*(0.6)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + (t==59)*(0.6)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2)))); 

% %low severity -- factor of 0.35
% Fun.otherems = @(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + ((t==1)*(0.35)*(20141.17*10/1000000000))...
%     + (t==1)*(0.35)*(-10/1000000000)*((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + ((t>1)&(t<=58))*(0.35)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+2))/95.63)/1.098).^2))))...
%     + (t==59)*(0.35)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2)))); 

%Boreal fire emissions vector
% Fun.otherems =@(t) Params.otherems0*(1-Params.gotherems).^((t-1)*Params.timestep/5) + (t==1)*(3325*10/1000000000)...
%        + (t==1)*(-10/1000000000)*sum(-181.99+0.0016*((5*(t-1)+2)^3)-0.304*((5*(t-1)+2)^2)+16.548*((5*(t-1)+2)))... 
%        + ((t>1)&(t<=3))*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)+1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)+1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)+1))))...
%        + (t==4)*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1))).^3)-0.304*(((5*(t-1)-3):(5*(t-1))).^2)+16.548*(((5*(t-1)-3):(5*(t-1)))))...
%        + (t==4)*(-10/1000000000)*sum(-181.99+0.0016*((5*(t-1)+1).^3)-0.304*((5*(t-1)+1).^2)+16.548*(5*(t-1)+1))...
%        + ((t>4)&(t<=17))*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)+1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)+1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)+1))))...
%        + (t==18)*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)-1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)-1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)-1))))...
%        + (t==18)*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)+1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)+1))))...
%        + ((t>18)&(t<=30))*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)+1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)+1))))...
%        + (t==31)*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)-1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)-1))));
%     %fire vector from Goulden et al. 2011, via Phillips et al. 2022
%     %switches from positive to negative emissions in year 15
    %units = gigatons C 


%Call trajectory function to calculate new consumption vector and welfare
[ C_fire1, Ynet_fire1, K_fire1, T_fire1, M_fire1, emsind_fire1, Tocean_fire1, Ygross_fire1, abatecost_fire1] = trajectory(out_policy, 1, Params.horizon, Fun, Params );

%ReCalculate marginal welfare over half the horizon
W1fire1 = sum(Params.discfactor(1:ceil(Params.horizon/2)).*Params.pop(1:ceil(Params.horizon/2)).*Fun.utility(C_fire1(1:ceil(Params.horizon/2)),Params.pop(1:ceil(Params.horizon/2)))); 
%Calculate partial derivative of welfare with respect to consumption
dWc = Params.discfactor(1)*Params.pop(1)*Fun.dutility_dC(C0(1),Params.pop(1));

%Scenario 1 - Calculate the discounted welfare difference (look at period 1)
SCfire1 = -((W1fire1 - W0)/dWc)*1000000000000*(258.8/218.1) ;  %convert from trillions of 2010 dollars to single 2020 dollars
    


%CPI conversion 2010 --> 2020: (258.8/218.1) = 1.1866


%Augmenting Carbon emissions and calculating SCF for IWG and Naive
%Approaches 
for tstep = 1:(ceil(Params.horizon/2))

    em_diff(tstep) =(tstep==1)*(20141.17*10/1000000000)...
        + (tstep==1)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)+2))/95.63)/1.098).^2))))...
        + ((tstep>1)&(tstep<=58))*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)-3):(5*(tstep-1)+1))/95.63)/1.098).^2))))...
        + (t==59)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2))));
    %fire vector from Law et al. 2003
    %units = gigatons C 

%     %moderate-severity -- factor of 0.6
%     em_diff(tstep) =(tstep==1)*(0.6)*(20141.17*10/1000000000)...
%         + (tstep==1)*(0.6)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)+2))/95.63)/1.098).^2))))...
%         + ((tstep>1)&(tstep<=58))*(0.6)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)-3):(5*(tstep-1)+1))/95.63)/1.098).^2))))...
%         + (t==59)*(0.6)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2))));
    
    %low-severity -- factor of 0.35
%     em_diff(tstep) =(tstep==1)*(0.35)*(20141.17*10/1000000000)...
%         + (tstep==1)*(0.35)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)+2))/95.63)/1.098).^2))))...
%         + ((tstep>1)&(tstep<=58))*(0.35)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((5*(tstep-1)-3):(5*(tstep-1)+1))/95.63)/1.098).^2))))...
%         + (t==59)*(0.35)*(-10/1000000000)*sum((-239.4+399.9*exp(-0.5*((log(((Params.timestep*(t-1)-2):(Params.timestep*(t-1)+1))/95.63)/1.098).^2))));

%Phillips curve (boreal)
%     em_diff(tstep) =(tstep==1)*(3325*10/1000000000); %...
%         + (t==1)*(-10/1000000000)*sum(-181.99+0.0016*((5*(t-1)+2)^3)-0.304*((5*(t-1)+2)^2)+16.548*((5*(t-1)+2)))... 
%        + ((t>1)&(t<=3))*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)+1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)+1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)+1))))...
%        + (t==4)*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1))).^3)-0.304*(((5*(t-1)-3):(5*(t-1))).^2)+16.548*(((5*(t-1)-3):(5*(t-1)))))...
%        + (t==4)*(-10/1000000000)*sum(-181.99+0.0016*((5*(t-1)+1).^3)-0.304*((5*(t-1)+1).^2)+16.548*(5*(t-1)+1))...
%        + ((t>4)&(t<=17))*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)+1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)+1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)+1))))...
%        + (t==18)*(-10/1000000000)*sum(-181.99+0.0016*(((5*(t-1)-3):(5*(t-1)-1)).^3)-0.304*(((5*(t-1)-3):(5*(t-1)-1)).^2)+16.548*(((5*(t-1)-3):(5*(t-1)-1))))...
%        + (t==18)*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)+1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)+1))))...
%        + ((t>18)&(t<=30))*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)+1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)+1))))...
%        + (t==31)*(-10/1000000000)*sum(29.439+0.0005*(((5*(t-1)):(5*(t-1)-1)).^2)-0.2628*(((5*(t-1)):(5*(t-1)-1))));
    %fire vector from Goulden et al. 2011, via Phillips et al. 2022
    %units = gigatons C 

    %Scenario 2 - Calculate the discounted consumption difference (look at period 1)
        %discount rate - Ramsey
    C_diff_ram(tstep) = (C_fire1(tstep) - C0(tstep))/prod(iR1(1:(tstep-1)).^Params.timestep) ;
        %discount rate - 2.5%
    C_diff_25(tstep) = (C_fire1(tstep) - C0(tstep))/(1.025^((tstep-1)*Params.timestep)) ;
            %discount rate - 3%
    C_diff_3(tstep) = (C_fire1(tstep) - C0(tstep))/(1.03^((tstep-1)*Params.timestep)) ;
                %discount rate - 5%
    C_diff_5(tstep) = (C_fire1(tstep) - C0(tstep))/(1.05^((tstep-1)*Params.timestep)) ;

    %Scenario 3 - Multiply SCC by emissions difference  by carbon -> CO2 conversion by discount factor
        %discount rate - Ramsey
    emSCC_ram(tstep) = (em_diff(tstep)*SCC(tstep)*1000000000*Params.co2_per_c)/prod(iR1(1:(tstep-1)).^Params.timestep) ;
        %discount rate - 2.5%
    emSCC_25(tstep) = (em_diff(tstep)*SCC(tstep)*1000000000*Params.co2_per_c)/(1.025^((tstep-1)*Params.timestep)) ;
        %discount rate - 3%
    emSCC_3(tstep) = (em_diff(tstep)*SCC(tstep)*1000000000*Params.co2_per_c)/(1.03^((tstep-1)*Params.timestep)) ;
        %discount rate - 5%
    emSCC_5(tstep) = (em_diff(tstep)*SCC(tstep)*1000000000*Params.co2_per_c)/(1.05^((tstep-1)*Params.timestep)) ;
end

%final Scenario 2 calculation
c_diff_ram_tot = -sum(C_diff_ram)*1000000000000*(258.8/218.1) ;  %convert from trillions of dollars to dollars
c_diff_25_tot = -sum(C_diff_25)*1000000000000*(258.8/218.1) ;  %convert from trillions of dollars to dollars
c_diff_3_tot = -sum(C_diff_3)*1000000000000*(258.8/218.1) ;  %convert from trillions of dollars to dollars
c_diff_5_tot = -sum(C_diff_5)*1000000000000*(258.8/218.1) ;  %convert from trillions of dollars to dollars
%Scenario 3 - emissions times SCC
emSCC_ram_tot = sum(emSCC_ram)*(258.8/218.1);
emSCC_25_tot = sum(emSCC_25)*(258.8/218.1);
emSCC_3_tot = sum(emSCC_3)*(258.8/218.1);
emSCC_5_tot = sum(emSCC_5)*(258.8/218.1);

%% Output Results to Excel File %%%%%%%%%%%%%%%

%Change Directory Back to Main Folder
cd(maindir);

%Run Output Script 
run OutputResults

save workspace_main.mat










