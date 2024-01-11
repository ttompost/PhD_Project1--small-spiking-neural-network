function net = N4_arch(eNum_perc, iNum_perc, vpmNum, timeParams, synParams, vpmParams, varySim)
% Function should be called with following inputs:
% (1) eNum_perc:    percentage of E neurons per layer
% (2) iNum_perc:    percentage of E neurons per layer
% (3) vpmNum:       number of VPM neurons
% (4) timeParams:   time in ms, [start stop step]
% (5) synParams:    synaptic statistics, synParams.synstr = [E_syn_str, I_syn_str]; synParams.pconn = [E_pconn, I_pconn];
% (6) vpmParams:    thalamic population parameters, [correlated_pop_rate, correlated_pop_input_prob, uncorrelated_pop_rate, uncorrelated_pop_input_prob, input_conductance]
% optional (7):     parameters to vary

% The function creates a thalamocortical model and runs one simulation with
% assigned temporal parameters.

if nargin < 7
    varySim = {};
end

%%%%%%%%%%%%%%%%%%
% I use this paper to set up number of neurons per layer:
% "The resulting neuron numbers per layer in the average [rat barrel] column
% were (...) 2039 ± 524 (L2), 3735 ± 905 (L3) and 4447 6 439 (L4)." 
% Meyer, Hanno S., Verena C. Wimmer, M. Oberlaender, Christiaan P.J. de Kock, Bert Sakmann, and Moritz Helmstaedter. 
% “Number and Laminar Distribution of Neurons in a Thalamocortical Projection Column of Rat Vibrissal Cortex.” 
% Cerebral Cortex (New York, NY) 20, no. 10 (October 2010): 2277–86. https://doi.org/10.1093/cercor/bhq067.

ctx_downscale = 20; % downscale total number of neurons per layer
l2e_num = floor((2039*eNum_perc)/ctx_downscale);
l2i_num = floor((2039*iNum_perc)/ctx_downscale);
l3e_num = floor((3735*eNum_perc)/ctx_downscale);
l3i_num = floor((3735*iNum_perc)/ctx_downscale);
l4e_num = floor((4447*eNum_perc)/ctx_downscale);
l4i_num = floor((4447*iNum_perc)/ctx_downscale);

%%%%%%%%%%%%%%%%%%
AMPA_downscale = 6; % scale i(AMPA) so that uEPSP is 1-2 mV
GABAa_downscale = 2; % scale i(GABAa) so that uIPSP is 4-5 mV

% uPSP reference values:
% Sun, Qian-Quan, John R. Huguenard, and David A. Prince. “Barrel Cortex Microcircuits: Thalamocortical 
% Feedforward Inhibition in Spiny Stellate Cells Is Mediated by a Small Number of Fast-Spiking Interneurons.” 
% The Journal of Neuroscience 26, no. 4 (January 25, 2006): 1219–30. https://doi.org/10.1523/JNEUROSCI.4727-04.2006.

%%%%%%%%%%%%%%%%%%
% Hodgking-Huxley neuron models were taken from here:
% Pospischil, Martin, Maria Toledo-Rodriguez, Cyril Monier, Zuzanna Piwkowska, Thierry Bal, Yves Frégnac, Henry Markram, 
% and Alain Destexhe. “Minimal Hodgkin–Huxley Type Models for Different Classes of Cortical and Thalamic Neurons.” 
% Biological Cybernetics 99, no. 4 (November 1, 2008): 427–41. https://doi.org/10.1007/s00422-008-0263-8.

somaDiameter_RS = 61.4; %um
somaDiameter_FS = 56.9; %um

% I suppose that neurons are spherical, so I take their diameter (from the
% paper above) and get their total surface (surface of the sphere formula:
% A=pi*D^2). using their surface will influence every parameter that
% depends on the cell size (e.g., capacitance, channel conductances etc.)

surf_RS = pi*somaDiameter_RS^2; % um2
surf_FS = pi*somaDiameter_FS^2; % um2

%%%%%%%%%%%%%%%%%%
% extract connectivity statistics per E and I neuron type from input
% structures

% synaptic strength
ss=getfield(synParams,'synstr'); % synaptic strength
synstr_e = ss(1);
synstr_i = ss(2);

% probability of connectivity
pc=getfield(synParams,'pconn'); % probability of connectivity
pconn_e = pc(1);
pconn_i = pc(2);

% extract simulation time from input vector
tspan = [timeParams(1) timeParams(2)];     % time span [begin end step], ms
dt = timeParams(3);                        % fixed time step, ms

%%%% define parameters for E and I neurons reported in Pospischil2008 (reference above)

% RS Pospischil2008 params // supposed E cells
rsparam = {...
    'Cm', 1*1e-4 * surf_RS,...
    'gleak',2.05*1e-5*1e3 * 1e-4 * surf_RS,... % mS/cm2 converted into mS/um2 * area in um2
    'Eleak',-70.3,...
    'gNa',56 * 1e-4 * surf_RS,...
    'Na_V1',-65,... % -56.2
    'gK',6 * 1e-4 * surf_RS,...
    'K_V1',-56.2,...
    'gM',7.5*1e-5*1e3 * 1e-4 * surf_RS,... % mS/cm2
    'mM_tau_max',608,... % ms
    };

% FS Pospischil2008 params // supposed I cells
fsparam = {...
    'Cm', 1 * 1e-4 * surf_FS,...
    'gleak',3.8*1e-5*1e3 * 1e-4 * surf_FS,... % mS/cm2
    'Eleak',-70.4,...
    'gNa',58 * 1e-4 * surf_FS,...
    'Na_V1',-57.9,...
    'gK',3.9 * 1e-4 * surf_FS,...
    'K_V1',-57.9,...
    'gM',7.87*1e-5*1e3 * 1e-4 * surf_FS,... % mS/cm2
    'mM_tau_max',502,... % ms
    };


%%%% make the network skeleton (define layers, populate with neurons, add intrinsic noise, connect populations)
% L2
L2.populations(1).name = 'L2E';
L2.populations(end).size = l2e_num;
L2.populations(end).equations = 'RS_Pospischil2008';
L2.populations(end).mechanism_list={'noise'};
L2.populations(end).parameters=[rsparam, 'noise_amp',100];

L2.populations(end+1).name = 'L2I';
L2.populations(end).size = l2i_num;
L2.populations(end).equations = 'FS_Pospischil2008';
L2.populations(end).mechanism_list={'noise'};
L2.populations(end).parameters=[fsparam, 'noise_amp',100];


% connect populations in L2 layer:
% calculate connectivity matrices: 
    % conmat=ConnMx(activity, layer, vary, Npre, Npost, Pconn, Nb, sf)

l2EEnetcon=ConnMx(0,'intrahom','none',l2e_num,l2e_num,ConnMxParam('l2ee','pconn', pconn_e),ConnMxParam('l2ee','nb', synstr_e),1); % l2E to l2E % removed autoptic
l2EInetcon=ConnMx(0,'intrahet','none',l2e_num,l2i_num,ConnMxParam('l2ei','pconn', pconn_e),ConnMxParam('l2ei','nb', synstr_e),1); % l2E to l2I
l2IEnetcon=ConnMx(0,'intrahet','none',l2i_num,l2e_num,ConnMxParam('l2ie','pconn', pconn_i),ConnMxParam('l2ie','nb', synstr_i),1); % l2I to l2E 
l2IInetcon=ConnMx(0,'intrahom','none',l2i_num,l2i_num,ConnMxParam('l2ii','pconn', pconn_i),ConnMxParam('l2ii','nb', synstr_i),1); % l2I to l2I % removed autoptic

% connect populations in L2 layer 
L2.connections(1).direction='L2E->L2E'; % E to E
L2.connections(end).mechanism_list={'iAMPA'};
L2.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l2EEnetcon};

L2.connections(end+1).direction='L2E->L2I'; % E to I
L2.connections(end).mechanism_list={'iAMPA'};
L2.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l2EInetcon}; 

L2.connections(end+1).direction='L2I->L2E'; % I to E
L2.connections(end).mechanism_list={'iGABAa'};
L2.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l2IEnetcon}; 

L2.connections(end+1).direction='L2I->L2I'; % I to I
L2.connections(end).mechanism_list={'iGABAa'};
L2.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l2IInetcon}; 
L2=dsCheckSpecification(L2); 

% L3
L3.populations(1).name = 'L3E';
L3.populations(end).size = l3e_num;
L3.populations(end).equations = 'RS_Pospischil2008';
L3.populations(end).mechanism_list={'noise'};
L3.populations(end).parameters=[rsparam, 'noise_amp',100];

L3.populations(end+1).name = 'L3I';
L3.populations(end).size = l3i_num;
L3.populations(end).equations = 'FS_Pospischil2008';
L3.populations(end).mechanism_list={'noise'};
L3.populations(end).parameters=[fsparam, 'noise_amp',100];

% modify layers with params from Huang,C.(2016)
    % l3 modifications: connection params
l3EEnetcon=ConnMx(0,'intrahom','none',l3e_num,l3e_num,ConnMxParam('l3ee','pconn', pconn_e),ConnMxParam('l3ee','nb', synstr_e),1); % l2E to l2E % removed autoptic
l3EInetcon=ConnMx(0,'intrahet','none',l3e_num,l3i_num,ConnMxParam('l3ei','pconn', pconn_e),ConnMxParam('l3ei','nb', synstr_e),1); % l2E to l2I
l3IEnetcon=ConnMx(0,'intrahet','none',l3i_num,l3e_num,ConnMxParam('l3ie','pconn', pconn_i),ConnMxParam('l3ie','nb', synstr_i),1); % l2I to l2E 
l3IInetcon=ConnMx(0,'intrahom','none',l3i_num,l3i_num,ConnMxParam('l3ii','pconn', pconn_i),ConnMxParam('l3ii','nb', synstr_i),1); % l2I to l2I % removed autoptic

% connect populations in L3 layer 
L3.connections(1).direction='L3E->L3E'; % E to E
L3.connections(end).mechanism_list={'iAMPA'};
L3.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l3EEnetcon};

L3.connections(end+1).direction='L3E->L3I'; % E to I
L3.connections(end).mechanism_list={'iAMPA'};
L3.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l3EInetcon}; 

L3.connections(end+1).direction='L3I->L3E'; % I to E
L3.connections(end).mechanism_list={'iGABAa'};
L3.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l3IEnetcon}; 

L3.connections(end+1).direction='L3I->L3I'; % I to I
L3.connections(end).mechanism_list={'iGABAa'};
L3.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l3IInetcon}; 
L3=dsCheckSpecification(L3); 


% L4
L4.populations(1).name = 'L4E';
L4.populations(end).size = l4e_num;
L4.populations(end).equations = 'RS_Pospischil2008';
L4.populations(end).mechanism_list={'noise'};
L4.populations(end).parameters=[rsparam, 'noise_amp',100];

L4.populations(end+1).name = 'L4I';
L4.populations(end).size = l4i_num;
L4.populations(end).equations = 'FS_Pospischil2008';
L4.populations(end).mechanism_list={'noise'};
L4.populations(end).parameters=[fsparam, 'noise_amp',100];

%     % l4 modifications: connection params
l4EEnetcon=ConnMx(0,'intrahom','none',l4e_num,l4e_num,ConnMxParam('l4ee','pconn', pconn_e),ConnMxParam('l4ee','nb', synstr_e)*0.5,1); 
l4EInetcon=ConnMx(0,'intrahet','none',l4e_num,l4i_num,ConnMxParam('l4ei','pconn', pconn_e),ConnMxParam('l4ei','nb', synstr_e),1); 
l4IEnetcon=ConnMx(0,'intrahet','none',l4i_num,l4e_num,ConnMxParam('l4ie','pconn', pconn_i),ConnMxParam('l4ie','nb', synstr_i),1); 
l4IInetcon=ConnMx(0,'intrahom','none',l4i_num,l4i_num,ConnMxParam('l4ii','pconn', pconn_i),ConnMxParam('l4ii','nb', synstr_i),1); 


% connect populations in L4 layer 
L4.connections(1).direction='L4E->L4E'; % E to E
L4.connections(end).mechanism_list={'iAMPA'};
L4.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l4EEnetcon};

L4.connections(end+1).direction='L4E->L4I'; % E to I
L4.connections(end).mechanism_list={'iAMPA'};
L4.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l4EInetcon}; 

L4.connections(end+1).direction='L4I->L4E'; % I to E
L4.connections(end).mechanism_list={'iGABAa'};
L4.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l4IEnetcon}; 

L4.connections(end+1).direction='L4I->L4I'; % I to I
L4.connections(end).mechanism_list={'iGABAa'};
L4.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l4IInetcon}; 
L4=dsCheckSpecification(L4); 

% MAKE A THALAMIC LAYER:
    % there are no Inhibitory population in Th 
    % this population presents neurons of VPM thalamic nucleus
    
% arbitrary proportions
corr_contribution = 0.5714;
uncorr_contribution = 0.4286;

vpmNum1=floor(vpmNum*corr_contribution);
vpmNum2=floor(vpmNum*uncorr_contribution);

VPMuncorr.populations(1).name='VPMuncorr';
VPMuncorr.populations(end).size=vpmNum2;
VPMuncorr.populations(end).equations={'dV/dt=@current'};
VPMuncorr.populations(end).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                   'CaBufferChing2010TC','iTChing2010TC',...
                                   'iHChing2010TC','iLeakChing2010TC',...
                                   'iKLeakChing2010TC'};
VPMuncorr.populations(end).parameters={'gH',0.005};
VPMuncorr.connections(1).direction='VPMuncorr->VPMuncorr';
VPMuncorr.connections(end).mechanism_list={'iPSTUncorr'};
VPMuncorr.connections(end).parameters={'g_esyn_uncorr',vpmParams(5)*1.5,...
                                'uncorr_rate',vpmParams(3),...
                                'uncorr_prob_cxn',vpmParams(4),...
                                'T',tspan(end),'jitter_stddev_corr',500,'tau_i_uncorr',25};

VPMcorr.populations(1).name='VPMcorr'; 
VPMcorr.populations(end).size=vpmNum1;
VPMcorr.populations(end).equations={'dV/dt=@current'};
VPMcorr.populations(end).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                   'CaBufferChing2010TC','iTChing2010TC',...
                                   'iHChing2010TC','iLeakChing2010TC',...
                                   'iKLeakChing2010TC'};
VPMcorr.populations(end).parameters={'gH',0.01};
VPMcorr.connections(1).direction='VPMcorr->VPMcorr';
VPMcorr.connections(end).mechanism_list={'iPSTCorr'};
VPMcorr.connections(end).parameters={'g_esyn_corr',vpmParams(5),...
                                'corr_rate',vpmParams(1),...
                                'corr_prob_cxn',vpmParams(2),...
                                'T',tspan(end),'jitter_stddev_corr',500,'tau_i_corr',14};

VPMcorr=dsCheckSpecification(VPMcorr); 
VPMuncorr=dsCheckSpecification(VPMuncorr); 

% make 4-layered network
Network=dsCombineSpecifications(L2,L3,L4,VPMcorr,VPMuncorr);

% make connectivity matrices for inter-layer connections: netcon=ConnMx(activity,layer,vary,Npre,Npost,va,Nb)
    % Th -> L4

thl4EEnetcon2=ConnMx(0,'trans','none',vpmNum2,l4e_num,ConnMxParam('thl4e','pconn', pconn_e)*uncorr_contribution,ConnMxParam('thl4e','nb', pconn_e),1); % uncorr
thl4EInetcon2=ConnMx(0,'trans','none',vpmNum2,l4i_num,ConnMxParam('thl4i','pconn', pconn_e)*uncorr_contribution,ConnMxParam('thl4i','nb', pconn_e),1); % uncorr
thl4EEnetcon1=ConnMx(0,'trans','none',vpmNum1,l4e_num,ConnMxParam('thl4e','pconn', pconn_e)*corr_contribution,ConnMxParam('thl4e','nb', pconn_e),1); % corr
thl4EInetcon1=ConnMx(0,'trans','none',vpmNum1,l4i_num,ConnMxParam('thl4i','pconn', pconn_e)*corr_contribution,ConnMxParam('thl4i','nb', pconn_e),1); % corr
%     % L4 -> L3
l4l3EEnetcon=ConnMx(0,'trans','none',l4e_num,l3e_num,ConnMxParam('l4el3e','pconn', pconn_e),ConnMxParam('l4el3e','nb', synstr_e),1);
l4l3EInetcon=ConnMx(0,'trans','none',l4e_num,l3i_num,ConnMxParam('l4el3i','pconn', pconn_e),ConnMxParam('l4el3i','nb', synstr_e),1);
l4l3IEnetcon=ConnMx(0,'trans','none',l4i_num,l3e_num,ConnMxParam('l4il3e','pconn', pconn_i),ConnMxParam('l4il3e','nb', synstr_i),1);
    % L3 -> L2
l3l2EEnetcon=ConnMx(0,'trans','none',l3e_num,l2e_num,ConnMxParam('l3el2e','pconn', pconn_e),ConnMxParam('l3el2e','nb', synstr_e),1);
l3l2EInetcon=ConnMx(0,'trans','none',l3e_num,l2i_num,ConnMxParam('l3el2i','pconn', pconn_e),ConnMxParam('l3el2i','nb', synstr_e),1);
l3l2IEnetcon=ConnMx(0,'trans','none',l3i_num,l2e_num,ConnMxParam('l3il2e','pconn', pconn_i),ConnMxParam('l3il2e','nb', synstr_i),1);
l3l2IInetcon=ConnMx(0,'trans','none',l3i_num,l2i_num,ConnMxParam('l3il2i','pconn', pconn_i),ConnMxParam('l3il2i','nb', synstr_i),1);

%-----------------------------------------------------------------------------------------------------------
% add feedforward connections (th->l4; l4->l3; l3->l2)
    % thalamic to l4
    
th_ampa_downscale = 2;
% AMPA_downscale is doubled so that uEPSP is ~1 mV, 
% because this paper (Bruno and Sakmann (2006) Science) shows
% that thalamocortical synapses have very low efficacy, but manage to
% trigger cortical response because they are synchronous.. they show
% majority of uEPSPs are ~ 0.5 mV (I tested this in
% VPM_activity_fitting.mlx)

Network.connections(end+1).direction='VPMcorr->L4E'; % ThE to l4E
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale*th_ampa_downscale,'netcon',thl4EEnetcon1};
 
Network.connections(end+1).direction='VPMcorr->L4I'; % ThE to l4I
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale*th_ampa_downscale,'netcon',thl4EInetcon1}; 

Network.connections(end+1).direction='VPMuncorr->L4E'; % ThE to l4E
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale*th_ampa_downscale,'netcon',thl4EEnetcon2};
 
Network.connections(end+1).direction='VPMuncorr->L4I'; % ThE to l4I
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale*th_ampa_downscale,'netcon',thl4EInetcon2};

    % l4 to l23 (l4 to l3)
Network.connections(end+1).direction='L4E->L3E'; % l4E to l3E
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l4l3EEnetcon};

Network.connections(end+1).direction='L4E->L3I'; % l4E to l3I
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l4l3EInetcon};

Network.connections(end+1).direction='L4I->L3E'; % l4I to l3E
Network.connections(end).mechanism_list={'iGABAa'};
Network.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l4l3IEnetcon}; 

    % l23 to l23 feedforward (l3 to l2)
Network.connections(end+1).direction='L3E->L2E'; % l3E to l2E
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l3l2EEnetcon}; 

Network.connections(end+1).direction='L3E->L2I'; % l3E to l2I
Network.connections(end).mechanism_list={'iAMPA'};
Network.connections(end).parameters={'downscaleFactor',AMPA_downscale,'netcon',l3l2EInetcon};

Network.connections(end+1).direction='L3I->L2E'; % l3I to l2E
Network.connections(end).mechanism_list={'iGABAa'};
Network.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l3l2IEnetcon}; 

Network.connections(end+1).direction='L3I->L2I'; % l3I to l2I
Network.connections(end).mechanism_list={'iGABAa'};
Network.connections(end).parameters={'downscaleFactor',GABAa_downscale,'netcon',l3l2IInetcon}; 

Network=dsCheckSpecification(Network); 

%SIMULATION:
net = dsSimulate(Network,'tspan',tspan,'dt',dt,'vary',varySim);
end