function net = SimN4(eNum_perc, iNum_perc, VPM_num, synParams, ang, curv, input_mode, seed_num, endSim, sparse_info)

% FINAL VERSION (20th March 2024) // Tompos, T. (t.tompos@neurophysiology.nl)

% Function should be called with following inputs:
% (1) eNum_perc:    percentage of E neurons per layer
% (2) iNum_perc:    percentage of E neurons per layer
% (3) VPM_num:      number of VPM neurons
% (4) synParams:    synaptic statistics, synParams.synstr = [E_syn_str, I_syn_str]; synParams.pconn = [E_pconn, I_pconn];
% (5) ang:          whisker angle data from Peron et al (2015) // it has to be structured so that the number of traces matches the number of VPM neurons
% (6) curv:         whisker curvature data from Peron et al (2015) // it has to be structured so that the number of traces matches the number of VPM neurons
% (7) input_mode:   pass 'exp' to stimulate the VPM with experimental data (this is where ang and curv have to be properly defined); if something else is passed, the script runs with synthetic thalamus (predefined parameters)
% (8) seed_num:     which rng() to use while generating connectivity matrices and/or thalamic activity
% (9) endSim:       length of the simulation, in ms

    if nargin < 10 || isempty(sparse_info)
        sparsify = 0;
        vary_conn = 'none';
        sparse_fact = 1;
    else
        sparsify = sparse_info{1};
        vary_conn = sparse_info{2};
        sparse_fact = sparse_info{3};

        fprintf('Sparsifying %i...', sparsify)
    end
    
    if ~exist('seed_num','var')
        ctx_seed = [];
        vpm_seed = [];
        rng('shuffle')
    else
        ctx_seed = seed_num.ctx_seed;
        vpm_seed = seed_num.vpm_seed;
    end

    simulation_dt = 0.01; % ms

    %%%%%%%%%%%%%%%%%%
    % "The resulting neuron numbers per layer in the average [rat barrel] column
    % were (...) 2039 ± 524 (L2), 3735 ± 905 (L3) and 4447 6 439 (L4)." 
    % Meyer, Hanno S., Verena C. Wimmer, M. Oberlaender, Christiaan P.J. de Kock, Bert Sakmann, and Moritz Helmstaedter. 
    % “Number and Laminar Distribution of Neurons in a Thalamocortical Projection Column of Rat Vibrissal Cortex.” 
    % Cerebral Cortex (New York, NY) 20, no. 10 (October 2010): 2277–86. https://doi.org/10.1093/cercor/bhq067.

    ctx_downscale = 20;
    l2e_num = floor((2039*eNum_perc)/ctx_downscale);
    l2i_num = floor((2039*iNum_perc)/ctx_downscale);
    l3e_num = floor((3735*eNum_perc)/ctx_downscale);
    l3i_num = floor((3735*iNum_perc)/ctx_downscale);
    l4e_num = floor((4447*eNum_perc)/ctx_downscale);
    l4i_num = floor((4447*iNum_perc)/ctx_downscale);

    %%%%%%%%%%%%%%%%%%
    % corticocortical synaptic parameters
    AMPA_scale = 0.17;
    GABAa_scale = 0.5;

    gAMPA = 0.1 * AMPA_scale;
    gGABAa = 0.25 * GABAa_scale;
    tauD_GABAa = 10;
    EGABAa = -80;

    IC = 1e-5;

    % uPSP reference values:
    % Sun, Qian-Quan, John R. Huguenard, and David A. Prince. “Barrel Cortex Microcircuits: Thalamocortical 
    % Feedforward Inhibition in Spiny Stellate Cells Is Mediated by a Small Number of Fast-Spiking Interneurons.” 
    % The Journal of Neuroscience 26, no. 4 (January 25, 2006): 1219–30. https://doi.org/10.1523/JNEUROSCI.4727-04.2006.

    %%% scaling of the thalamocortical EPSPs
    gAMPA_th = gAMPA*2;

    % thalamocortical uEPSP should be decreased to half of what the corticocortical uEPSP amplitude is, 
    % because this paper (Bruno and Sakmann (2006) Science) shows
    % that thalamocortical synapses have very low efficacy, but manage to
    % trigger cortical response because they are synchronous.. they show
    % majority of uEPSPs are ~ 0.5 mV (I tested this in VPM_activity_fitting.mlx)

    % ^ actually, I did the opposite.. uEPSP(TC) = 2 x uEPSP(CC)

    %%%%%%%%%%%%%%%%%%
    % neuron model parameter are taken from:
    % Pospischil, Martin, Maria Toledo-Rodriguez, Cyril Monier, Zuzanna Piwkowska, Thierry Bal, Yves Frégnac, Henry Markram, 
    % and Alain Destexhe. “Minimal Hodgkin–Huxley Type Models for Different Classes of Cortical and Thalamic Neurons.” 
    % Biological Cybernetics 99, no. 4 (November 1, 2008): 427–41. https://doi.org/10.1007/s00422-008-0263-8.

    %%%%%%%%%%%%%%%%%%

    % synaptic strength
    ss=getfield(synParams,'synstr'); % synaptic strength
    synstr_e = ss(1);
    synstr_i = ss(2);

    % probability of connectivity
    pc=getfield(synParams,'pconn'); % probability of connectivity
    pconn_e = pc(1);
    pconn_i = pc(2);

    % % % Adjusted Pospischil2008
    rsparam = {...
        'Cm', 1,...                      % uF/cm2
        'gleak',0.85*1e-5 * 1e3,...      % mS/cm2
        'Eleak',-69,...                  % mV
        'gNa',0.086 * 1e3,...            % mS/cm2
        'Na_V1',-64,...                  % mV
        'gK',0.006 * 1e3,...             % mS/cm2
        'K_V1',-64,...                   % mV
        'gM',18 * 1e-5 * 1e3,...         % mS/cm2
        'mM_tau_max',808,...             % ms
        };

    fsparam = {...
        'Cm', 1,...                     % uF/cm2
        'gleak',3.8*1e-5 * 1e3,...      % mS/cm2
        'Eleak',-70.4,...               % mV
        'gNa',0.078 * 1e3,...           % mS/cm2
        'Na_V1',-58,...                 % mV
        'gK',0.0039 * 1e3,...           % mS/cm2
        'K_V1',-58,...                  % mV
        'gM',5 * 1e-5 * 1e3,...         % mS/cm2
        'mM_tau_max',602,...            % ms
        };

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

    l2EEnetcon=ConnMx(ActSyn('l2ee',sparsify),'intrahom',vary_conn,l2e_num,l2e_num,ConnMxParam('l2ee','pconn', pconn_e),ConnMxParam('l2ee','nb', synstr_e),sparse_fact,ctx_seed); % l2E to l2E % removed autaptic
    l2EInetcon=ConnMx(ActSyn('l2ei',sparsify),'intrahet',vary_conn,l2e_num,l2i_num,ConnMxParam('l2ei','pconn', pconn_e),ConnMxParam('l2ei','nb', synstr_e),sparse_fact,ctx_seed); % l2E to l2I
    l2IEnetcon=ConnMx(ActSyn('l2ie',sparsify),'intrahet',vary_conn,l2i_num,l2e_num,ConnMxParam('l2ie','pconn', pconn_i),ConnMxParam('l2ie','nb', synstr_i),sparse_fact,ctx_seed); % l2I to l2E 
    l2IInetcon=ConnMx(ActSyn('l2ii',sparsify),'intrahom',vary_conn,l2i_num,l2i_num,ConnMxParam('l2ii','pconn', pconn_i),ConnMxParam('l2ii','nb', synstr_i),sparse_fact,ctx_seed); % l2I to l2I % removed autaptic

    % connect populations in L2 layer 
    L2.connections(1).direction='L2E->L2E'; % E to E
    L2.connections(end).mechanism_list={'iAMPA_temp'};
    L2.connections(end).parameters={'gAMPA',gAMPA,'netcon',l2EEnetcon,'IC',IC};

    L2.connections(end+1).direction='L2E->L2I'; % E to I
    L2.connections(end).mechanism_list={'iAMPA_temp'};
    L2.connections(end).parameters={'gAMPA',gAMPA,'netcon',l2EInetcon,'IC',IC}; 

    L2.connections(end+1).direction='L2I->L2E'; % I to E
    L2.connections(end).mechanism_list={'iGABAa_temp'};
    L2.connections(end).parameters={'gGABA',gGABAa,'netcon',l2IEnetcon,'IC',IC,'tauD',tauD_GABAa, 'EGABAa',EGABAa}; 

    L2.connections(end+1).direction='L2I->L2I'; % I to I
    L2.connections(end).mechanism_list={'iGABAa_temp'};
    L2.connections(end).parameters={'gGABA',gGABAa,'netcon',l2IInetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 
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
    l3EEnetcon=ConnMx(ActSyn('l3ee',sparsify),'intrahom',vary_conn,l3e_num,l3e_num,ConnMxParam('l3ee','pconn', pconn_e),ConnMxParam('l3ee','nb', synstr_e),sparse_fact,ctx_seed); % l2E to l2E % removed autoptic
    l3EInetcon=ConnMx(ActSyn('l3ei',sparsify),'intrahet',vary_conn,l3e_num,l3i_num,ConnMxParam('l3ei','pconn', pconn_e),ConnMxParam('l3ei','nb', synstr_e),sparse_fact,ctx_seed); % l2E to l2I
    l3IEnetcon=ConnMx(ActSyn('l3ie',sparsify),'intrahet',vary_conn,l3i_num,l3e_num,ConnMxParam('l3ie','pconn', pconn_i),ConnMxParam('l3ie','nb', synstr_i),sparse_fact,ctx_seed); % l2I to l2E 
    l3IInetcon=ConnMx(ActSyn('l3ii',sparsify),'intrahom',vary_conn,l3i_num,l3i_num,ConnMxParam('l3ii','pconn', pconn_i),ConnMxParam('l3ii','nb', synstr_i),sparse_fact,ctx_seed); % l2I to l2I % removed autoptic

    % connect populations in L3 layer 
    L3.connections(1).direction='L3E->L3E'; % E to E
    L3.connections(end).mechanism_list={'iAMPA_temp'};
    L3.connections(end).parameters={'gAMPA',gAMPA,'netcon',l3EEnetcon,'IC',IC};

    L3.connections(end+1).direction='L3E->L3I'; % E to I
    L3.connections(end).mechanism_list={'iAMPA_temp'};
    L3.connections(end).parameters={'gAMPA',gAMPA,'netcon',l3EInetcon,'IC',IC}; 

    L3.connections(end+1).direction='L3I->L3E'; % I to E
    L3.connections(end).mechanism_list={'iGABAa_temp'};
    L3.connections(end).parameters={'gGABA',gGABAa,'netcon',l3IEnetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 

    L3.connections(end+1).direction='L3I->L3I'; % I to I
    L3.connections(end).mechanism_list={'iGABAa_temp'};
    L3.connections(end).parameters={'gGABA',gGABAa,'netcon',l3IInetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 
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
    l4EEnetcon=ConnMx(ActSyn('l4ee',sparsify),'intrahom',vary_conn,l4e_num,l4e_num,ConnMxParam('l4ee','pconn', pconn_e),ConnMxParam('l4ee','nb', synstr_e),sparse_fact,ctx_seed); 
    l4EInetcon=ConnMx(ActSyn('l4ei',sparsify),'intrahet',vary_conn,l4e_num,l4i_num,ConnMxParam('l4ei','pconn', pconn_e),ConnMxParam('l4ei','nb', synstr_e),sparse_fact,ctx_seed); 
    l4IEnetcon=ConnMx(ActSyn('l4ie',sparsify),'intrahet',vary_conn,l4i_num,l4e_num,ConnMxParam('l4ie','pconn', pconn_i),ConnMxParam('l4ie','nb', synstr_i),sparse_fact,ctx_seed); 
    l4IInetcon=ConnMx(ActSyn('l4ii',sparsify),'intrahom',vary_conn,l4i_num,l4i_num,ConnMxParam('l4ii','pconn', pconn_i),ConnMxParam('l4ii','nb', synstr_i),sparse_fact,ctx_seed); 


    % connect populations in L4 layer 
    L4.connections(1).direction='L4E->L4E'; % E to E
    L4.connections(end).mechanism_list={'iAMPA_temp'};
    L4.connections(end).parameters={'gAMPA',gAMPA,'netcon',l4EEnetcon,'IC',IC};

    L4.connections(end+1).direction='L4E->L4I'; % E to I
    L4.connections(end).mechanism_list={'iAMPA_temp'};
    L4.connections(end).parameters={'gAMPA',gAMPA,'netcon',l4EInetcon,'IC',IC}; 

    L4.connections(end+1).direction='L4I->L4E'; % I to E
    L4.connections(end).mechanism_list={'iGABAa_temp'};
    L4.connections(end).parameters={'gGABA',gGABAa,'netcon',l4IEnetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 

    L4.connections(end+1).direction='L4I->L4I'; % I to I
    L4.connections(end).mechanism_list={'iGABAa_temp'};
    L4.connections(end).parameters={'gGABA',gGABAa,'netcon',l4IInetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 
    L4=dsCheckSpecification(L4); 

    % MAKE A THALAMIC LAYER:
        % there are no Inhibitory population in Th 
        % this population presents neurons of VPM thalamic nucleus

    input_amp_scaling = 0.08;
    try
        if size(ang,2) ~= VPM_num % data has to be: [timepoints, neurons]
            ang = ang';
            curv = curv';
        end    
        endSim = length(ang(:,1))*simulation_dt; 
    catch ME
        VPM_num = 70;
        if nargin < 9
            endSim = 1000;
        end
    end

    if strcmp(input_mode, 'exp')
        VPM.populations(1).name='VPM';
        VPM.populations(end).size=VPM_num;
        VPM.populations(end).equations={'dV/dt=@current+Iapp'};
        VPM.populations(end).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                           'CaBufferChing2010TC','iTChing2010TC',...
                                           'iHChing2010TC','iLeakChing2010TC',...
                                           'iKLeakChing2010TC', 'WhiskerData'};
        VPM.populations(end).parameters={'gH',0.005,'EndTime',endSim,'Iapp',1.5,...'scale_a',angle_contribution,
                            'signal_scale',input_amp_scaling,...
                            'ConvAng',ang,'ConvCurv',curv,...
                            };
        VPM=dsCheckSpecification(VPM); 

        % make 4-layered network
        Network=dsCombineSpecifications(L2,L3,L4,VPM);

        % make connectivity matrices for inter-layer connections: netcon=ConnMx(activity,layer,vary,Npre,Npost,va,Nb)
        % Th -> L4
        thl4EEnetcon=ConnMx(ActSyn('thl4e',sparsify),'trans',vary_conn,VPM_num,l4e_num,ConnMxParam('thl4e','pconn', pconn_e),ConnMxParam('thl4e','nb', pconn_e),sparse_fact,ctx_seed); % corr
        thl4EInetcon=ConnMx(ActSyn('thl4i',sparsify),'trans',vary_conn,VPM_num,l4i_num,ConnMxParam('thl4i','pconn', pconn_e),ConnMxParam('thl4i','nb', pconn_e),sparse_fact,ctx_seed); % corr

        % add feedforward connections (th->l4; l4->l3; l3->l2)
            % thalamic to l4
        Network.connections(end+1).direction='VPM->L4E'; % ThE to l4E
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EEnetcon, 'IC', IC};

        Network.connections(end+1).direction='VPM->L4I'; % ThE to l4I
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EInetcon, 'IC', IC}; 

    else
        % thalamic params
        CorrPoissR =    9;   % in Hz, firing rate of correlated poisson spike train // 7
        CorrPoissIn =   0.3; % input probability to thalamic cells // 0.2
        UncorrPoissR =  5;  % in Hz, firing rate of uncorrelated poisson spike train // 5
        UncorrPoissIn = 0.1; % input probability to thalamic cells // 0.1
        gSyn= 0.5; % synaptic conductance for the input

        % group thalamic variables
        vpmParams=[CorrPoissR,CorrPoissIn,UncorrPoissR,UncorrPoissIn,gSyn];
        vpmNum = 70;

        corr_contribution = 0.65;
        uncorr_contribution = 0.35;
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
                                        'T',endSim,'jitter_stddev_corr',5/simulation_dt,'tau_i_uncorr',25};

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
                                        'T',endSim,'jitter_stddev_corr',5/simulation_dt,'tau_i_corr',14};

        VPMcorr=dsCheckSpecification(VPMcorr); 
        VPMuncorr=dsCheckSpecification(VPMuncorr); 

        % make 4-layered network
        Network=dsCombineSpecifications(L2,L3,L4,VPMcorr,VPMuncorr);

        % make connectivity matrices for inter-layer connections: netcon=ConnMx(activity,layer,vary,Npre,Npost,va,Nb)
            % Th -> L4
        thl4EEnetcon2=ConnMx(ActSyn('thl4e',sparsify),'trans',vary_conn,vpmNum2,l4e_num,ConnMxParam('thl4e','pconn', pconn_e)*uncorr_contribution,ConnMxParam('thl4e','nb', pconn_e),sparse_fact,ctx_seed); % uncorr
        thl4EInetcon2=ConnMx(ActSyn('thl4i',sparsify),'trans',vary_conn,vpmNum2,l4i_num,ConnMxParam('thl4i','pconn', pconn_e)*uncorr_contribution,ConnMxParam('thl4i','nb', pconn_e),sparse_fact,ctx_seed); % uncorr
        thl4EEnetcon1=ConnMx(ActSyn('thl4e',sparsify),'trans',vary_conn,vpmNum1,l4e_num,ConnMxParam('thl4e','pconn', pconn_e)*corr_contribution,ConnMxParam('thl4e','nb', pconn_e),sparse_fact,ctx_seed); % corr
        thl4EInetcon1=ConnMx(ActSyn('thl4i',sparsify),'trans',vary_conn,vpmNum1,l4i_num,ConnMxParam('thl4i','pconn', pconn_e)*corr_contribution,ConnMxParam('thl4i','nb', pconn_e),sparse_fact,ctx_seed); % corr

        %-----------------------------------------------------------------------------------------------------------
        % add feedforward connections (th->l4; l4->l3; l3->l2)
            % thalamic to l4
        Network.connections(end+1).direction='VPMcorr->L4E'; % ThE to l4E
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EEnetcon1,'IC',IC};

        Network.connections(end+1).direction='VPMcorr->L4I'; % ThE to l4I
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EInetcon1,'IC',IC}; 

        Network.connections(end+1).direction='VPMuncorr->L4E'; % ThE to l4E
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EEnetcon2,'IC',IC};

        Network.connections(end+1).direction='VPMuncorr->L4I'; % ThE to l4I
        Network.connections(end).mechanism_list={'iAMPA_temp'};
        Network.connections(end).parameters={'gAMPA',gAMPA_th,'netcon',thl4EInetcon2,'IC',IC};

        %     % thalamic to l3
        % thl3EEnetcon=ConnMx(0,'trans','none',vpmNum1,l3e_num,(ConnMxParam('thl4e','pconn', pconn_e)*corr_contribution)/2,ConnMxParam('thl4e','nb', pconn_e)/2,1); % corr
        % thl3EInetcon=ConnMx(0,'trans','none',vpmNum1,l3i_num,(ConnMxParam('thl4i','pconn', pconn_e)*corr_contribution)/2,ConnMxParam('thl4i','nb', pconn_e)/2,1); % corr    
        % 
        % Network.connections(end+1).direction='VPMcorr->L3E'; % ThE to l4E
        % Network.connections(end).mechanism_list={'iAMPA_temp'};
        % Network.connections(end).parameters={'gAMPA',gAMPA*AMPA_scale_th,'netcon',thl3EEnetcon,'IC',IC};
        %  
        % Network.connections(end+1).direction='VPMcorr->L3I'; % ThE to l4I
        % Network.connections(end).mechanism_list={'iAMPA_temp'};
        % Network.connections(end).parameters={'gAMPA',gAMPA*AMPA_scale_th,'netcon',thl3EInetcon,'IC',IC}; 
    end

    %     % L4 -> L3
    l4l3EEnetcon=ConnMx(ActSyn('l4el3e',sparsify),'trans',vary_conn,l4e_num,l3e_num,ConnMxParam('l4el3e','pconn', pconn_e),ConnMxParam('l4el3e','nb', synstr_e),sparse_fact,ctx_seed);
    l4l3EInetcon=ConnMx(ActSyn('l4el3i',sparsify),'trans',vary_conn,l4e_num,l3i_num,ConnMxParam('l4el3i','pconn', pconn_e),ConnMxParam('l4el3i','nb', synstr_e),sparse_fact,ctx_seed);
    l4l3IEnetcon=ConnMx(ActSyn('l4il3e',sparsify),'trans',vary_conn,l4i_num,l3e_num,ConnMxParam('l4il3e','pconn', pconn_i),ConnMxParam('l4il3e','nb', synstr_i),sparse_fact,ctx_seed);
        % L3 -> L2
    l3l2EEnetcon=ConnMx(ActSyn('l3el2e',sparsify),'trans',vary_conn,l3e_num,l2e_num,ConnMxParam('l3el2e','pconn', pconn_e),ConnMxParam('l3el2e','nb', synstr_e),sparse_fact,ctx_seed);
    l3l2EInetcon=ConnMx(ActSyn('l3el2i',sparsify),'trans',vary_conn,l3e_num,l2i_num,ConnMxParam('l3el2i','pconn', pconn_e),ConnMxParam('l3el2i','nb', synstr_e),sparse_fact,ctx_seed);
    l3l2IEnetcon=ConnMx(ActSyn('l3il2e',sparsify),'trans',vary_conn,l3i_num,l2e_num,ConnMxParam('l3il2e','pconn', pconn_i),ConnMxParam('l3il2e','nb', synstr_i),sparse_fact,ctx_seed);
    l3l2IInetcon=ConnMx(ActSyn('l3il2i',sparsify),'trans',vary_conn,l3i_num,l2i_num,ConnMxParam('l3il2i','pconn', pconn_i),ConnMxParam('l3il2i','nb', synstr_i),sparse_fact,ctx_seed);

        % l4 to l23 (l4 to l3)
    Network.connections(end+1).direction='L4E->L3E'; % l4E to l3E
    Network.connections(end).mechanism_list={'iAMPA_temp'};
    Network.connections(end).parameters={'gAMPA',gAMPA,'netcon',l4l3EEnetcon,'IC',IC};

    Network.connections(end+1).direction='L4E->L3I'; % l4E to l3I
    Network.connections(end).mechanism_list={'iAMPA_temp'};
    Network.connections(end).parameters={'gAMPA',gAMPA,'netcon',l4l3EInetcon,'IC',IC};

    Network.connections(end+1).direction='L4I->L3E'; % l4I to l3E
    Network.connections(end).mechanism_list={'iGABAa_temp'};
    Network.connections(end).parameters={'gGABA',gGABAa,'netcon',l4l3IEnetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 

        % l23 to l23 feedforward (l3 to l2)
    Network.connections(end+1).direction='L3E->L2E'; % l3E to l2E
    Network.connections(end).mechanism_list={'iAMPA_temp'};
    Network.connections(end).parameters={'gAMPA',gAMPA,'netcon',l3l2EEnetcon,'IC',IC}; 

    Network.connections(end+1).direction='L3E->L2I'; % l3E to l2I
    Network.connections(end).mechanism_list={'iAMPA_temp'};
    Network.connections(end).parameters={'gAMPA',gAMPA,'netcon',l3l2EInetcon,'IC',IC};

    Network.connections(end+1).direction='L3I->L2E'; % l3I to l2E
    Network.connections(end).mechanism_list={'iGABAa_temp'};
    Network.connections(end).parameters={'gGABA',gGABAa,'netcon',l3l2IEnetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 

    Network.connections(end+1).direction='L3I->L2I'; % l3I to l2I
    Network.connections(end).mechanism_list={'iGABAa_temp'};
    Network.connections(end).parameters={'gGABA',gGABAa,'netcon',l3l2IInetcon,'IC',IC,'tauD',tauD_GABAa,'EGABAa',EGABAa}; 

    Network=dsCheckSpecification(Network); 



    %SIMULATION:
    tspan = [0+simulation_dt endSim];     % time span [begin end], ms
%     ms_to_s = 1e-3;
%     disp(['Simulation time: ' num2str(round(length([0:simulation_dt:endSim])*simulation_dt*ms_to_s,3)) 's. (seed ' num2str(ctx_seed) ')'])

    % this part controls if random processes during the simulation (e.g., the thalamic activity) are truly random or fixed-random (with a specific seed) 
    if ~isempty(vpm_seed) % fixed random
        net = dsSimulate(Network,'tspan',tspan,'dt',simulation_dt, 'random_seed', vpm_seed);
    else
        rng('shuffle') % truly random
        net = dsSimulate(Network,'tspan',tspan,'dt',simulation_dt);
    end

end