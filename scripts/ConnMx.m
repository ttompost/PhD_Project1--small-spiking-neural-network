function Mx=ConnMx(activity, layer, vary, Npre, Npost, Pconn, Nb, sf, seed)
% this function calculates the connectivity matrix with various input
% conditions:
    % (1) activity:
            % true/1:  scaling factor is applied
            % false/0: scaling factor is not applied
    % (2) layer: 
            % intrahom: projections INSIDE of the layer, same cell type; avoids autoptic connections
            % intrahet: projections INSIDE of the layer, different cell type
            % trans: projections BETWEEN the layers
    % (3) vary:
            % bouts: vary Nb while connectivity probability (Pconn) is fixed
            % probs: vary connectivity probability (Pconn) while Nb is fixed
    % (4) Npre = number of presynaptic neurons
    % (5) Npost = number of postsynaptic neuron
    % (6) Pconn = probability of connection (max=1)
    % (7) Nb = number of boutons
    % (8) sf = scaling factor 
    
    
%%%%% precondition:
% when I test the behaviour of the network with n(cells)=1 in each layer, I need:
%   #1 Ne->Ne: (intra-hom) = 0; (trans) = 1;
%   #2 Ne->Ni: (intra-het) = 1; (trans) = 1;
%   #3 Ni->Ne: (intra-het) = 1; (trans) = 1;
%   #4 Ni->Ni: (intra-hom) = 0; (trans) = 1;

% that is why I have 2 'intra' cases:
%   #1 hom: homogeneous/same celly type (e.g. Ne->Ne)
%   #2 het: heterogeneous/different cell type (e.g.Ne->Ni)
   
%% choose params
m = Npre; % desired row size
n = Npost; % desired column size

if ~exist('seed','var') || isempty(seed)
    rng('shuffle')
else
    rng(seed);
%     fprintf('Using seed %i for ConnMx.\n',seed);
end

switch activity 
    case 1 % apply scaling factor
    case 0 % do not apply scaling factor
        sf=1;
    otherwise
        disp('Invalid input for "activity" in ConnMx. Can not determine scaling fact.')
end%activity 

switch vary
    case 'probs' %vary connectivity probability
        Pconn=Pconn*sf;
    case 'bouts' % vary bouton number (i.e. synaptic strength)
        Nb=Nb*sf;
    otherwise
%         disp('Invalid input for "vary" in ConnMx. Probs or Bouton num. are not varied.')
end%vary

if Nb == 0; disp('Connx weight is zero, meaning nothing is connected.'); end

% make the connectivity matrix
f = abs(Pconn); % fraction of Nb's to place randomly in the result
Mx = zeros(m,n); % pre-allocate result
k = round(f*m*n); % number of Nb's to place in result
Mx(randperm(m*n,k)) = Nb; % randomly put Nb's in result using linear indexing

% additional conditions
switch layer
    case 'intrahom' % this one avoids autoptic connections
          if m==n
              for ii=1:m
                  Mx(ii,ii)=0; % avoids autoptic connections
              end
          end
    case 'intrahet' 
    case 'trans'
    otherwise
        disp('Invalid input for "layer" in ConnMx. Autaptic connx are present.')
end%layer

if isempty(Mx); disp('Connx Mx is empty, exiting simulation..'); end
end%end function

