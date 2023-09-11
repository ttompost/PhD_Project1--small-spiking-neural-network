function ParamValue=ConnMxParam(projection, type, scale)
% function outputs value for either the number of boutons (nb) or connectivity
% probability (pconn) for a specific projection

    switch type
        case 'pconn'
            fval=1;
            pconnstr=scale;
            synstr=1;
        case 'nb'
            fval=2;
            synstr=scale;
            pconnstr=1;
        otherwise
            disp('Error in providing input variables in ConnMxParam')
    end

paramset = 'new';
% 'old'
% old values come from independent literature review; sources reported in
% my master's thesis. they were fitted to optimize network activity, so
% they have to be used with: (uncomment to use)
% synstr=0.10; % this is the reduction factor of the synaptic strength

% 'new'
% new values come from the blue brain database (recommended)
% source: https://bbp.epfl.ch/nmc-portal/welcome.html
  
p.old.l2ee=[0.10*pconnstr, 3*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.old.l2ei=[0.65*pconnstr, 10*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.old.l2ie=[0.60*pconnstr, 6*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.old.l2ii=[0.55*pconnstr, 10*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

p.new.l2ee=[0.056*pconnstr, 2.8*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.new.l2ei=[0.051*pconnstr, 8.1*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.new.l2ie=[0.059*pconnstr, 17*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.new.l2ii=[0.051*pconnstr, 15*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

p.old.l3ee=[0.10*pconnstr, 3*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.old.l3ei=[0.65*pconnstr, 10*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.old.l3ie=[0.60*pconnstr, 6*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.old.l3ii=[0.55*pconnstr, 10*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

p.new.l3ee=[0.056*pconnstr, 2.8*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.new.l3ei=[0.051*pconnstr, 8.1*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.new.l3ie=[0.059*pconnstr, 17*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.new.l3ii=[0.051*pconnstr, 15*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

p.old.l4ee=[0.06*pconnstr, 3*synstr]; % old: [0.06*pconnstr, 3*synstr] / new: [0.076*pconnstr, normrnd(3.3,1.4)*synstr]
p.old.l4ei=[0.43*pconnstr, 8*synstr]; % old: [0.43*pconnstr, 8*synstr] / new: [0.042*pconnstr, normrnd(7.9,3)*synstr]
p.old.l4ie=[0.44*pconnstr, 14*synstr]; % old: [0.44*pconnstr, 14*synstr] / new: [0.063*pconnstr, normrnd(16,6.2)*synstr]
p.old.l4ii=[0.55*pconnstr, 10*synstr]; % old: [0.55*pconnstr, 10*synstr] / new: [0.062*pconnstr, normrnd(14,6)*synstr]

p.new.l4ee=[0.076*pconnstr, 3.3*synstr]; % old: [0.06*pconnstr, 3*synstr] / new: [0.076*pconnstr, normrnd(3.3,1.4)*synstr]
p.new.l4ei=[0.042*pconnstr, 7.9*synstr]; % old: [0.43*pconnstr, 8*synstr] / new: [0.042*pconnstr, normrnd(7.9,3)*synstr]
p.new.l4ie=[0.063*pconnstr, 16*synstr]; % old: [0.44*pconnstr, 14*synstr] / new: [0.063*pconnstr, normrnd(16,6.2)*synstr]
p.new.l4ii=[0.062*pconnstr, 14*synstr]; % old: [0.55*pconnstr, 10*synstr] / new: [0.062*pconnstr, normrnd(14,6)*synstr]

p.old.thl4e=[0.43*pconnstr, 7*synstr]; 
p.old.thl4i=[0.5*pconnstr, 9*synstr]; 

p.new.thl4e=[(0.43/10)*pconnstr, 7*synstr];   % no data from new database so i divide pconn with 10
p.new.thl4i=[(0.5/10)*pconnstr, 9*synstr];    % no data from new database so i divide pconn with 10

p.old.l4el3e=[0.12*pconnstr, 5*synstr]; % old: [0.12*pconnstr, 5*synstr] / new: [0.0058*pconnstr, normrnd(2.4,0.78)*synstr]
p.old.l4el3i=[0.2*pconnstr, 3*synstr]; % old: [0.20*pconnstr, 3*synstr] / new: [0.0033*pconnstr, normrnd(5.7,1.9)*synstr]
p.old.l4il3e=[0.44*pconnstr, 15*synstr]; % old: [0.44*pconnstr, 15*synstr] / new: [0.017*pconnstr, normrnd(14,5.5)*synstr]
% p.old.l4il3i=[0.44*pconnstr, 15*synstr]; % new: [0.014*pconnstr, normrnd(15,6.9)*synstr]  % does not exist tho?

p.new.l4el3e=[0.0058*pconnstr, 2.4*synstr]; % old: [0.12*pconnstr, 5*synstr] / new: [0.0058*pconnstr, normrnd(2.4,0.78)*synstr]
p.new.l4el3i=[0.0033*pconnstr, 5.7*synstr]; % old: [0.20*pconnstr, 3*synstr] / new: [0.0033*pconnstr, normrnd(5.7,1.9)*synstr]
p.new.l4il3e=[0.017*pconnstr, 14*synstr]; % old: [0.44*pconnstr, 15*synstr] / new: [0.017*pconnstr, normrnd(14,5.5)*synstr]
% p.new.l4il3i=[0.014*pconnstr, 15*synstr]; % new: [0.014*pconnstr, normrnd(15,6.9)*synstr]  % does not exist tho?

p.old.l4el2e=[(0.12/2)*pconnstr, (5/2)*synstr]; 
p.old.l4el2i=[(0.2/2)*pconnstr, (3/2)*synstr]; 
p.old.l4il2e=[(0.44/2)*pconnstr, (15/2)*synstr]; 

p.new.l4el2e=[0.0058/2*pconnstr, (2.4/2)*synstr]; 
p.new.l4el2i=[0.0033/2*pconnstr, (5.7/2)*synstr]; 
p.new.l4il2e=[0.017/2*pconnstr, (14/2)*synstr]; 

p.old.l3el2e=[0.10*pconnstr, 3*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.old.l3el2i=[0.65*pconnstr, 10*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.old.l3il2e=[0.60*pconnstr, 6*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.old.l3il2i=[0.55*pconnstr, 10*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

p.new.l3el2e=[0.056*pconnstr, 2.8*synstr]; % old: [0.10*pconnstr, 3*synstr] / new: [0.056*pconnstr, normrnd(2.8,1.2)*synstr]
p.new.l3el2i=[0.051*pconnstr, 8.1*synstr]; % old: [0.65*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(8.1,3.1)*synstr]
p.new.l3il2e=[0.059*pconnstr, 17*synstr]; % old: [0.60*pconnstr, 6*synstr] / new: [0.059*pconnstr, normrnd(17,6.5)*synstr]
p.new.l3il2i=[0.051*pconnstr, 15*synstr];% old: [0.55*pconnstr, 10*synstr] / new: [0.051*pconnstr, normrnd(15,7)*synstr]

    ParamValue=p.(paramset).(projection)(fval);
end