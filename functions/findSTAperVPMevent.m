function events = findSTAperVPMevent(vpm,ctx, Bin, Spikes, Win, tstart, tstop)
   
    % get spike times
    VPMactivity = getSpikeTimes(vpm, [tstart:0.01:tstop]);
    CTXactivity = getSpikeTimes(ctx, [tstart:0.01:tstop]);

    % find psth and peaks
    [vpm_psth, vpm_edges] = plotPsthAndRaster('psth',VPMactivity,Bin,0,[]);
    [~, vpm_events_location] = findpeaks(vpm_psth,'minpeakheight',Spikes,'minpeakdistance',round(Win/Bin),'sortstr','descend');
    time_vec_vpm = vpm_edges(vpm_events_location); % shift some bins before
    
    [ctx_psth, ctx_edges] = plotPsthAndRaster('psth',CTXactivity,Bin,0,[]);
    pre = round(30/Bin);
    post = round(30/Bin);

    % find psth part per event
    events = [];
    for ev = 1:length(time_vec_vpm)
        try
            events{ev,1} = vpm_psth(nearest(vpm_edges, time_vec_vpm(ev)-pre):nearest(vpm_edges, time_vec_vpm(ev))+post);
            events{ev,2} = ctx_psth(nearest(vpm_edges, time_vec_vpm(ev)-pre):nearest(vpm_edges, time_vec_vpm(ev))+post);
        catch ME
%             events{ev,1} = vpm_psth(nearest(vpm_edges, time_vec_vpm(ev)-pre):nearest(vpm_edges, time_vec_vpm(ev))+post);
%             events{ev,2} = ctx_psth(nearest(vpm_edges, time_vec_vpm(ev)-pre):nearest(vpm_edges, time_vec_vpm(ev))+post);
        end
    end
end
