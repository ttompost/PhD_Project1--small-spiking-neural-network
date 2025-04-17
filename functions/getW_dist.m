function getW_dist(model)

    names = fieldnames(model.model.fixed_variables); 
    netcons = cellfun(@(x) endsWith(x, '_netcon'), names, 'UniformOutput',0);

    netcon_idx = find(cell2mat(netcons));

    w_data = [];
    for ii = 1:length(netcon_idx)
        idx = netcon_idx(ii);
        thisNetcon = model.model.fixed_variables.(names{idx});

        w_data(end+1,1) = max(unique(thisNetcon));

        if length(unique(thisNetcon)) > 2
            disp('More than 1 W found.')
            break
        end
    end
    histogram(w_data, 0:0.2:50)
end