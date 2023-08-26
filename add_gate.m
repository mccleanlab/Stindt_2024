function data_out = add_gate(data_in, gates,x)

% [DATA_OUT] = ADDGATE(DATA_IN,GATE) applys a predefined set of gates to
% tabular flow cytometry measurements DATA_IN to return a table DATA_OUT
% with labels for each event indicating whether they lie within a given
% gate. It appends a logical column for each gate where each event is
% labeled true if it falls within the gate and false if it does not. It
% also appends a logical column GATE_NET indicating whether events lie
% within all applied gates.
%
% This function does not discard events that fall outside any gates, it
% simply labels events with logical labels that can be applied when
% plotting or to select subsets of events. The set of gates can be drawn
% using drawGate() of autoGate() or loaded from .mat files if they were
% previously drawn using either of these functions.

% data_in = load_fcs('map','plate','folder',folder_list{1});
% gates = [BactGate; BactDbls];

% for i = 1:size(gates,1)
%     gatenames{i,1} = cellstr(inputname(2));
% end
gatename = inputname(2);

data_out = data_in;

for f = 1:numel(data_in)
    
    fcsdat = data_out(f).fcsdat;
    fcshdr = data_out(f).fcshdr;
    
    for nGate=1:size(gates,1)
        % Reset variables
        xdata = [];
        ydata = [];
        idx = [];
        idx0=[];
        
        channels2gate = gates{nGate,1};
        gate_name = strrep(channels2gate,'-','');
%         gate_name = strcat('gate_',gate_name{1},'_',gate_name{2});
%         gate_net = 'gate_net';
        gate_name = strcat(gatename,'_',gate_name{1},'_',gate_name{2});
        gate_net = strcat(gatename,'_net');
        
        cX = find(strcmp({fcshdr.par.name},channels2gate(1))==1);
        cY = find(strcmp({fcshdr.par.name},channels2gate(2))==1);
        
        channels2scale = gates{nGate,2};
        if strcmp(channels2scale(1),'log')
            xdata = fcsdat{:,cX};
            xdata(xdata<=0) = nan;
            xdata = log10(xdata);
        else
            xdata = fcsdat{:,cX};
        end
        
        if strcmp(channels2scale(2),'log')
            ydata = fcsdat{:,cY};
            ydata(ydata<=0) = nan;
            ydata = log10(ydata);
        else
            ydata = fcsdat{:,cY};
        end
        
        for i = 1:numel(gates{nGate,3})
            gatePts = gates{nGate,3}{1,i};
            if size(gates(nGate,:),2)==5
                if strcmp(gates{nGate,5},'Subtract')
                    idx0(:,i) = i*double(~inpolygon(xdata, ydata, gatePts(:,1), gatePts(:,2)));
                else
                    idx0(:,i) = i*double(inpolygon(xdata, ydata, gatePts(:,1), gatePts(:,2)));
                end
            else
                idx0(:,i) = i*double(inpolygon(xdata, ydata, gatePts(:,1), gatePts(:,2)));
            end

        end
        
        idx = max(idx0,[],2);
        idx = logical(idx);
        
        gatenamelist{nGate} = gate_name;
        gateidxlist{nGate} = idx;
        data_out(f).fcsdat.(gate_name) = idx;
    end
    
    idxnet = all([gateidxlist{:,:}],2);
    data_out(f).fcsdat.(gate_net) = idxnet;
end



