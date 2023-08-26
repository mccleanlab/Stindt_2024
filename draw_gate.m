function gate_out = draw_gate(channels_to_gate, channels_to_scale)

% [GATEOUT] = DRAW_GATE(CHANNELS_TO_GATE,CHANNELS_TO_SCALE) allows one to
% draw gates on the channel pairs stored in CHANNELS_TO_GATE using the
% scales (linear/log) set in CHANNELS_TO_SCALE. Any arbitrary number of
% channel pairs can be used for gating. The .fcs file used to define the
% gates is selected by UI prompt. The gates stored in the cell array
% GATE_OUT, which is also exported as a .mat file for easy reuse. This
% function does not apply the gate to your measurements. That must be done
% later using the function add_gate.m

% EXAMPLE: 
% channels_to_gate = {'FSC-A', 'SSC-A'; 'FSC-A', 'FSC-H'}
% channels_to_scale = {'linear','linear';'log','linear'};
% gate_out = draw_gate(channels_to_gate,channels_to_scale) 

% The above code will allow you to draw a gate first on a plot of FSCA-A vs
% SCC-A (with linear scaling) then on a plot of FSC-A vs FSC-H (where
% FSC-A, arbitrarily, is log scaled and FSC-H is linear scaled).

% Note: when drawing gates press ENTER to close a polygon. Once the polygon is
% closed you may edit the gate further by moving vertices or
% adding/deleting vertices via the right-click menu. Once you are happy
% with your gate, press enter to show the points within the gate and ENTER
% again to save the gate. Continue this process until all channel pairs are
% gated. 

% The following is only to troubleshoot function script, keep commented
% out, generally...
% clearvars -except b
% 
% channels_to_gate = {'FSC-A', 'SSC-A','Include'}; % Specify pairs of channels for gating
% channels_to_scale = {'hyperlog','hyperlog'}; % Specify scale for each pair of channels


% Select .fcs file on which to draw gate
[files, folder] = uigetfile('.fcs','Select a .fcs file to gate');

% Correct filename format if only one file loaded
if ischar(files)==1
    files = {files};
else
    files = files';
end

% Load data from .fcs file
[fcsdat, fcshdr] = fca_readfcs([folder files{1}]);

for n_gate = 1:size(channels_to_gate,1)
    
    % Get index of channels used to gate x and y axes from channel names
    cX = find(strcmp({fcshdr.par.name},channels_to_gate(n_gate,1))==1);
    cY = find(strcmp({fcshdr.par.name},channels_to_gate(n_gate,2))==1);
    
    % Get x measurments and scale if needed
    if strcmp(channels_to_scale(n_gate,1),'log')
        xdata = fcsdat(:,cX);
        xdata(xdata<=0) = nan;
        xdata = log10(xdata);
    elseif strcmp(channels_to_scale(n_gate,1),'arcsinh')
        xdata = fcsdat(:,cX);
        xdata = asinh(xdata);
    elseif strcmp(channels_to_scale(n_gate,1),'hyperlog')
        xdata = fcsdat(:,cX);
        maxx = max(abs(xdata));
        d = ceil(abs(real(log(min(xdata)-max(xdata))))); %Number of (real) decades in data, rounded up
        for i = 1:length(xdata)
            if xdata(i)<0
                xdata(i) = real(sqrt(-10^(-d*xdata(i)/maxx) + (b*d*xdata(i)/maxx)+1));
            else
                xdata(i) = real(sqrt(10^(d*xdata(i)/maxx) + (b*d*xdata(i)/maxx)-1));
            end
        end
    else
        xdata = fcsdat(:,cX);
    end
    
    % Get y measurements and scale if needed
    if strcmp(channels_to_scale(n_gate,2),'log')
        ydata = fcsdat(:,cY);
        ydata(ydata<=0) = nan;
        ydata = log10(ydata);
    elseif strcmp(channels_to_scale(n_gate,2),'arcsinh')
        ydata = fcsdat(:,cY);
        ydata = asinh(ydata);
    elseif strcmp(channels_to_scale(n_gate,2),'hyperlog')
        ydata = fcsdat(:,cY);
        d = ceil(abs(real(log(min(ydata)-max(ydata))))); %Number of (real) decades in data, rounded up
        maxy = max(abs(ydata));
        for i = 1:length(ydata)
            if ydata(i)<0
                ydata(i) = real(sqrt(-10^(-d*ydata(i)/maxy) + (b*d*ydata(i)/maxy)+1));
            else
                ydata(i) = real(sqrt(10^(d*ydata(i)/maxy) + (b*d*ydata(i)/maxy)-1));
            end
        end
    else
        ydata = fcsdat(:,cY);
    end
    
    % Exclude events removed by previous gate
    if n_gate>1
        idxp = gate{n_gate-1,4};
        if strcmp(channels_to_gate(n_gate-1,3),'Subtract')
            xdata = xdata(~idxp);
            ydata = ydata(~idxp);
        else
            xdata = xdata(idxp);
            ydata = ydata(idxp);
        end
    end
    
    % Set number of bins based on data and flag events on edge
    n0 = length(xdata);
    edge_idx = any([xdata==max(xdata), ydata==max(ydata)],2);
    if n0 > 300
        nbins = round(n0/100);
    else
        nbins = n0;
    end
    
    % Plot measurements to be gated
    clear g
    g = gramm('x',xdata(~edge_idx),'y',ydata(~edge_idx));
    g.stat_bin2d('nbins',[nbins nbins],'geom','image');
    g.no_legend();
%     figure('Position',[100 100 800 800])
    g.set_names('x',channels_to_gate{1},'y',channels_to_gate{1});
%     g.axe_property('XLim',[0 10],'YLim',[0 10],'XScale','log','YScale','log')
    g.draw();
    
    % Set x and y limits
    xlim auto
    ylim auto
    
    % Delete unwanted GRAMM axis to enable drawing on top of plot
    ax = findall(gcf, 'type', 'axes');
    delete(ax(1))
    
    % Show events beyond plot edge
    hold on;
    scatter(xdata(edge_idx),ydata(edge_idx),'r.')
    xlabel([fcshdr.par(cX).name newline '(' num2str(100*sum(edge_idx)/numel(xdata)) '% of events beyond plot edges)']);
    ylabel(fcshdr.par(cY).name);
    
    % Draw gate over measurements
    gate_vertices = drawpolygon('Color',[255 94 105]./255);
    pause
    gate_vertices = gate_vertices.Position;
    gate_vertices = [gate_vertices; gate_vertices(1,:)];
    
    % Get index of events within gate
    idx = inpolygon(xdata, ydata, gate_vertices(:,1),gate_vertices(:,2));
    
    % Show events within gate on plot
    nf = length(xdata(idx));
    scatter(xdata,ydata,1,'filled','MarkerFaceColor',[255 94 105]./255) ;
    xlabel([fcshdr.par(cX).name newline '(' num2str(n0) ' events before gate, ' num2str(nf) ' events after gate)'])
    pause
    
    % Save gate
    gate{n_gate,1} = channels_to_gate(n_gate,:);
    gate{n_gate,2} = channels_to_scale(n_gate,:);
    gate{n_gate,3} = {gate_vertices};
    gate{n_gate,4} = idx;
    gate{n_gate,5} = {files(1)};
    if size(channels_to_gate(n_gate,:),2)==3 %Check that there's a field for inclusion/subtraction
        gate{n_gate,6} = channels_to_gate(n_gate,3);
    else
        gate{n_gate,6} = 'Include';
    end

    
    close all
end

% Get all gates and discard per gate indices (unneeded)
gate_out = gate;
gate_out(:,4) = [];

% Save gate
% [~, filename_out, ~] = fileparts(files{1});
% save(['gate_' filename_out '.mat'],'gate_out')