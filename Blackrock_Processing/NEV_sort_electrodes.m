function NEV_out=NEV_sort_electrodes(NEV)
% NEV_out=NEV_sort_electrodes(NEV)
% Adds another field called NEV.ElectrodeSpikes to NEV.Data, which takes the spike data from
% NEV.Data.Spikes and sorts it into each electrode (which blackrock doesn't
% do for some reason)
% David Xing 2-22-2016

%initialize
NEV_out=NEV;

%I'm assuming that the neural recording electrodes are the
%electrodes labeled with 'elec#' in the ElectrodeLabel field
Electrodes=cat(2,NEV.ElectrodesInfo(find(cellfun(@(x) strcmpi(x(1:4)','elec'),{NEV.ElectrodesInfo.ElectrodeLabel}))).ElectrodeID);
nElectrodes=length(Electrodes);

for whichElectrode=1:nElectrodes
    %indices for each particular electrode
    ElecInd{whichElectrode}=find(NEV.Data.Spikes.Electrode==Electrodes(whichElectrode));
    
    %copy over fields from Spikes to ElectrodeSpikes
    NEV_out.Data.ElectrodeSpikes(whichElectrode).ElectrodeID=Electrodes(whichElectrode);
    NEV_out.Data.ElectrodeSpikes(whichElectrode).LabelNumber=str2num(NEV.ElectrodesInfo(whichElectrode).ElectrodeLabel(5:end)');
    NEV_out.Data.ElectrodeSpikes(whichElectrode).TimeStamp=NEV.Data.Spikes.TimeStamp(ElecInd{Electrodes(whichElectrode)});
    NEV_out.Data.ElectrodeSpikes(whichElectrode).Unit=NEV.Data.Spikes.Unit(ElecInd{Electrodes(whichElectrode)});
    NEV_out.Data.ElectrodeSpikes(whichElectrode).Waveform=NEV.Data.Spikes.Waveform(:,ElecInd{Electrodes(whichElectrode)});
    
end

%check sum
foundunits=sum(cellfun(@length,ElecInd));
allunits=length(NEV.Data.Spikes.TimeStamp);

if (foundunits~=allunits)
    warning(['Not all spike units accounted for by the electrodes! # Electrodes: ' num2str(nElectrodes),...
        ', Total # spike units: ' num2str(allunits) ', # Spike units accounted for by elecrodes: ' num2str(foundunits)])
end


% 
