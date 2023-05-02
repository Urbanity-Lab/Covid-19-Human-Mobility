clear all;
close all;
MobilityData = readtable('2020_US_Region_Mobility_Report.csv');
whos MobilityData
MobilityData.Properties.VariableNames
idxFlorida = MobilityData.sub_region_1 == 'Florida' ;
FloridaMobility = MobilityData(idxFlorida,:)
FloridaMobility = FloridaMobility{:,:};
FloridaMobility = FloridaMobility(1,:);