%odor structure maker
%this is a quick database of odors, mainly for organizing colors and molar
%concentration given the final dillution.
%Future version of the database will be an h5 file taken from a google
%spreadsheet.
%Rinberg lab, June 2014, ZK.
clear odor;
odor=struct('name',{},'color',{},'colorCode',{},'plotSymbol',{},'n',{},'ec50',{},'molarFactor',{},'neatMolarity',{},'density',{},'molecularWeight',{});
nOdor=0;


nOdor=nOdor+1;
odor(nOdor).name        = {'ethyl_tiglate','ethyl-tiglate','ethyl tiglate'};
odor(nOdor).color       = 'g';
odor(nOdor).colorCode   = [0 1 0];
odor(nOdor).plotSymbol      = '-^';  
odor(nOdor).n           = nOdor;
odor(nOdor).ec50        = 1.4;
odor(nOdor).neatMolarity= 3.68e-5; %Molarity of saturated vapor in headspace of vial with pure odor
odor(nOdor).density     = 0.92; %density (g/ml)
odor(nOdor).molecularWeight = 128.17; %molecular weight

odor(nOdor).molarFactor = 165.09e-6; % factor for concentration to molarity (WRONG WAY OF DOING THE CALCULATION:C in headspace in um)

nOdor=nOdor+1;
o.name        = {'2-hydroxyacetophenone'};
o.color       = 'b';
o.colorCode   = [0.7 0 1];
o.plotSymbol      = '-^';
o.n           = nOdor;
o.ec50        = 17.5;
o.neatMolarity= 8.03e-7; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.13; %density (g/ml)
o.molecularWeight = 136.15; %molecular weight

o.molarFactor = 3.78e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'methyl_salicylate','methyl salicylate'};
o.color       = 'y';
o.plotSymbol      = '-^';
o.colorCode   = [1 1 0];
o.n           = nOdor;
o.ec50        = 2.2;
o.neatMolarity= 3.51e-7; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.17; %density (g/ml)
o.molecularWeight = 152.15; %molecular weight

o.molarFactor = 1.87e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'acetophenone'};
o.color       = 'b';
o.colorCode   = [0 0 1];
o.plotSymbol      = '-^';
o.n           = nOdor;
o.ec50        = 17.5;
o.neatMolarity= 5.15e-6; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.03; %density (g/ml)
o.molecularWeight = 120.15; %molecular weight

o.molarFactor = 24.52e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'menthone'};
o.color       = 'b';
o.colorCode   = [0 0.7 1];
o.plotSymbol      = '-^';
o.n           = nOdor;
o.ec50        = nan;
o.neatMolarity= 3.74e-6; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = .9; %density (g/ml)
o.molecularWeight = 154.25; %molecular weight

o.molarFactor = 19.9e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'benzaldehyde'};
o.color       = 'y';
o.plotSymbol      = '-^';
o.colorCode   = [1 1 0];
o.n           = nOdor;
o.ec50        = nan;
o.neatMolarity= 1.86e-5; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.04; %density (g/ml)
o.molecularWeight = 106.12; %molecular weight

o.molarFactor = 68.3e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'2-4-dimethylacetophenone'};
o.color       = 'y';
o.plotSymbol  = '-o';
o.colorCode   = [0.3 1 0.7];
o.n           = nOdor;
o.ec50        = 2.2;
o.neatMolarity= 3.51e-7; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.17; %density (g/ml)
o.molecularWeight = 152.15; %molecular weight

o.molarFactor = 4.36e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'4-methyl-acetophenone','4-methyl_acetophenone','4-methylacetophenone'};
o.color       = 'm';
o.plotSymbol      = '-^';
o.colorCode   = [1 0 1];
o.n           = nOdor;
o.ec50        = 2.2;
o.neatMolarity= 2.57e-6; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.01; %density (g/ml)
o.molecularWeight = 134.18; %molecular weight

o.molarFactor = 11.89e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'pinene','a-pinene','a_pinene'};
o.color       = 'c';
o.plotSymbol      = ':v';
o.colorCode   = [0 1 1];
o.n           = nOdor;
o.ec50        = nan;
o.neatMolarity= nan; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = nan; %density (g/ml)
o.molecularWeight = nan; %molecular weight

o.molarFactor = 258.8e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'1-octanol'};
o.color       = ':r';
o.plotSymbol      = 'v';
o.colorCode   = [1 0 0];
o.n           = nOdor;
o.ec50        = nan;
o.neatMolarity= 1.67e-6; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = .83; %density (g/ml)
o.molecularWeight = 130.23; %molecular weight
o.molarFactor = 7.63e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'eugenol'};
o.colorCode   = [0 0 0.7];
o.plotSymbol      = ':v';
o.n           = nOdor;
o.ec50        = nan;
o.neatMolarity= 3.8e-7; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = 1.07; %density (g/ml)
o.molecularWeight = 164.2; %molecular weight
o.molarFactor = 1.23e-6;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'2-Hexanone'};
o.plotSymbol      = ':v';
o.colorCode   = [0 0.7 0];
o.n           = nOdor;
o.neatMolarity= 1.8e-4; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = .81; %density (g/ml)
o.molecularWeight = 100.16; %molecular weight
o.molarFactor = 632.01e-6;
o.ec50        = nan;
odor(nOdor)=o;

nOdor=nOdor+1;
o.name        = {'1-octanal'};
o.plotSymbol      = ':v';
o.colorCode   = [0.7 0 0];
o.n           = nOdor;
o.neatMolarity= nan; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = nan; %density (g/ml)
o.molecularWeight = nan; %molecular weight
o.molarFactor = 108.97e-6;
o.ec50        = nan;
odor(nOdor)=o;

% last odor is always for the unlisted ones and is orange.
nOdor=nOdor+1;
o.name       = {'unlisted'};
o.plotSymbol      = '--o';
o.colorCode  = [0.7 0.7 0.7];
o.neatMolarity= nan; %Molarity of saturated vapor in headspace of vial with pure odor
o.density     = nan; %density (g/ml)
o.molecularWeight = nan; %molecular weight
o.n          = nOdor;
o.ec50        = nan;
o.molarFactor=1;
odor(nOdor)=o;


%quick and dirty, save as a matlab structure in the current folder
if ~exist(fullfile(get_local_disk,'stimuli'))
    mkdir(get_local_disk,'stimuli')
end
save(fullfile(get_local_disk,'stimuli','odorDb.mat'),'odor');
fprintf('Odor database (%d odors) saved in %s\n',nOdor,fullfile(get_local_disk,'stimuli','odorDb.mat'));