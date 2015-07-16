%quick and dirty
%from cellsArray get a unit meta and run the visualize response
cp = cell_passport_tools();
%open cellsArray
fishy_cells = {'KPawakeM72_014_b_012','KPawakeM72_014_c_012', 'ZKawakeM72_004_h_002', 'ZKawakeM72_004_i_002', ...
    'ZKawakeM72_013_e_011', 'ZKawakeM72_012_a_001'}

fishy_cells = {'KPawakeM72_019_a_001', 'KPawakeM72_019_b_001', 'KPawakeM72_024_b_002', 'KPawakeM72_024_c_002', ...
    'KPawakeM72_817_e_001', 'KPawakeM72_817_f_001', 'ZKawakeM72_022_d_001'}

%fishy_cells = {'ZKawakeM72_004_h_002'}
%fishy_select = cellfun(@(x) find(strcmpi(x, {cellsArray.Id})), fishy_cells)
%select cellsarray
cellsForPass = cellsArray(([cellsArray.light]==1))
%cellsForPass = cellsArray(fishy_select)
for i=1:numel(cellsForPass)
    a_unit=cellsForPass(i);
    try
    pa = cp.make_passport(a_unit);
    catch me
        warning('Could not finish passport for cell %s',a_unit.Id)
        continue
    end
end


%Load all the cell files
mice = {'ZKawakeM72','KPawakeM72'};
ffn = file_names();
%mice = {'KPawakeM72'};
%load all the cells into an array of cells
disp(wrap_message('Getting cells and trial structures for Neil','*'));
cellFilesArray = [];
for is=1:numel(mice)
    cs=mice{is};
    cellFiles=dir(fullfile(ffn.fold_exp_data, [cs '*cell.mat']));
    allCells = arrayfun(@(x) load(fullfile(ffn.fold_exp_data,x.name)),cellFiles,'UniformOutput',false);
    cellFilesArray = [cellFilesArray [allCells{:}]]; %#ok<AGROW>
end

