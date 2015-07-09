%quick and dirty
%from cellsArray get a unit meta and run the visualize response
cp = cell_passport_tools();
%open cellsArray

%select cellsarray
cellsForPass = cellsArray(([cellsArray.light]==1))
for i=1:numel(cellsForPass)
    a_unit=cellsForPass(i);
    try
    pa = cp.make_passport(a_unit);
    catch me
        continue
    end
end
