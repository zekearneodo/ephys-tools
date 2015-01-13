%**************************************************************************
%
%  This example shows how to read and write compound
%  datatypes to a dataset.  The program first writes
%  compound structures to a dataset with a dataspace of DIM0,
%  then closes the file.  Next, it reopens the file, reads
%  back the data, and outputs it to the screen.
%
%  This file is intended for use with HDF5 Library verion 1.6
%
%**************************************************************************

fileName       = 'h5ex_t_cmpd.h5';
DATASET        = 'DS1';
DIM0           = 4;
SDIM           = 32;
dims = DIM0;

%
% Initialize data.
%
wdata.serial_no   =int32([1153 ; 1184 ; 1027  ;    1313]);
wdata.temperature =[53.23; 55.12; 130.55; 1252.89];
wdata.pressure    =[24.57; 22.95;  31.23;   84.11];
wdata.texts      = ['Parting'; 'is such'; 'sweet  '; 'sorrow.'];

 
%% Create a new file using the default properties.
%
file = H5F.create (fileName, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

%
% Create the compound datatype for memory.
%
dtypenames=fields(wdata);
mtypetag={'H5T_NATIVE_INT','H5T_NATIVE_DOUBLE','H5T_NATIVE_DOUBLE','H5T_C_S1'};
ftypetag={'H5T_STD_I64BE','H5T_IEEE_F64BE','H5T_IEEE_F64BE','H5T_FORTRAN_S1'}

%make the compound for the memory (with native types) 
for i=1:numel(mtypetag)
    %native types for the memory data type
    %standard types for the file data type
    %types
    mtype(i)  = H5T.copy(mtypetag{i});
    ftype(i) = H5T.copy(ftypetag{i});
    %sizes
    if ischar(wdata.(dtypenames{i}))
        dim = SDIM;
    else
        dim = 1;
    end
    msz(i)    = dim*H5T.get_size(mtype(i));
    fsz(i)    = dim*H5T.get_size(ftype(i));
end
moffset = [0 cumsum(msz(1:end-1))];
foffset = [0 cumsum(fsz(1:end-1))];

%create the compound for memory and file
struct_mtype = H5T.create('H5T_COMPOUND',sum(msz));
struct_ftype = H5T.create('H5T_COMPOUND',sum(fsz));

for i=1:numel(mtypetag)
    H5T.insert(struct_mtype,dtypenames{i},moffset(i),mtype(i))
    H5T.insert(struct_ftype,dtypenames{i},foffset(i),ftype(i))
end



%
% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
%
space = H5S.create_simple (2,( size(wdata.(dtypenames{1}))), []);

%
% Create the dataset and write the compound data to it.
%
dset = H5D.create (file, DATASET, struct_ftype, space, 'H5P_DEFAULT');
H5D.write (dset, struct_mtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', wdata);

%
% Close and release resources.
%
H5D.close (dset);
H5S.close (space);
H5T.close (struct_mtype);
H5F.close (file);

