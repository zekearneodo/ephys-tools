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

dims = DIM0;

%
% Initialize data.
%
wdata.serial_no   =int32([1153 ; 1184 ; 1027  ;    1313]);
wdata.temperature =[53.23; 55.12; 130.55; 1252.89];
wdata.pressure    =[24.57; 22.95;  31.23;   84.11];

%
%% Create a new file using the default properties.
%
file = H5F.create (fileName, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

%
% Create the compound datatype for memory.
%
memtype = H5T.create ('H5T_COMPOUND', H5ML.sizeof (wdata));

H5T.insert (memtype, 'serial_no',H5ML.hoffset(wdata(1), 'serial_no'), 'H5T_NATIVE_INT');
H5T.insert (memtype, 'temperature',H5ML.hoffset (wdata(1), 'temperature'),'H5T_NATIVE_DOUBLE');
H5T.insert (memtype, 'pressure',H5ML.hoffset (wdata(1), 'pressure'), 'H5T_NATIVE_DOUBLE');

%
% Create the compound datatype for the file.  Because the standard
% types we are using for the file may have different sizes than
% the corresponding native types, we must manually calculate the
% offset of each member.
%
filetype = H5T.create ('H5T_COMPOUND', 8 + 8 + 8);
H5T.insert (filetype, 'serial_no', 0, 'H5T_STD_I64BE');
H5T.insert (filetype, 'temperature', 8 ,'H5T_IEEE_F64BE');
H5T.insert (filetype, 'pressure', 8 + 8,'H5T_IEEE_F64BE');

%
% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
%
space = H5S.create_simple (1,fliplr( dims), []);

%
% Create the dataset and write the compound data to it.
%
dset = H5D.create (file, DATASET, filetype, space, 'H5P_DEFAULT');
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', wdata);

%
% Close and release resources.
%
H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5F.close (file);


%
%% Now we begin the read section of this example.  Here we assume
% the dataset has the same name and rank, but can have any size.
%

%
% Open file and dataset.
%
file = H5F.open (fileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
dset = H5D.open (file, DATASET);

%
% Get dataspace and allocate memory for read buffer.
%
space = H5D.get_space (dset);
[numdims dims maxdims] = H5S.get_simple_extent_dims (space);
dims = fliplr(dims');

%
% Read the data.
%
rdata=H5D.read (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');

%
% Output the data to the screen.
%
for i=1: dims(1)
    fprintf ('%s[%d]:\n', DATASET, i);
    fprintf ('Serial number   : %d\n', rdata.serial_no(i));
    fprintf ('Temperature (F) : %f\n', rdata.temperature(i));
    fprintf ('Pressure (inHg) : %f\n\n', rdata.pressure(i));
end

%
% Close and release resources.
%
H5D.close (dset);
H5S.close (space);
H5T.close (memtype);
H5F.close (file);
