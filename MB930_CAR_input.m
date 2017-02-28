
% READ IN MB930 CAR data

MB930_infile_name_csv = 'HPM-20 input files\MB930_core_CAR_data.csv';
MB930_CAR_time_series = importdata(MB930_infile_name_csv);

MB930_depth = MB930_CAR_time_series.data(:,1);  % cm
MB930_age = MB930_CAR_time_series.data(:,2);  % yr
MB930_age_kyr = MB930_CAR_time_series.data(:,3); %kyr
MB930_CAR = MB930_CAR_time_series.data(:,4);  % g C/m2/y

% READ IN MB930 age-depth data

MB930_infile_name_csv = 'HPM-20 input files\MB930_core_age_depth_data.csv';
MB930_age_depth_time_series = importdata(MB930_infile_name_csv);
