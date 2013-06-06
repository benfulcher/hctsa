# IMPORTS ALL TIME-SERIES DATA INTO MYSQL TABLES

# Import time series distribution codes from INP_dcodes.csv into the TimeSeriesDistributionCodes table
LOAD DATA LOCAL INFILE 'INP_dcodes.csv' INTO TABLE TimeSeriesDistributionCodes
FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' ESCAPED BY '"' LINES TERMINATED BY '\r'
IGNORE 1 LINES
(Dcode_id, Description);

# Import source information from INP_sources.csv into the TimeSeriesSource table
LOAD DATA LOCAL INFILE 'INP_sources.csv' INTO TABLE TimeSeriesSource
FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' ESCAPED BY '"' LINES TERMINATED BY '\r'
IGNORE 1 LINES
(SourceName, Link, Description, ProcessingNote, AddedBy, Dcode_id, Parent);

# Import all category information from INP_categories.csv into the TimeSeriesCategories table
LOAD DATA LOCAL INFILE 'INP_categories.csv' INTO TABLE TimeSeriesCategories
FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' ESCAPED BY '"' LINES TERMINATED BY '\r'
IGNORE 1 LINES
(CategoryName, Description, ParentString);

# Import all time-series information from INP_ts.csv into the TimeSeries table
LOAD DATA LOCAL INFILE 'INP_ts.csv' INTO TABLE TimeSeries
FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' ESCAPED BY '"' LINES TERMINATED BY '\r'
IGNORE 1 LINES
(FileName, Keywords, Quantity, Unit, SamplingRate, Description, SourceString, CategoryString);
