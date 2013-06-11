function CreateString = SQL_TableCreateString(whattable)
% Determines the appropriate mySQL CREATE TABLE statement to use to create a given table, identified by the input string, whattable
% Ben Fulcher, May 2013

switch whattable
case 'Operations'
    CreateString = ['CREATE TABLE Operations ' ...
        '(m_id INTEGER NOT NULL AUTO_INCREMENT, ' ... % Unique integer identifier
        'OpName VARCHAR(255), ' ... % Unique name for the operation
        'MasterLabel VARCHAR(255), ' ... % Label of master code
        'Pointer TINYINT(1), ' ... % perhaps redundant given MasterLabel
        'Code VARCHAR(255), ' ... % Code to execute, or Master to retrieve from
        'Keywords VARCHAR(255), ' ... % Comma separated keyword metadata
        'CanDistribute INTEGER UNSIGNED, ' ... % Code for whether safe to distribute
        'LicenseType INTEGER UNSIGNED, ' ... % License for the code
        'Stochastic TINYINT(1), ' ... % Boolean identifier: is it a stochastic algorithm?
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, ' ... % Last modified
        'PRIMARY KEY (m_id))']; % sets primary key as m_id

case 'TimeSeries'
    CreateString = ['CREATE TABLE TimeSeries ' ...
        '(ts_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY, ' ... % Unique integer identifier
        'Filename VARCHAR(255) NOT NULL, ' ... % FileName of the time series
        'Keywords VARCHAR(255), ' ... % Comma-delimited keywords assigned to the time series
        'Length INTEGER UNSIGNED, ' ... % Length of the time series
        'Quantity VARCHAR(255), ' ... % The quantity measured in the time series
        'Unit VARCHAR(25), ' ... % The physical unit of the quantity measured
        'SamplingRate VARCHAR(25), ' ... % Sampling rate of the time series
        'Description TEXT, ' ... % More information about this specific time series
        'SourceString VARCHAR(255), ' ... % Name of source
        'Source_id INTEGER, '  ... % Source ID
        'CategoryString VARCHAR(255), ' ... % Name of category
        'Category_id INTEGER, ' ... % Category ID
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, ' ... % Time stamp of when the time series was last modified/inserted
        'FOREIGN KEY (Source_id) REFERENCES TimeSeriesSource(Source_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (Category_id) REFERENCES TimeSeriesCategories(Category_id) ON DELETE CASCADE ON UPDATE CASCADE)'];

case 'MasterOperations'
    CreateString = ['CREATE TABLE MasterOperations ' ...
        '(Master_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY, ' ... % Unique integer identifier
        'MasterLabel VARCHAR(255), ' ... % Name given to master code file
        'MasterCode VARCHAR(255), ' ... % Code to execute
        'NPointTo INTEGER UNSIGNED, ' ... % Number of children
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP)']; % Time stamp of when entry was last modified
        
case 'MasterPointerRelate'
    CreateString = ['CREATE TABLE MasterPointerRelate ' ...
        '(Master_id INTEGER, ' ... % Unique integer identifier
        'm_id INTEGER, ' ...
        'FOREIGN KEY (m_id) REFERENCES Operations(m_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (Master_id) REFERENCES MasterOperations(Master_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
        
case 'OperationKeywords'
    CreateString = ['CREATE TABLE OperationKeywords ' ...
        '(mkw_id INTEGER NOT NULL AUTO_INCREMENT PRIMARY KEY, ' ...
        'Keyword VARCHAR(255), ' ...
        'NumOccur INTEGER, ' ...
        'PercentageCalculated FLOAT, ' ...
        'PercentageGood FLOAT, ' ...
        'MeanCalcTime FLOAT)'];
        
case 'OpKeywordsRelate'
    CreateString = ['CREATE TABLE OpKeywordsRelate ' ...
        '(mkw_id INTEGER,' ...
        'm_id INTEGER, '  ...
        'FOREIGN KEY (mkw_id) REFERENCES OperationKeywords (mkw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (m_id) REFERENCES Operations (m_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
        
case 'TimeSeriesKeywords'
    CreateString = ['CREATE TABLE TimeSeriesKeywords ' ...
        '(tskw_id INTEGER AUTO_INCREMENT PRIMARY KEY, ' ... % Unique identifier for each keyword
        'Keyword varchar(50), ' ... % The keyword
        'NumOccur INTEGER, ' ... % Number of time series with this keyword
        'MeanLength INTEGER, ' ... % Mean length of time series with this keyword
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP)'];

case 'TimeSeriesCategories' % Hierarchical categories to assign different sources of time-series data
    CreateString = ['CREATE TABLE TimeSeriesCategories ' ...
        '(Category_id INTEGER AUTO_INCREMENT PRIMARY KEY, ' ...
        'CategoryName VARCHAR(255) UNIQUE, ' ... % Name of the category
        'Description TEXT, ' ... % Description of the category
        'ParentString VARCHAR(255), ' ... % Temporary string to label the CategoryName
        'Parent_id INTEGER, ' ... % Later assign proper id tags for faster indexing, joins, etc.
        'NMembers INTEGER UNSIGNED, ' ... % Automatically fill the parent_id from the string CategoryName
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP)'];

case 'TimeSeriesDistributionCodes'
    CreateString = ['CREATE TABLE TimeSeriesDistributionCodes ' ...
        '(Dcode_id INTEGER PRIMARY KEY, ' ... % A unique code for each distribution type
        'Description VARCHAR(255))']; % Short description of each distribution type
            
case 'TimeSeriesSource' % Table defining all the sources of time series in the database
    CreateString = ['CREATE TABLE TimeSeriesSource ' ...
        '(Source_id INTEGER AUTO_INCREMENT PRIMARY KEY, ' ...
        'SourceName VARCHAR(255), ' ... % Where the time series was sourced from
        'Link TEXT, ' ... % URL (or other) reference to time series source
        'Description TEXT, ' ... % Description of the source
        'ProcessingNote TEXT, ' ... % Notes on how time series from this source were processed
        'AddedBy VARCHAR(255), ' ... % Who added this source to the database ...
        'Dcode_id INTEGER, ' ... % Code for whether the data is safe to distribute
        'Parent VARCHAR(255), ' ... % Parent source
        'SourceParent_id INTEGER, ' ... % Parent source id
        'NMembers INTEGER UNSIGNED, ' ... % Number of time series from this source
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, ' ... % Add a time stamp of modifications
        'FOREIGN KEY (Dcode_id) REFERENCES TimeSeriesDistributionCodes(Dcode_id) ON DELETE CASCADE ON UPDATE CASCADE)'];

case 'TimeSeriesSourceRelate' % Links the time series with their sources
    CreateString = ['CREATE TABLE TimeSeriesSourceRelate ' ...
        '(ts_id INTEGER, ' ...
        'Source_id INTEGER, ' ...
        'FOREIGN KEY (ts_id) REFERENCES TimeSeries(ts_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (Source_id) REFERENCES TimeSeriesSource(Source_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
    
case 'TsKeywordsRelate'
    CreateString = ['CREATE TABLE TsKeywordsRelate ' ...
        '(tskw_id INTEGER, ' ...
        'ts_id INTEGER, ' ...
        'FOREIGN KEY (tskw_id) REFERENCES TimeSeriesKeywords(tskw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (ts_id) REFERENCES TimeSeries(ts_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
    
case 'Results'
    CreateString = ['CREATE TABLE Results ' ...
        '(ts_id integer, ' ...
        'm_id INTEGER, ' ...
        'Output DOUBLE, ' ...
        'QualityCode INTEGER UNSIGNED, ' ...
        'CalculationTime FLOAT, ' ...
        'LastModified TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP, ' ...
        'FOREIGN KEY (ts_id) REFERENCES TimeSeries(ts_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
        'FOREIGN KEY (m_id) REFERENCES Operations(m_id) ON DELETE CASCADE ON UPDATE CASCADE, '...
        'PRIMARY KEY(ts_id,m_id))'];


otherwise
    error(['Unknown table ' whattable]) 
end

end