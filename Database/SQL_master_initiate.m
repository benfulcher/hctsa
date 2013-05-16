function SQL_master_initiate
    % Create all the table
    
    createstring = {};
    name = {};
    
    createstring{end+1} = ['CREATE TABLE Operations (m_id integer not null auto_increment, MasterLabel varchar(255), OpName varchar(255), Pointer tinyint(1), ' ...
    'Code varchar(255), Keywords varchar(255), Stochastic tinyint(1), Normalize tinyint(2),  LastModified datetime, ' ...
    'PRIMARY KEY (m_id))'];
    name{end+1} = 'Operations';
    
    createstring{end+1} = ['CREATE TABLE MasterOperations ' ...
    '(Master_id integer not null auto_increment, MasterLabel varchar(255), ' ...
    'MasterCode varchar(255), NPointTo integer unsigned, LastModified datetime, ' ...
    'PRIMARY KEY (Master_id))'];
    name{end+1} = 'MasterOperations';

    createstring{end+1} = ['CREATE TABLE TimeSeries (ts_id integer not null auto_increment, Filename varchar(255), Keywords varchar(255), ' ...
    'Length integer unsigned, SamplingRate float, LastModified datetime, PRIMARY KEY (ts_id))'];
    name{end+1} = 'TimeSeries';

    createstring{end+1} = ['CREATE TABLE MasterPointerRelate (Master_id integer, m_id integer, ' ...
    'FOREIGN KEY (m_id) REFERENCES Operations(m_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
    'FOREIGN KEY (Master_id) REFERENCES MasterOperations(Master_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
    name{end+1} = 'MasterPointerRelate';

    createstring{end+1} = ['CREATE TABLE OperationKeywords ' ...
    '(mkw_id integer not null auto_increment primary key, Keyword varchar(255), ' ...
    'NumOccur integer, PercentageCalculated float, PercentageGood float, MeanCalcTime float)'];
    name{end+1} = 'OperationKeywords';

    createstring{end+1} = ['CREATE TABLE mkwFileRelate (mkw_id integer, m_id integer, '  ...
    'FOREIGN KEY (mkw_id) REFERENCES OperationKeywords (mkw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
    'FOREIGN KEY (m_id) REFERENCES Operations (m_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
    name{end+1} = 'mkwFileRelate';

    createstring{end+1} = ['CREATE TABLE TimeSeriesKeywords (tskw_id integer auto_increment primary key, Keyword varchar(50), ' ...
    'NumOccur integer, PercentageCalculated float, PercentageGood float, MeanCalcTime float, MeanLength integer)'];
    name{end+1} = 'TimeSeriesKeywords';

    createstring{end+1} = ['CREATE TABLE tskwFileRelate (tskw_id integer, ts_id integer, ' ...
    'FOREIGN KEY (tskw_id) REFERENCES TimeSeriesKeywords(tskw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
    'FOREIGN KEY (ts_id) REFERENCES TimeSeries(ts_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
    name{end+1} = 'tskwFileRelate';

    createstring{end+1} = ['CREATE TABLE Results (ts_id integer, m_id integer, Output double, Quality integer unsigned, CalculationTime float, LastModified datetime, ' ...
    'FOREIGN KEY (ts_id) REFERENCES TimeSeries (ts_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
    'FOREIGN KEY (m_id) REFERENCES Operations (m_id) ON DELETE CASCADE ON UPDATE CASCADE, '...
    'PRIMARY KEY(ts_id,m_id) )'];
    name{end+1} = 'Results';
    
    %createstring{end+1} = 'CREATE INDEX index_ts_id on Results(ts_id)';
    %name{end+1} = 'Results (index ts_id)';
    %createstring{end+1} = 'CREATE INDEX index_m_id on Results(m_id)';
    %name{end+1} = 'Results (index m_id)';
    
    %createstring{end+1} = 'ALTER TABLE Results ADD CONSTRAINT ts_id_m_id UNIQUE(ts_id,m_id)';
    %name{end+1} = 'Results (unique constraint on ts_id and m_id)';
    
    [dbc, dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)
    
    for j = 1:length(createstring)
        mysql_dbexecute(dbc,createstring{j});
        fprintf(1,'Created table: %s\n',name{j});
    end
    
    
