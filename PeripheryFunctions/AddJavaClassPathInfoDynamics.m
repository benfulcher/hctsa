function AddJavaClassPathInfoDynamics()
    % Adds the path required for the infodynamics toolkit
    % Required for the parallel workers, which don't always inherit the proper
    % javaclasspath

    computeDir = which('TS_Compute');
    splits = regexp(computeDir,filesep);
    hctsaDir = computeDir(1:splits(end-1));
    javaaddpath(fullfile(hctsaDir,'Toolboxes','infodynamics-dist','infodynamics.jar'));
end
