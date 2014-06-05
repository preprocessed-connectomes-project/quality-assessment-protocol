def gen_file_map(sink_dir, resource_paths=None):
    """
    Generate the mapping of different resources (e.g., functional_mni or sca) 
    to the file paths of those resources across subjects/scans. Note that this
    will only work with functional paths that have a scan_id in their path and
    will fail for any anatomical resources since those are assumed to only have
    one scan.
    
    Parameters
    ----------
    sink_dir : string
        Path to CPAC output or sink directory containing different pipeline 
        directories each with preprocessed data for individual subjects
    resource_paths : None or string or list (default = None)
        Path to file(s) in each subject's CPAC directory containing path information.
        Can include glob-style * or ? that will be expanded and can be a string or list.
        If this is None, then it will autoset to `os.path.join(sink_dir, 'pipeline_*', '*', 'path_files_here', '*.txt')`.
    
    Returns
    -------
    subject_infos : dict
        Mapping of output types (e.g., functional_mni) to information on each subject's associated output.
        Information on each subject's output is a list of length 4 and this includes pipeline_id, subject_id, 
        scan_id, and subject_path.
    """
    import os
    from glob import glob
    from CPAC.pipeline.cpac_group_runner import split_folders
    
    if resource_paths is None:
        resource_paths = os.path.join(sink_dir, 'pipeline_*', '*', 'path_files_here', '*.txt')
    
    if isinstance(resource_paths, str):
        resource_paths = [resource_paths]
    
    # Look through resource paths and expand
    # Then read in file paths
    subject_paths = []
    for resource_path in resource_paths:
        files = glob(os.path.abspath(resource_path))
        for file in files:
            path_list = open(file, 'r').readlines()
            subject_paths.extend([s.rstrip('\r\n') for s in path_list])
    
    if len(subject_paths) == 0:
        raise Exception("No subject paths found based on given resource_paths")
    
    # Remove any duplicate paths
    set_subject_paths   = set(subject_paths)
    subject_paths       = list(set_subject_paths)
    
    # Setup an mapping between resource and path info
    # so a dictionary with a list as the default value
    from collections import defaultdict
    analysis_map    = defaultdict(list)
    
    # Parse each subject path into 4 relevant parts: 
    # pipeline_id, subject_id, scan_id, and subject_path
    for subject_path in subject_paths:
        # Removes the base path
        if subject_path.find(sink_dir) == -1:
            print "WARNING: Couldn't find sink_dir: %s in subject's path: %s" % (sink_dir, subject_path)
        rs_path     = subject_path.replace(sink_dir, "", 1)
        rs_path     = rs_path.lstrip('/')
        
        # Split the path into a list (of folders)
        folders     = split_folders(rs_path)
                
        # If there aren't at least 4 sub-folders, then something is amiss
        if len(folders) < 3:
            raise Exception("Incorrect subject path, need 3-4 but only %i sub-folders found: %s" % (len(folders), subject_path))
        
        # Extract the desired elements
        pipeline_id = folders[0]
        subject_id  = folders[1]
        resource_id = folders[2]    # e.g., functional_mni, falff, etc
        if len(folders) == 3:
            scan_id = ""
        else:
            scan_id = folders[3]
        
        # Add to the mappings (seperate one for group analysis)
        # Note that the key is actually a tuple of the resource_id and key
        key         = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))
    
    return analysis_map

