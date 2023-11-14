def mafft_paths(algorithm):
    return {'linsi':r'/davidb/local/bin/linsi'}.get(algorithm)

def tree_path(algorithm):
    return {'fast_tree':'FastTree'}.get(algorithm)