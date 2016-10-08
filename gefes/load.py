b"""This module needs Python 2.7.x"""

# Built-in modules #
import os, glob

# Internal modules #
from gefes import samples, projects, _projects
from gefes.groups.projects import Project
from gefes.groups.samples  import Sample

###############################################################################
def load(json_dir_path, raw_files_must_exist=True):
    """Will load all the JSON files found in the given json_dir_path.
    If all the files are from the same project, return that project."""

    # Expand the tilda #
    if "~" in json_dir_path: json_dir_path = os.path.expanduser(json_dir_path)
    if not json_dir_path.endswith('/'): json_dir_path += '/'

    # Find all JSON files #
    json_paths  = glob.glob(json_dir_path + '*.json')
    if not json_paths: raise Exception("Did not find any json files at '%s'" % json_dir_path)

    # Load all JSON files #
    new_samples = [Sample(j, raw_files_must_exist) for j in json_paths]

    # Check names are unique #
    names = [s.short_name for s in new_samples]
    assert len(names) == len(set(names))

    # Add them to the master list #
    samples.extend(new_samples)

    # Compose into projects #
    proj_names = sorted(list(set([s.project_short_name for s in new_samples])))
    new_projs  = [Project(name, [s for s in new_samples if s.project_short_name==name]) for name in proj_names]
    _projects.extend(new_projs)

    # Link the samples to their project #
    for s in new_samples: s.project = projects[s.project_short_name]

    # Return project #
    proj_name = set(s.project_short_name for s in new_samples)
    if len(proj_name) == 1: return projects[proj_name.pop()]