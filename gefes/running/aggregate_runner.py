# Built-in modules #
import platform

# Internal modules #
import gefes

# First party modules #
import plumbing
from plumbing.runner import Runner

# Third party modules #
hostname = platform.node()

###############################################################################
class AggregateRunner(Runner):
    """Will run stuff on an aggregate"""
    modules = [gefes, plumbing]
    default_time = '7-00:00:00'

    default_steps = [
        {'assembly.run':     {}},
    ]

    @property
    def job_name(self): return "agg_%s" % self.parent.name

    def extra_slurm_params(self):
        # Standard cases #
        params = {'partition':  'node',
                  'constraint': 'mem512GB'}
        # Special cases #
        if self.parent.project.name == 'test':
            params['time'] = '01:00:00'
            if hostname.startswith('milou'): params['partition'] = 'devel'
        # Return result #
        return params

    def command(self, steps):
        from gefes.groups.projects import Project
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.parent.name]
        if isinstance(self.parent, Project): command += ["aggregate = gefes.projects[name]"]
        else: command += ["aggregate = getattr(gefes.groups.favorites, name)"]
        command += ["aggregate.load()"]
        command += ["aggregate.runner.run(steps)"]
        return command