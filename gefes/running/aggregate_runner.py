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
    default_time = '2-00:00:00'

    default_steps = [
        {'assembly.run':     {}},
    ]

    @property
    def job_name(self): return "agg_%s" % self.parent.name

    @property
    def extra_slurm_params(self):
        # Standard cases #
        params = {}
        # Special cases #
        if self.parent.project.name == 'test':
            params['time'] = '00:15:00'
            if hostname.startswith('milou'): params['qos'] = 'short'
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