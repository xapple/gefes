# Built-in modules #

# Internal modules #
import gefes
from plumbing.runner import Runner

# Third party modules #

###############################################################################
class AggregateRunner(Runner):
    """Will run stuff on an aggregate"""
    default_time = '2-00:00:00'

    default_steps = [
        {'assembly.run':     {}},
    ]

    @property
    def job_name(self): return "agg_%s" % self.parent.name

    def command(self, steps):
        from gefes.groups.projects import Project
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.parent.name]
        if isinstance(self.parent, Project): command += ["aggregate = gefes.projects[name]"]
        else: command += ["aggregate = getattr(gefes.groups.favorites, name)"]
        command += ["aggregate.load()"]
        command += ["aggregate.runner.run(steps)"]
        return command

    def run_slurm(self, steps=None, **kwargs):
        return Runner.run_slurm(self, steps=steps, module=gefes, **kwargs)