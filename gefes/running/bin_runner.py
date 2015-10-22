# Built-in modules #
import platform

# Internal modules #
import gefes

# First party modules #
from plumbing.runner import Runner

# Third party modules #

# Constants #
hostname = platform.node()

###############################################################################
class BinRunner(Runner):
    """Will run stuff on a bin."""
    modules      = [gefes]
    default_time = '1-00:00:00'

    default_steps = [
        {'taxonomy.run': {}},
    ]

    @property
    def job_name(self): return "bin_%s" % self.parent.name

    @property
    def extra_slurm_params(self): return {}

    def command(self, steps):
        command =  ["steps    = %s"   % steps]
        command += ["p_name   = '%s'" % self.parent.assembly.samples[0].project.name]
        command += ["a_name   = '%s'" % self.parent.assembly.name]
        command += ["b_name   = '%s'" % self.parent.name]
        command += ["project  = gefes.projects[p_name].load()"]
        command += ["assembly = [a for a in project.assemblies.values() if a.name==a_name][0]"]
        command += ["bin      = [b for b in assembly.results.binner.results.bins if b.name==s_name][0]"]
        command += ["bin.runner.run(steps)"]
        return command