# Built-in modules #

# Internal modules #
import gefes
from plumbing.runner import Runner

# Third party modules #

# Constants #

###############################################################################
class SampleRunner(Runner):
    """Will run stuff on an aggregate"""
    default_time = '02:00:00'

    default_steps = [
        {'pair.fwd.fastqc.run':                     {}},
        {'pair.rev.fastqc.run':                     {}},
        {'quality_checker.run':                     {}},
        {'clean.fwd.fastqc.run':                    {}},
        {'clean.rev.fastqc.run':                    {}},
        {'clean.fwd.graphs.length_dist.plot':       {}},
        {'clean.rev.graphs.length_dist.plot':       {}},
        {'report.generate':                         {}},
    ]

    @property
    def job_name(self): return "gefes_%s" % self.parent.name

    def command(self, steps):
        command =  ["steps = %s" % steps]
        command += ["name = '%s'" % self.parent.name]
        command += ["sample = [s for s in gefes.samples if s.name==name][0]"]
        command += ["sample.load()"]
        command += ["sample.runner.run(steps)"]
        return command

    def run_slurm(self, steps=None, **kwargs):
        return Runner.run_slurm(self, steps=steps, module=gefes, **kwargs)