# Built-in modules #
import platform

# Internal modules #
import gefes

# First party modules #
import plumbing
from plumbing.runner import Runner

# Third party modules #

# Constants #
hostname = platform.node()

###############################################################################
class SampleRunner(Runner):
    """Will run stuff on a sample"""
    modules = [gefes, plumbing]
    default_time = '2-00:00:00'

    default_steps = [
        {'pair.fwd.fastqc.run':                     {}},
        {'pair.rev.fastqc.run':                     {}},
        {'quality_checker.run':                     {}},
        {'clean.fwd.fastqc.run':                    {}},
        {'clean.rev.fastqc.run':                    {}},
        {'clean.fwd.graphs.length_dist.plot':       {'x_log': True, 'y_log': True}},
        {'clean.rev.graphs.length_dist.plot':       {'x_log': True, 'y_log': True}},
        {'report.generate':                         {}},
    ]

    @property
    def job_name(self): return "gefes_%s" % self.parent.name

    @property
    def extra_slurm_params(self):
        # Standard cases #
        params = {'partition': 'core',
                  'cores'    : 1}
        # Special cases #
        if self.parent.project.name == 'test':
            params['time'] = '00:15:00'
            if hostname.startswith('milou'): params['qos'] = 'short'
        # Return result #
        return params

    def command(self, steps):
        command =  ["steps = %s" % steps]
        command += ["s_name = '%s'" % self.parent.name]
        command += ["p_name = '%s'" % self.parent.project.name]
        command += ["project = [p for p in gefes.projects if p.name==p_name][0]"]
        command += ["project.load()"]
        command += ["sample  = [s for s in project        if s.name==s_name][0]"]
        command += ["sample.load()"]
        command += ["sample.runner.run(steps)"]
        return command