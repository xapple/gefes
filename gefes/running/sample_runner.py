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
    modules      = [gefes, plumbing]
    default_time = '2-00:00:00'

    default_steps = [
        {'pair.fwd.fastqc.run':                     {}},
        {'pair.rev.fastqc.run':                     {}},
        {'quality_checker.run':                     {}},
        {'clean.fwd.fastqc.run':                    {}},
        {'clean.rev.fastqc.run':                    {}},
        {'clean.fwd.graphs.length_dist.plot':       {'x_log': True, 'y_log': True}},
        {'clean.rev.graphs.length_dist.plot':       {'x_log': True, 'y_log': True}},
        {'clean.fwd.kraken.run':                    {}},
        {'clean.rev.kraken.run':                    {}},
        {'assembly.run':                            {}},
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
            if hostname.startswith('milou'):
                params['time'] = '00:15:00'
                params['qos'] = 'short'
            if hostname.startswith('taito'):
                params['time'] = '00:30:00'
                params['partition'] = 'test'
        # Return result #
        return params

    def command(self, steps):
        command =  ["steps   = %s"   % steps]
        command += ["p_name  = '%s'" % self.parent.project.name]
        command += ["s_name  = '%s'" % self.parent.name]
        command += ["project = gefes.projects[p_name].load()"]
        command += ["sample  = [s for s in project if s.name==s_name][0]".load()]
        command += ["sample.runner.run(steps)"]
        return command