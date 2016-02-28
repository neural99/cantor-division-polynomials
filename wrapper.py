import sage.interfaces.sage0
import sage.interfaces.quit

import os

class CanHyperCurveWrapper:
    def __init__(self, instance_limit, algo):
        self.instance_limit = instance_limit
        self.algo = algo
        self.sage_instance = None
        self.create_new_instance()

    def reap_zombies(self):
        while True:
            try:
                pid, exit_code = os.waitpid(-1, 0)
                #print "Reaping " + str(pid)
            except OSError:
                break
        
    def compute_div_points(self, N, f, p):
        self.counter += 1
        if self.counter > self.instance_limit:
            self.create_new_instance()

        if self.algo == 'cantor':
            string_output = self.sage_instance.eval('H=CanHyperCurve(2, ' + str(f) + '); len(list(H.division_points(' + str(N) + ',' + str(p) + '))) - 1')
        else:
            string_output = self.sage_instance.eval('H=CanHyperCurve(2, ' + str(f) + '); len(list(H.division_points_naive(' + str(N) + ',' + str(p) + '))) - 1')
            
        return string_output
        
    def create_new_instance(self):
        sage.interfaces.quit.expect_quitall(verbose=False)
        self.reap_zombies()

        if self.sage_instance:
            del self.sage_instance

        self.sage_instance = sage.interfaces.sage0.Sage()
        self.sage_instance('load("CanHyperCurve.sage")')
        self.counter = 0
