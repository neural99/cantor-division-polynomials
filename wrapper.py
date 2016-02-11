import sage.interfaces.sage0

class CanHyperCurveWrapper:
    def __init__(self, instance_limit):
        self.instance_limit = instance_limit
        self.sage_instance = None
        self.create_new_instance()
        
    def compute_div_points(self, N, f, p):
        self.counter += 1
        if self.counter > self.instance_limit:
            self.create_new_instance()

        string_output = self.sage_instance.eval('H=CanHyperCurve(2, ' + str(f) + '); len(list(H.division_points(' + str(N) + ',' + str(p) + '))) - 1')
        #string_output = self.sage_instance.eval('H=CanHyperCurve(2, ' + str(f) + ');')
        return string_output
        
    def create_new_instance(self):
        if self.sage_instance:
            del self.sage_instance

        self.sage_instance = sage.interfaces.sage0.Sage()
        self.sage_instance('load("CanHyperCurve.sage")')
        self.counter = 0
