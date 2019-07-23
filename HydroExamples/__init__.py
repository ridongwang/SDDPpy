
'''
Objects for Hydro scheduling examples
'''
def import_SDDP():
    import os 
    import sys
    hydro_path = os.path.dirname(os.path.realpath(__file__))
    parent_path= os.path.abspath(os.path.join(hydro_path, os.pardir))
    print(parent_path)
    sys.path.append(parent_path)
class Turbine():
    def __init__(self, flowknots, powerknots):
        self.flowknots = flowknots
        self.powerknots = powerknots
class Reservoir():
    def __init__(self, minlevel, maxlevel, initial, turbine, s_cost, inflows):
        self.min = minlevel
        self.max = maxlevel
        self.initial = initial
        self.turbine = turbine
        self.spill_cost = s_cost
        self.inflows = inflows 

dro_radii = [b*(10**c) for c in [-1,0,1,2] for b in [1,2,3,4,5,6,7,8,9]]#[1,4,7]]