# Material stream class
# need to add cooprop to get thermo_physical properties

import CoolProp as CP

from bisect import bisect_left

class material_stream:
    '''
    A material stream class that stores state variables of a material stream.
    The class also returns the thermodynamic and physical properties of the stream.
    The streams using this class can be mixed and split.

    Inputs: composition; state
    composition: a dictionary object containing name of component as key and its mass frac as value.
    e.g. composition = {'Water': 0.10, 'Ethanol': 0.90}
    state: a dictionary contaning temperature (temp, K), pressure (pres, Pa), mole flowrate (flow, mole/s).
    e.g. state = {'temp': 298.15, 'pres': 101.325, 'mole': 5}

    '''

    def __init__(self, composition, state):
        self.composition=composition
        self.state=state
        comp_name='&'.join(list(self.composition.keys()))
        comp_values=list(self.composition.values())

        # defining thermodynamic state for CoolProp
        system = CP.AbstractState("HEOS",comp_name)
        system.set_mole_fractions(comp_values)
        system.update(CP.PT_INPUTS, self.state['pres'], self.state['temp'])
        #self.tp_para = [system.keyed_output(k) for k in [CP.iQ, CP.iDmass , CP.iDmolar, CP.iHmass, CP.iHmolar, CP.iSmolar, CP.iCpmolar, CP.iviscosity, CP.iconductivity]]
        
        # Getting property parameter from CoolProp
        self.Phase=system.keyed_output(CP.iPhase)   # vapor fraction
        self.Q=system.keyed_output(CP.iQ)           # vapor fraction
        self.Dmass=system.keyed_output(CP.iDmass)   # mass density
        self.Dmolar=system.keyed_output(CP.iDmolar) # molar density
        self.Hmass=system.keyed_output(CP.iHmass)   # specific enthalpy mass
        self.Hmolar=system.keyed_output(CP.iHmolar)     # specific enthalpy molar
        self.Cpmass=system.keyed_output(CP.iCpmass)     # specific heat mass
        self.Cpmolar=system.keyed_output(CP.iCpmolar)   # specific heat molar
        try:
            self.Vis=system.keyed_output(CP.iviscosity)     # viscosity
        except: 
            self.Vis=[]
        try:
            self.cond=system.keyed_output(CP.iconductivity)     # thermla conductivity
        except:
            self.cond=[]     # thermla conductivity
        

    def no_component(self):
        return len(self.composition)
    
    def __repr__(self):
        return 'Stream containing {!r}'.format(list(self.composition.keys()))
    
    def takeClosest(self,myList, myNumber):
        """
        Assumes myList is sorted. Returns closest value to myNumber.

        If two numbers are equally close, return the smallest number.
        """
        self.myList=myList
        self.myNumber=myNumber

        pos = bisect_left(myList, myNumber)
        if pos == 0:
            return myList[0]
        if pos == len(myList):
            return myList[-1]
        before = myList[pos - 1]
        after = myList[pos]
        if after - myNumber < myNumber - before:
            return after
        else:
            return before

      
    def __add__(self, other):
        self_flow=self.state['flow']
        other_flow=other.state['flow']
        new_flow=self_flow+other_flow
        
        self_comp=self.composition
        other_comp=other.composition
        try:
            self_comp = {key: comp * self_flow for key, comp in self_comp.items()}
            other_comp = {key: comp * other_flow for key, comp in other_comp.items()}
            new = self_comp.copy()
            new.update(other_comp)
            
            for key in new.keys():
                first=self_comp.get(key,0)
                second=other_comp.get(key,0)
                new[key]=(first+second)/(new_flow)
            
            # State calculation
            new_state = self.state.copy()
            new_state['flow'] = new_flow # total flow
            # pressure = lowest pressure among streams
            if self.state['pres'] <= other.state['pres']:
                new_state['pres'] = self.state['pres']
            else:
                new_state['pres'] = other.state['pres']

            # temperature from energy balance--isenthalpic mixing
            # Note: heat of mixing is not considered here
            
            # Specific molar enthalpy of mixture
            new_Hmolar = (self.Hmolar*self.state['flow'] + other.Hmolar*other.state['flow'])/new_flow
            
            # State of mixture
            new_comp_name='&'.join(list(new.keys()))
            new_comp_values=list(new.values())
            new_system = CP.AbstractState("HEOS",new_comp_name)
            new_system.set_mole_fractions(new_comp_values)
            # Make phase envelop
            new_system.build_phase_envelope("dummy")
            PE = new_system.get_phase_envelope_data()
            val=self.takeClosest(PE.hmolar_liq, new_Hmolar)
            ind=PE.hmolar_liq.index(val)
            new_state['temp']=PE.T[ind]
            # State with enthalpy and pressure as inputs
            #new_system.update(CP.HmolarP_INPUTS , new_Hmolar, new_state['pres'])
            #New temp
            #new_state['temp'] = new_system.keyed_output(CP.iT)
        except:
            ValueError

        return material_stream(new, new_state)
        # return Polynomial(*(x + y for x, y in zip(self.coeffs, other.coeffs)))


if __name__ == '__main__':
    # stream 1
    comp1 = {'Water': 0.70, 'Acetone': 0.30}
    state1 = {'temp': 391, 'pres': 4e5, 'flow': 2}
    
    # stream 2
    comp2 = {'Water': 0.70, 'Ethanol': 0.30}
    state2 = {'temp': 300, 'pres': 1.01e5, 'flow': 15}
    # stream 3
    comp3 = {'Water': 1.0}
    state3 = {'temp': 320, 'pres': 101325, 'flow': 2}

    # instantiate the material_stream class
    MS100 = material_stream(composition=comp1, state=state1)
    #print(help(MS100))
    #print('Number of component: ',MS100.no_component())
    print(repr(MS100))
    print('Density: {0:.3f} kg/m3 {1:.3f} mol/m3'.format(MS100.Dmass, MS100.Dmolar))

    #print(MS100.composition)
        
    # add streams
    print('\nAdd two streams')
    MS200 = material_stream(composition=comp2, state=state2)
    MS_add2 = MS100+MS200
    print('New composition:\n', MS_add2.composition)
    print('New state: \n', MS_add2.state)
    
    # Adding three streams
    print('\nAdd three streams')
    MS300 = material_stream(composition=comp3, state=state3)
    MS_add3 = MS100 + MS200 + MS300
    print('New composition:\n', MS_add3.composition)
    print('New state: \n', MS_add3.state)