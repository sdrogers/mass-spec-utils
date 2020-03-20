from molmass import Formula
import re

PROTON_MASS = 1.00727645199076
ELECTRON_MASS = Formula('H').isotope.mass - PROTON_MASS
adduct_rules = {
    '[M+H]+': {
        'mass_transform': Formula('H').isotope.mass - ELECTRON_MASS,
        'charge': 1
    },
    '[M+2H]2+': {
        'mass_transform': 2*(Formula('H').isotope.mass - ELECTRON_MASS),
        'charge': 2
    },
    '[M+2H+Na]3+' : {
        'mass_transform': 2*(Formula('H').isotope.mass - ELECTRON_MASS) + 
                          Formula('Na').isotope.mass - ELECTRON_MASS,
        'charge': 3
    },
    '[M+H+Na2]3+' : { 
        'mass_transform': (Formula('H').isotope.mass - ELECTRON_MASS) + 
                          2*(Formula('Na').isotope.mass - ELECTRON_MASS),
        'charge': 3
    },
    '[M+NH4]+' : {
        'mass_transform': (Formula('NH4').isotope.mass - ELECTRON_MASS),
        'charge': 1
    },
    '[M+H+NH4]2+': {
        'mass_transform': (Formula('H').isotope.mass - ELECTRON_MASS) + 
                          (Formula('NH4').isotope.mass - ELECTRON_MASS),
        'charge': 2
    },

    # -ve from this point
    '[M-H]-': {
        'mass_transform': -(Formula('H').isotope.mass - ELECTRON_MASS),
        'charge': -1
    },
    '[M-2H]2-': {
        'mass_transform': -2*(Formula('H').isotope.mass - ELECTRON_MASS),
        'charge': -2
    }

}


def mass2ion(mass,adduct_name):
    return (mass + adduct_rules[adduct_name]['mass_transform'])/abs(adduct_rules[adduct_name]['charge'])

def ion2mass(ion_mass,adduct_name):
    return ion_mass * abs(adduct_rules[adduct_name]['charge']) - adduct_rules[adduct_name]['mass_transform']

def get_transform_list():
    return list(adduct_rules.keys())

def get_positive_transform_list():
    return [a for a,v in adduct_rules.items() if v['charge'] > 0]

def get_negative_transform_list():
    return [a for a,v in adduct_rules.items() if v['charge'] < 0]



def adduct_string_parser(adduct_string):
    # Step 1, access the charge
    
    tokens = adduct_string.split(']')
    adduct_info = tokens[0]
    charge = tokens[1]

    if charge == '+':
        charge = 1
    elif charge == '-':
        charge = -1
    elif '+' in charge:
        charge = charge.replace('+','')
        charge = int(charge)
    elif '-' in charge:
        charge = charge.replace('-','')
        charge = -int(charge)
    else:
        charge = int(charge)

    adduct_info = adduct_info[2:] # strip the '[M'

    plus_pos = [('+',m.start()) for m in re.finditer("\+",adduct_info)]
    minus_pos = [('-',m.start()) for m in re.finditer("\-",adduct_info)]
    
    delim_positions = sorted(plus_pos + minus_pos,key = lambda x: x[1])
    
    mass_shift = 0
    for i,(t,d) in enumerate(delim_positions):
        if i == len(delim_positions) - 1:
            end_pos = len(adduct_info)
        else:
            end_pos = delim_positions[i+1][1]
        adduct = adduct_info[d+1:end_pos]

        # check for a pre-multiplier 
        # e.g. [M+2H]2+
        # we need to add the H twice.
        # Using Formula('2H') already removes an electron
        r = re.search('^[0-9]+',adduct)
        if not r is None:
            # get the multiplier
            multiplier = int(r.group())
            # extract the remaining adduct
            adduct = adduct[r.span()[1]:]
        else:
            multiplier = 1

        # get the mass
        formula = Formula(adduct)        
        adduct_mass = formula.isotope.mass

        # subtract if necessary
        if t == '-':
            adduct_mass = -adduct_mass

        # cumulative mass shift
        mass_shift += multiplier * adduct_mass
    
    # remove / add enough electrons
    mass_shift -= charge*ELECTRON_MASS

    return (mass_shift,charge)    


if __name__ == '__main__':
    print("E mass: ",ELECTRON_MASS)

    a_list = get_positive_transform_list()
    print("Positive transformations:")
    for a in a_list:
        print(a)

    a_list = get_negative_transform_list()
    print("Negative transformations:")
    for a in a_list:
        print(a)

    for adduct in adduct_rules:
        print()
        print("Transform 100 with {}: {}".format(adduct,mass2ion(100,adduct)))
        print("Transform back: ",ion2mass(mass2ion(100,adduct),adduct))


    for adduct in adduct_rules:
        print()
        print(adduct)
        print(adduct_string_parser(adduct))
        print(adduct_rules[adduct])

   
