from molmass import Formula

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
