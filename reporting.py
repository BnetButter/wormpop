from config import load_constants

constants = load_constants()

import collections

def CreateCounter():
    counter = 0
    mass_counter = 0
    def wrapper(func):
        def _wraps(self, *args, **kwargs):
            nonlocal counter
            nonlocal mass_counter
            val = func(self, *args, **kwargs)
            counter += 1
            mass_counter += self.mass
            return val
        return _wraps

    def reporter():
        nonlocal counter
        nonlocal mass_counter
        new_val = counter
        new_mass = mass_counter
        counter = 0
        mass_counter = 0
        return new_val, new_mass

    return wrapper, reporter

HatchSet, HatchGet = CreateCounter()

LarvaToDauerSet, LarvaToDauerGet = CreateCounter()
LarvaToAdultSet, LarvaToAdultGet = CreateCounter()

LarvaToL1ArrestSet, LarvaToL1ArrestGet = CreateCounter()
L1ArrestToLarvaSet, L1ArrestToLarvaGet = CreateCounter()

def get_column_default(column):
    if column.default is None:
        return None
    if isinstance(column.default, DefaultClause):
        if isinstance(column.default.arg, expression.Function):
            return None
        return column.default.arg
    return column.default.arg

def create_death_counter():
    """
    Create the dictionary needed to count mass and individuals that die
    """

    def create_cause_of_death():
        return {
            "arrested_development": 0,
            "starvation": 0,
            "old_age": 0,
            "culled": 0,
            "bag": 0,
            "end_of_simulation": 0,
        }

    return {
        "Egg": create_cause_of_death(),
        "Larva": create_cause_of_death(),
        "Adult": create_cause_of_death(),
        "Dauer": create_cause_of_death(),
        "Parlad": create_cause_of_death(),
        "L1Arrest": create_cause_of_death(),
    }

def CreateDeathCounter():
    counter = create_death_counter()
    mass_counter = create_death_counter()

    def die_wrapper(func):
        def die_fn(self, cause_of_death):
            counter[self.__class__.__name__][cause_of_death] += 1
            mass_counter[self.__class__.__name__][cause_of_death] += self.mass
            return func(self, cause_of_death)
        return die_fn
    
    def reporter():
        nonlocal counter
        nonlocal mass_counter
        tmp_counter = counter
        tmp_mass_counter = mass_counter
        counter = create_death_counter()
        mass_counter = create_death_counter()
        return tmp_counter, tmp_mass_counter
    return die_wrapper, reporter

die_wrapper, die_reporter = CreateDeathCounter()