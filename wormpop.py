#@Author: Matt Mosley
#2021-10-19
#%%

"""
USAGE: wormpop [--parameters=<string>] [ --database=<string> ] [ --name=<string> ] [--directory=<string>]
"""

import pathlib
import math
import numpy
from numpy import random
import csv
import json
import docopt
import sqlite3
import sys

args = docopt.docopt(__doc__)
parameters = args["--parameters"]
database = args["--database"] if args["--database"] is not None else ":memory:"
name = args["--name"] if args["--name"] is not None else "Simulation"
directory = args["--directory"] if args["--directory"] is not None else "Simulation"

# Constants:

import json

# Read the JSON file

if parameters:
    with open(parameters, 'r') as file:
        param = json.load(file)
else:
    param = json.load(sys.stdin)

# Simulation time details
SIMULATION_LENGTH = param['SIMULATION_LENGTH']  # 800 timesteps = 100 days
TIMESTEP = param['TIMESTEP']  # 3 hr per timestep, 8 timesteps per day

# Initial conditions
STARTING_WORMS = param['STARTING_WORMS']
STARTING_STAGE = param['STARTING_STAGE']  # 'egg'
EGGMASS = param['EGGMASS']  # Nanograms

# Adult constants
MIN_ADULT_MASS = param['MIN_ADULT_MASS']  # ng
MIN_ADULT_AGE = param['MIN_ADULT_AGE']  # Minimum age in hours to transition to adult
MAX_ADULT_AGE = param['MAX_ADULT_AGE']  # Max age in hours to transition (larvae past this age die of "arrested development")

# Larva constants
STANDARD_LARVA_MASS = param['STANDARD_LARVA_MASS']  # ~228 ng
LARVAL_STARVE_PROB = param['LARVAL_STARVE_PROB']  # Chance to cheat death by starvation, though starvation is probabilistic

# Dauer constants
MIN_DAUER_MASS = param['MIN_DAUER_MASS']  # ~137 ng
MAX_DAUER_MASS = param['MAX_DAUER_MASS']  # ~456 ng
DAUER_THRESHOLD = param['DAUER_THRESHOLD']  # Concentration (mg/mL) that scales probability of dauering = 250,000 ng total food available
DAUER_RATE = param['DAUER_RATE']  # Number of days at 0 food concentration for a larva to have a 50% chance of dauering/starving
DAUER_EXIT_PROB = param['DAUER_EXIT_PROB']  # Chance per timestep to exit dauer, based on empirical data

# Bag constants
BAG_THRESHOLD = param['BAG_THRESHOLD']  # mg/mL (=2500 ng)
BAG_RATE = param['BAG_RATE']
BAG_EFFICIENCY = param['BAG_EFFICIENCY']  # Efficiency with which somatic mass of parlads can be converted to dauers

# Food constants
STARTING_FOOD = param['STARTING_FOOD']  # 10 mg = 1x10^7 ng
FEEDING_AMOUNT = param['FEEDING_AMOUNT']  # 10 mg added per feeding schedule

# Scheduling constants
FEEDING_SCHEDULE = param['FEEDING_SCHEDULE']  # Frequency of adding food (default = 24 hr)
CULLING_SCHEDULE = param['CULLING_SCHEDULE']  # Frequency of culling (default = 24 hr)
PERCENT_CULL = param['PERCENT_CULL']  # Percent of "media" culled at each culling interval

# Metabolic constants
COST_OF_LIVING = param['COST_OF_LIVING']  # Percent biomass consumed per timestep through metabolism
METABOLIC_EFFICIENCY = param['METABOLIC_EFFICIENCY']  # Percent food converted to worm or egg mass after consumption

# Culling percentages for each stage
EGG_CULL_PERCENT = param['EGG_CULL_PERCENT']
LARVA_CULL_PERCENT = param['LARVA_CULL_PERCENT']
DAUER_CULL_PERCENT = param['DAUER_CULL_PERCENT']
ADULT_CULL_PERCENT = param['ADULT_CULL_PERCENT']
PARLAD_CULL_PERCENT = param['PARLAD_CULL_PERCENT']

# You can now use these constants in your simulation code


# Constants for logistic growth formula:
Kr = 1.78027908103543
Ks = 2.13217438144031
bm = 0.000129574732510846
bn = 0.000248858348410121

# Constants for progeny production formula:
eggFood = [0.07,0.13,0.25,0.50,1.0,4.0]
eggN = [1.319189633,2.201468416,3.484973388,3.640140648,4.725008562,4.156785526]
eggScale = [3.481473101,2.024212254,1.111015285,0.897214743,0.739552921,0.749568986]
eggM = [3.272545277,4.817111781,7.326346126,15.73194284,8.617245216,16.06551326]

# Constants for Gompertz lifespan determination:
gompertzN = 3 # Shape parameter, higher numbers = more square lifespan curve
gompertzLS = 21 * 24 # Roughly average lifespan (days). Equivalent to the 168 timestep value used in the paper.
gompertzA = gompertzLS * (math.exp(gompertzN) - 1)
gompertzTau = 0.85 * (gompertzLS / gompertzN)

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



class Simulation:
    """Totality of the environment

    Hopefully a useful way to keep track of both food and worms. The actual work of the simulation will be run with functions from here.

    In the future, I'll make this executable either as a python script from shell, but currently it's best run interactively.

    Example:
    > ipython
    > import wormpop
    > simulation = wormpop.Simulation(output_location='output') # create simulation object and tell it to place outputs in the directory 'output' in the current working directory
    > simulation.run() # Run simulation with default number of timesteps


    """
    def __init__(self, output_location, number_worms=STARTING_WORMS, starting_stage=STARTING_STAGE, starting_food=STARTING_FOOD, length=SIMULATION_LENGTH, report_individuals=False, connection=None):
        self.worms = Worms()
        self.worms.initialize_worms(number_worms, starting_stage)
        self.dead = Dead_worms()
        self.food = starting_food
        self.food_concentration = self.food / 5 / 1e6
        self.food_history = [self.food_concentration] # Used to keep track of how much food each worm has seen
        self.path = pathlib.Path(output_location)
        self.timestep = 0
        self.time = 0
        self.length = length
        self.report_individuals = report_individuals
        self.connection = connection
        self.bulk_data = []
    
    def iterate_once(self):
        """Meat and potatoes algorithm of the simulation.
        
        At each timestep:
        1) The clock advances/worm age is updated
        2) Culling/adding food
        3) Worm appetite is calculated
        4) Worms eat and grow accordingly
        5) Pay cost of living
        6) Worms undergo checks dependent on their age, stage, and food-availability
        7) Outcome of the timestep is recorded 

        Could mess with the order of this a bit as well. Unclear to me whether worms should pay cost of living "up front" or after eating.
        """
        
        # Advance clock:
        self.timestep += 1
        self.time = self.timestep * TIMESTEP
        
        # Age worms
        self.worms.ageup()

        # Cull/add bacteria, if applicable:
        if self.time % CULLING_SCHEDULE == 0: self.cull(PERCENT_CULL)
        if self.time % FEEDING_SCHEDULE == 0: self.food += FEEDING_AMOUNT

        # Calculate appetite
        self.food_concentration = self.food / 1e6 / 5 # Convert from nanograms to mg/mL
        self.food_history.append(self.food_concentration)
        self.worms.compute_appetite(self.food_concentration) # For simplicity, worms only detect environment once at the start of each time step

        # Feed worms, grow worms
        amount_consumed = self.worms.eat(self.food)
        self.food -= amount_consumed

        # Metabolic upkeep
        self.worms.tax()
        
        # Run Checks
        self.worms.make_checks(self.food_history) # Using food concentration detected before feeding so you're only starving if you didn't get enough to eat

        # Report outcome of timestep
        self.report()


    def cull(self, percent):
        """Periodic culling
        Removes a set percentage of the "media," e.g. 10% of all worms and food to simulate prediation.

        TODO stage specific culling
        
        """

        pct_cull = percent / 100

        self.worms.cull(pct_cull)
        self.food -= self.food * pct_cull

    
    def report(self, header=False):
        """Generates file to keep track of simulation progress.

        Things to keep track of: timestep, hours/days since start, food mass, food concentration, total number of worms, number of each stage, total mass of worms, mass of each stage,
        number dead, causes of death, average age. 
        
        Not yet implemented: average lifespan, mass allocation (growth vs. eggs), rates of transition.

        Realizing this might be faster to just start by building lists/dicts of worms of each stage, rather than iterating over multiple times, but let's see how this does.

        TODO Save individual life histories
        TODO Rates
        TODO Run parameters
        TODO separation of culled worms, non-culled dead worms as food?

        """

        # Individual reporting:

        attributes = [
            'name', 'age', 'stage', 'mass', 'current_egg_progress', 'eggs_laid', 
            'sensed_food', 'appetite', 'growth_mass', 'desired_egg_mass', 
            'actual_egg_mass', 'maintenance', 'portion', 'p_starve', 
            'p_awaken', 'p_death', 'note'
        ]

        # Create the table if it doesn't exist
        create_table_sql = '''
        CREATE TABLE IF NOT EXISTS worms (
            Worm_Name TEXT,
            Timestep INTEGER, 
            Age_hours REAL,
            Stage TEXT,
            Mass REAL,
            Egg_Mass REAL,
            Eggs_Laid INTEGER,
            Available_Food REAL,
            Total_Appetite REAL,
            Desired_Growth REAL,
            Desired_Eggs REAL,
            Actual_Egg_Investment REAL,
            Metabolic_Cost REAL,
            Amount_Eaten REAL,
            Chance_of_Starvation REAL,
            Chance_of_Dauer_Awakening REAL,
            Chance_of_Death REAL,
            Notes TEXT
        )
        '''
        if self.report_individuals:
          
            cursor = self.connection.cursor()
            cursor.execute(create_table_sql)

        
            for w in self.worms:
                w.note = 'Born at timestep {}'.format(self.timestep) if not hasattr(w, 'note') else w.note

                for a in attributes:
                    if not hasattr(w, a):
                        setattr(w, a, '')

                reportlist = [getattr(w, a) for a in attributes]

                # Adding current timestep before the rest of the reportlist attributes
                data_to_insert = [w.name, self.timestep] + reportlist[1:]
                self.bulk_data.append(data_to_insert)
                w.note = ''

            if self.timestep % 10 == 0:
                insert_sql = '''
                INSERT INTO worms (
                    Worm_Name, Timestep, Age_hours, Stage, Mass, Egg_Mass, Eggs_Laid, 
                    Available_Food, Total_Appetite, Desired_Growth, Desired_Eggs, 
                    Actual_Egg_Investment, Metabolic_Cost, Amount_Eaten, 
                    Chance_of_Starvation, Chance_of_Dauer_Awakening, 
                    Chance_of_Death, Notes
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                '''
                cursor.executemany(insert_sql, self.bulk_data)
                self.bulk_data.clear()  # Clear the data list after inserting
                self.connection.commit()
                        

        self.dead.extend([w for w in self.worms if w.stage == 'dead'])
        self.worms[:] = [w for w in self.worms if w.stage != 'dead']

        # Group reporting:
        
        if header:
            with open(self.summary_path, 'w+') as file:
                file.write('\t'.join(['Timestep','Time (hours)','Time (days)', 'Food Mass (ng)', 'Food Conc (mg/mL)', 'Number Worms', 'Number Eggs','Number Larvae', 'Number Dauer',
                'Number Adults', 'Number Parlads','Number Dead','Total Worm Mass (ng)','Egg Mass','Larva Mass','Dauer Mass','Adult Mass','Parlad Mass','Dead Mass','Eggs Laid',
                'Died of old age', 'Died of starvation','Died of bagging','Died of predation','Died of arrested development'])+'\n')

        stages = ['egg','larva','dauer','adult','parlad']
        current_stages = numpy.array([w.stage for w in self.worms])
        stagecounts = [numpy.count_nonzero(current_stages==stage) for stage in stages]
        
        stagemasses = [numpy.sum(numpy.array([w.mass for w in self.worms if w.stage == stage])) for stage in stages]

        n_alive = len(self.worms) # Keeping parlads in the counts for now
        n_dead = len(self.dead)
        mass_alive = numpy.sum(numpy.array([w.mass for w in self.worms]))
        dead_mass = numpy.sum(numpy.array([w.mass for w in self.dead]))

        eggs_laid = numpy.sum(numpy.array([w.eggs_laid for w in self.worms]))

   
        causes_of_death = ['old_age','starvation','bag','culled','arrested_development']

        current_deaths = numpy.array([w.cause_of_death for w in self.dead])
        deathcounts = [numpy.count_nonzero(current_deaths==cause) for cause in causes_of_death]

        #if len(self.dead) > 0: Taking this out for now since it slows things down and isn't that useful
        #    avg_life = numpy.mean(numpy.array([w.lifespan for w in self.dead]))
        #else: 
        #    avg_life = ''
        
        reportlist = [self.timestep, self.time, self.time / 24, self.food_concentration * 5e6, self.food_concentration]

        reportlist.append(n_alive)
        reportlist.extend(stagecounts)
        reportlist.append(n_dead)
        reportlist.append(mass_alive)
        reportlist.extend(stagemasses)
        reportlist.append(dead_mass)
        reportlist.append(eggs_laid)
        reportlist.extend(deathcounts)
        #reportlist.append(avg_life)

        with open(self.summary_path,'a+') as file:
            file.write('\t'.join(map(str, reportlist)) +'\n')
        
        # Report transitions
        egg_to_larva, egg_to_larva_mass = HatchGet()
        larva_to_adult, larva_to_adult_mass = LarvaToAdultGet()
        larva_to_dauer, larva_to_dauer_mass = LarvaToDauerGet()
        adult_to_bag, adult_to_bag_mass = AdultToBagGet()
        dauer_to_larva, dauer_to_larva_mass = DauerToLarvaGet()
        death_metrics = die_ind, die_mass = die_reporter()

        if header:
            with open(self.stage_transition, "w") as fp:
                writer = csv.writer(fp, delimiter="\t")
                writer.writerow([
                    "Timestep",
                    "egg_to_larva", "egg_to_larva_mass", 
                    "larva_to_adult", "larva_to_adult_mass", 
                    "larva_to_dauer","larva_to_dauer_mass",
                    "adult_to_bag","adult_to_bag_mass",
                    "dauer_to_larva", "darva_to_larva_mass"])
            
            with open(self.death_transition, "w") as fp:
                
                writer = csv.writer(fp, delimiter="\t")
                fields = ["Timestep"]
                for _class, value in die_ind.items():
                    for cause_of_death in value.keys():
                        for metric in [ "ind", "mass" ]:
                            fields.append(f"{_class}-{cause_of_death}-{metric}")
                writer.writerow(fields)
                    
                
        with open(self.stage_transition, "a+") as fp:
            writer = csv.writer(fp, delimiter="\t")
            writer.writerow([self.timestep, 
                    egg_to_larva, egg_to_larva_mass,
                    larva_to_adult, larva_to_adult_mass,
                    larva_to_dauer, larva_to_dauer_mass,
                    adult_to_bag, adult_to_bag_mass,
                    dauer_to_larva, dauer_to_larva_mass,
            ])
        
        with open(self.death_transition, "a+") as fp:
            writer = csv.writer(fp, delimiter="\t")
            fields = [self.timestep]
            for _class, value in die_ind.items():
                for cause_of_death in value.keys():
                    for i, _ in enumerate([ "ind", "mass" ]):
                        metric = death_metrics[i]
                        fields.append(metric[_class][cause_of_death])

            writer.writerow(fields)
                    

    def run(self):
        """Run function
        
        Run simulation for set number of timesteps (three hour increments)
        Standard length of 100 days means running for 800 timesteps

        """
        self.path.mkdir(exist_ok=True)
        self.summary_path = self.path / 'summary.tsv'
        self.death_transition = self.path / 'death_transitions.tsv'
        self.stage_transition = self.path / 'stage_transitions.tsv'
    
        if self.report_individuals:
            self.individual_path = self.path / 'indivduals'
            self.individual_path.mkdir(exist_ok=True)
        
        with open(self.path / 'parameters.json', "w") as fp:
            json.dump(param, fp, indent=4)

        self.report(header=True) # Initial conditions/header for output file

        for i in range(1, self.length):
            self.iterate_once()
            if self.timestep % 10 == 0:
                print('{} Timesteps, Food = {} mg/mL, {} Worms Alive, {} Worms Dead'.format(self.timestep, round(self.food_concentration,2), len(self.worms), len(self.dead)))

class Worms(list):
    """Class for holding all worms in the simulation

    Contains methods for things that apply to the whole population, e.g. "Do x to all worms at once."

    Worms object is a list, so it can be iterated through, indexed, and appended to like any other list.
    """

    def initialize_worms(self, number_worms, starting_stage):
        """Using this instead of a standard __init__ function so I can easily build and rebuild list.

        Currently only with identical eggs, larvae, and dauers, but could potentially start with mixed population 
        of randomized ages, masses, etc.
        """
        stagedict = {'egg' : Egg,
                     'larva' : Larva,
                     'dauer' : Dauer}

        assert starting_stage in stagedict, "Only 'egg', 'larva', and 'dauer' may currently be used as starting stage"
        
        for i in range(number_worms):
            name = 'worm_' + str(i + 1)
            self.append(stagedict[starting_stage](name))
        self.total_worm_number = number_worms
    
    def ageup(self):
        [w.ageup() for w in self]

    def tax(self):
        [w.tax() for w in self]

    def cull(self, percent_chance):
        """Each living worm has chance of getting culled at each culling interval.
        """

        for w in self:
            if w.stage == 'dead':
                pass
            else:
                w.cull_maybe()
     
    #@profile
    def compute_appetite(self, food_concentration):
        """Appetite based on growth mass + egg mass + cost of living
        
        Only larvae and adults actually eat and grow, and only adults lay eggs.

        A little confused here since the paper phrases appetite as "the amount of food
        [a worm] would eat if food were plentiful," but both growth mass and egg mass are 
        based on current food availability? Maybe I'm misunderstanding something.

        Probably what is meant by this is something closer to "the amount a worm would eat
        if it had all the food in the environment to itself," which would make sense since 
        a worm presumably has knowledge of the food concentration and its desire to grow,
        but it has less knowledge of how much it will need to share that food (in the model at
        least, since there are still crowd sensing mechanisms in the real world).
        
        """

        for w in self:
            w.sensed_food = food_concentration
            w.get_growth_mass(food_concentration)
            w.get_egg_mass(food_concentration)
            w.get_maintenance()
            w.appetite = (w.growth_mass + w.desired_egg_mass + w.maintenance) / METABOLIC_EFFICIENCY # Previous model only adjusts growth and egg mass by efficiency,
                                                                                                     # so this is a change I am making. Will be good to compare

        self.summed_appetite = numpy.sum(numpy.array([w.appetite for w in self]))


    def eat(self, bacterial_mass):
        """Worms eat as much as they can based on their growth requirements and appetites of other worms.

        Confused about how portion is handled in the previous model, since portion is calculated and then a second restriction:
        (portion*appetite) / (portion + appetite) appears to be applied. I think this is to keep worms from consuming all the 
        available food. If we know empirically that worms grow at a certain rate in a certain concentration, then I think it makes
        the most sense to assume they eat at least that much bacteria, though.
        
        Returns total amount consumed
        """
        if bacterial_mass > self.summed_appetite:
            for w in self:
                w.portion = w.appetite
                w.eat(w.portion)
            return self.summed_appetite
        
        else:
            for w in self:
                w.portion = (w.appetite / self.summed_appetite) * bacterial_mass
                w.eat(w.portion)
            return bacterial_mass

    def make_checks(self, food_history):
        """Runs checks applicable to each worm

        Since this is the only way for new worms to enter the simulation, each check function returns an empty list if there are no new 
        worms, or a list of class objects of the appropriate worm sublcass (Egg or Dauer, for instance). The new arrivals are then appended
        to the Worms object.
        """
        
        current_food = food_history[-1]
        prev_food = food_history[-2]

        new_arrivals = [w.make_checks(current_food, prev_food) for w in self]
        flat_list = [w for new in new_arrivals for w in new]
        
        self += [w('worm_' + str(self.total_worm_number + i + 1)) for i, w in enumerate(flat_list)]
        self.total_worm_number += len(flat_list)



class Dead_worms(list):
    """Testing moving dead worms into this object instead of keeping them with the others to more easily keep track of living worms.

    For the purposes of bookkeeping, parlads are considered "alive" in that they aren't added to this list until they burst. Their lifespan,
    however, is still determined as the moment at which they starve and bag.
    """

    #TODO: decide if this class is worth keeping
    # Might be useful for writing out invdividuals only after they die

    def get_causes_of_death(self):
        """Return a dictionary keyed by cause of death for all dead worms at given timepoint
        """
        causes_of_death = ['old_age','starvation','bag','culled','arrested_development']
        deathcounts = {}
        for cause in causes_of_death:
            deathcounts[cause] = numpy.count_nonzero([w.cause_of_death == cause for w in self if hasattr(w, 'cause_of_death')])

        return deathcounts

    def get_lifespans(self):
        lifespans = [w.lifespan for w in self]
        return lifespans



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
        }

    return {
        "Egg": create_cause_of_death(),
        "Larva": create_cause_of_death(),
        "Adult": create_cause_of_death(),
        "Dauer": create_cause_of_death(),
        "Parlad": create_cause_of_death(),
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


class Worm:
    """Individual in simulation/Parent class for other worm states

    Each individual will behave according to globally defined rules (rates of transition, food availability, etc) when 
    the simulation is run.
    
    Methods in this class are either inherited or overwritten by subclasses. E.g. only larvae and adults need to eat,
    so they get their own methods for calculating appetite, whereas eggs, dauers, parlads, and dead	worms inherit 
    the dummy methods of this parent class.

    Needs parameters for:
    Stage (egg, larva, dauer, parlad (bag), and adult)
    Mass (ng)
    Life history
        Age (hr)
        Transitions (at least for if dauer has already occured)
        Origin? (e.g. born from laid egg or parlad?) * Not currently implemented
        etc

    Each subclass also has its own list of checks to be made at each timestep, e.g. if an egg is ready to hatch or if a
    larva transitions to dauer. After these checks are made, methods for transitions are called if applicable.

    TODO add reporting for individual worms
    TODO add counter for number of transitions called to get wt rates as in previous model

    """

    CULL_PERCENT = 10

    def __init__(self, name):
        self.name = name
    
    def cull_maybe(self):
        roll = random.rand()
        if roll <= self.CULL_PERCENT / 100:
            self.die('culled')
    
    def ageup(self):
        self.age += TIMESTEP

    @die_wrapper
    def die(self, cause_of_death):
        self.__class__ = Dead
        self.__init__(self.name, cause_of_death)

    def tax(self):
        pass

    def get_growth_mass(self, food_concentration):
        self.growth_mass = 0

    def get_egg_mass(self, food_concentration):
        self.desired_egg_mass = 0

    def get_maintenance(self):
        self.maintenance = 0

    def eat(self, amount):
        pass

    def make_checks(self, current_food, prev_food):
        return []


HatchSet, HatchGet = CreateCounter()

class Egg(Worm):
    """First stage
    
    For the sake of simplicity, eggs are considered "worms".
    Eggs are set at a mass of 65 ng by default and hatch after 15 hours (5 timesteps)
    
    """
    def __init__(self, name):
        self.mass = EGGMASS
        self.stage = 'egg'
        self.age = 0
        self.eggs_laid = 0
        self.egg_age = 0 # Eggs hatch after 15 hours, and larvae are born at age 0
        super(Egg, self).__init__(name)

    def ageup(self):
        self.egg_age += TIMESTEP

    def make_checks(self, current_food, prev_food):
        """Only check an egg needs to make is if it's time to hatch
        """
        if self.egg_age >= 15:
            self.hatch()

        return []
    
    @HatchSet
    def hatch(self):
        """After 5 timesteps, an egg becomes a larva
        """
        self.__class__ = Larva
        self.__init__(self.name)

LarvaToDauerSet, LarvaToDauerGet = CreateCounter()
LarvaToAdultSet, LarvaToAdultGet = CreateCounter()

class Larva(Worm):
    """Second stage

    Larvae eat, grow, test their environment to see if they dauer or starve, and potentially become adults after
    a set time and if in a specific mass range.

    """
    def __init__(self, name):
        self.stage = 'larva'
        self.can_dauer = False
        self.eggs_laid = 0
        # The below "if" statements account for situations like if the simulation is being started with larvae,
        # or if a worm is re-entering larvahood after having been a dauer, but still remembers its larval age.
        if not hasattr(self, 'mass'): self.mass = STANDARD_LARVA_MASS
        if not hasattr(self, 'age'): self.age = 0
        if not hasattr(self, 'larval_age'): self.larval_age = 0
        self.p_awaken = ''
        super(Larva, self).__init__(name)

    def ageup(self):
        self.age += TIMESTEP
        self.larval_age += TIMESTEP

    def tax(self):
        self.mass -= self.mass * COST_OF_LIVING

    def get_maintenance(self):
        self.maintenance = self.mass * COST_OF_LIVING

    def get_growth_mass(self, food_conc):
        """Logistic growth formula:

        dx/dt = Kx(1 - bx)

        K = Kr * tanh(Ks * [food])
        b = bm + (bn/[food])
        x = mass
        t = time in days

        Important parameters are current mass and available food.

        Weird discrepancy here between the manuscript and the code: b is "bm - bn/food" in the manuscript vs the "+" above.
        Decided to hew closer to the code, which seems to produce numbers in line with those in the paper.

        """
        if food_conc > 0:
            
            dt = TIMESTEP / 24 # dt = timestep length in days

            K = Kr * math.tanh(Ks * food_conc)
            b = bm + (bn / food_conc)
            dx = K * self.mass * (1 - (b * self.mass)) * dt
            self.growth_mass = dx
        else:
            self.growth_mass = 0

    def eat(self, amount):
        self.mass += amount * METABOLIC_EFFICIENCY

    #@profile
    def make_checks(self, current_food, prev_food):
        """Larvae check if they starve, dauer, advance to adulthood, or fail to hit required adult before max transition age.
        
        dauer threshold defaults = 0.05 mg/mL while in mass range 137-456 ng.
        
        was doing it this way before:

        if food_conc <= 0.05:
            self.starve_count += 1
        else:
            self.starve_count = 0

        if self.starve_count >= 2 and self.can_dauer:
            self.dauer()
        
        Turns out that may be an overly simplistic way of doing it. While I interpreted the paper as saying "< 0.05 mg/mL food for two consecutive timepoints,"
        the code treats it as an exponential probability function. Thus, the probability of a worm to transition to dauer after two timepoints with 0 food is 50%.
        The dauer threshold concentration is just a scaling factor that affects this probability (higher values mean a greater fraction of worms dauering and starving, 
        in other words a lower threshold).
        """

        if self.larval_age <= 3:
            prev_food = current_food # Wasn't born when last food concentration check occurred, only knows current concentration

        if self.mass > MIN_DAUER_MASS and self.mass < MAX_DAUER_MASS and not hasattr(self, 'has_dauered'):
            self.can_dauer = True
        else:
            self.can_dauer = False


        # Check dauer/starvation:
        self.p_starve = (1 / DAUER_RATE) * math.exp(-0.5 * (current_food + prev_food) / DAUER_THRESHOLD)
        roll = random.rand() # Random number between 0 and 1

        if roll < self.p_starve:

            if self.can_dauer:
                self.dauer()
                return []

            else: # Unable to dauer -> starve
                newroll = random.rand()
                if newroll < LARVAL_STARVE_PROB: # Chance of larvae to cheat death
                    self.die('starvation')
                    return []
                else:
                    # Previous code has stipulation here that larave who cheat death lose 1/10th of their body mass, which I may add
                    self.note = 'Cheated death by starvation'
                    pass

        # Check maturity:
        if self.larval_age >= MIN_ADULT_AGE and self.larval_age <= MAX_ADULT_AGE and self.mass >= MIN_ADULT_MASS:
            self.molt()
            return []

        # Larvae die by sticking around too long:
        if self.larval_age > MAX_ADULT_AGE:
            # Before it was determined that worms who never hit adult mass starve if they reach this age. 
            # I guess it's hard to figure out what to do with these worms. Setting a different cause of death 
            # so I can at least see what % of worms die this way.
            self.die('arrested_development')
        
        return []

    @LarvaToDauerSet
    def dauer(self):
        """Enter dauer diapause
        """
        self.__class__ = Dauer
        self.__init__(self.name)
    
    @LarvaToAdultSet
    def molt(self):
        """Mature to adult
        """
        self.__class__ = Adult
        self.__init__(self.name)


DauerToLarvaSet, DauerToLarvaGet = CreateCounter()


class Dauer(Worm):
    """Third stage.

    Dauers don't eat or grow, maintaining the same mass as when they enter dauer.
    
    In the previous model, dauers can eventually die of "attrition," which has an adjustable timescale,
    if they never see enough nutrients to return to a larval state. This is not yet implemented here.
    
    Worms can only enter dauer once.

    TODO alternative diapause states (L1)
    TODO Dauer attrition based on data
    TODO Dauer pheromone

    """
    def __init__(self, name):
        self.stage = 'dauer'
        self.has_dauered = True
        self.eggs_laid = 0
        if not hasattr(self, 'mass'): self.mass = STANDARD_LARVA_MASS
        if not hasattr(self, 'age'): self.age = 15 # assuming 30 hours from parlad bagginng -> 15 hours for eggs to hatch, 15 hours for dauers to develop
        self.p_starve = ''
        super(Dauer, self).__init__(name)

    def make_checks(self, current_food, prev_food):
        """Dauers check if conditions are safe to exit dauer. Currently I'm treating dauers as immportal

        TODO: Add dauer attrition
        """

        self.p_awaken = DAUER_EXIT_PROB * math.sqrt(0.5 * (current_food + prev_food) * 5e6) # Converted back to ng here for convenience
                                                                                            # Could also just use converted dauer exit probability (3.24e-5 * sqrt(5e6) = 0.0724486)
        roll = random.rand()
        if roll < self.p_awaken:
            self.exit_dauer()

        return []
    
    @DauerToLarvaSet
    def exit_dauer(self):
        self.__class__ = Larva
        self.__init__(self.name)


AdultToBagSet, AdultToBagGet = CreateCounter()

class Adult(Worm):
    """Fourth stage

    Adults continue to eat and grow using the same logarithmic function as larvae. However, adults 
    apportion intaken nutrients between somatic and germline mass, which accumulates until an egg can be laid.
    As I currently have it, eggs are laid as soon as germline mass accumulates 65 ng.

    Adults can also starve, at which point they bag and become a "parlad."

    Adult starvation is currently probabilistic like dauer entry for larvae, which is distinct from the previous model which
    used the "two time points below x threshold" method for starvation.

    Finally, adults die of old age according to a gompertz hazard function.

    TODO Fertility span
    TODO Apportion of somatic vs germ mass - fixed value or just based on growth and egg-laying curves?

    """
    def __init__(self, name):
        self.stage = 'adult'
        self.adult_age = 0
        self.total_egg_mass = 0
        self.min_somatic_mass = (self.mass + MIN_ADULT_MASS) / 2
        self.note = 'Min adult mass set to {}'.format(self.min_somatic_mass)
        self.bag_rate = BAG_RATE
        self.bag_threshold = BAG_THRESHOLD
        super(Adult, self).__init__(name)

    def ageup(self):
        self.age += TIMESTEP
        self.adult_age += TIMESTEP

    def tax(self):
        self.mass -= self.mass * COST_OF_LIVING

    def get_maintenance(self):
        self.maintenance = self.mass * COST_OF_LIVING


    def get_growth_mass(self, food_conc):
        """Logistic growth formula:

        Same as in larvae. See above for docstring.

        """

        if food_conc > 0:
            dt = TIMESTEP / 24 # dt = timestep length in days

            K = Kr * math.tanh(Ks * food_conc)
            b = bm + (bn / food_conc)
            dx = K * self.mass * (1 - (b * self.mass)) * dt
            self.growth_mass = dx
        else:
            self.growth_mass = 0

    def get_egg_mass(self, food_conc):
        """Mass of eggs (desired to be) produced on a given day. Used in determining appetite.
        Based on empirical measurement

        eggs/day = M * x^n * exp(-x / x0)

        M = 17.33 * x^(-0.166) + 30.08
        n = -1.86 * x^(-0.299) + 5.908
        x0 = 0.233 * x^(-0.762) + 0.581
        
        x = age in days (?)
        
        ---

        The previous code is written fairly differently to how the above functions are described in the paper,
        so I'll err on the side of the code for now and see if I get reasonable numbers. 

        In the code as written, M, n, and x0 are referenced from a table according to food concentration, from which
        eggs/day is calculated. So I'm unclear on how the 9 constant parameters in the functions above translate to those tables.

        Converted everything here to work in hours and nanograms.
        Updates egg_mass attribute to desired egg mass generated for the current timestep (this can then be added to total egg mass
        to calculate progeny production).

        """

        x = self.adult_age / 24
        food_avail = food_conc

        if food_avail > 1: # Return egg formula for high food conditions
            i = 5
            eggs = eggM[i] * math.exp(eggN[i] * math.log(x) - x/eggScale[i])

        elif food_avail in eggFood: # Edge case if food conc is exactly one of the measured concentrations (mostly for testing)
            i = numpy.where([idx == food_avail for idx in eggFood])[0][0]
            eggs = eggM[i] * math.exp(eggN[i] * math.log(x) - x/eggScale[i])

        else: # Interpolate between egg values for two closest measured food conditions
            if food_avail < eggFood[1]: # Account for values below minimum measured food concentration
                i = 0
                j = 1
            else:
                for idx in range(len(eggFood) - 1): # Choose nearest measured concentrations
                    if food_avail > eggFood[idx] and food_avail < eggFood[idx+1]:
                        i = idx
                        j = idx + 1
            ri = eggM[i] * math.exp(eggN[i] * math.log(x) - x/eggScale[i])
            rj = eggM[j] * math.exp(eggN[j] * math.log(x) - x/eggScale[j])
            p = (food_avail - eggFood[i]) / (eggFood[j] - eggFood[i]) # Weighting factor

            eggs = math.exp(math.log(ri)*(1-p)+math.log(rj)*p) # Interpolate logarithmically between the two values

        self.desired_egg_mass = eggs * EGGMASS

    def eat(self, amount):
        self.mass += amount * METABOLIC_EFFICIENCY
        self.convert_mass()

    def convert_mass(self):
        """Allocates consumed mass to eggsf

        Only allowed to allocate a set amount of mass down to min somatic mass, set as a value between mass at adulthood
        and the minimum possible adult mass. In nutrient rich conditions this shouldn't be an issue, but in times of 
        scarcity this will mimic limitation in progeny production (in addition to the above fertility calculation).

        Total egg mass is then checked at each timestep to see if an egg has been "completed."
        
        The previous code has something called "fertconst," which seems to control amount of energy/mass converted to eggs each timestep.
        I think it probably makes more sense to use the mass of eggs generated based on the curve, since we know that empirically.

        """
        if (self.mass - self.desired_egg_mass) >= self.min_somatic_mass:
            self.actual_egg_mass = self.desired_egg_mass
        elif self.mass > self.min_somatic_mass:
            self.actual_egg_mass = self.mass - self.min_somatic_mass
            self.note = 'Dipped into fat stores'
        else:
            self.actual_egg_mass = 0

        self.mass -= self.actual_egg_mass
        self.total_egg_mass += self.actual_egg_mass

    #@profile
    def make_checks(self, current_food, prev_food):
        """Adults check if they are ready to lay an egg, if they starve and turn into a parlad, or if they die of old age.

        Interestingly, the code for bagging in the previous model is deterministic (like how I had dauer entry coded above),
        in that two timesteps below the threshold automatically triggers bagging. I wonder why it was decided not to make this
        probabilistic in the same way that dauering and larval starving is? I think I'm going to code it as probabilistic for now. 
        
        Not entirely sure how to determine fertility yet, so I'm currently treating any adults that starve as parlads, despite the
        fact that many of them will likely not still be in an egg-producing mode... Will circle back to this.

        #TODO try deterministic vs probabilistic model

        """

        eggs = []

        eggs_available = self.total_egg_mass // EGGMASS
        self.current_egg_progress = self.total_egg_mass - (self.eggs_laid * EGGMASS)

        if eggs_available > self.eggs_laid: # Check if enough egg mass has been added to lay a new egg
            new_eggs = self.lay_egg(eggs_available - self.eggs_laid)
            eggs.extend(new_eggs)

        self.p_starve = (1 / self.bag_rate) * math.exp(-0.5 * (current_food + prev_food) / self.bag_threshold)
        roll = random.rand()
        if roll < self.p_starve:
            self.bag()
            return eggs

        self.p_death = (math.exp(self.age / gompertzTau) - 1) / gompertzA # Probability of dying at given time
        roll = random.rand()
        if roll < self.p_death:
            self.die('old_age')

        return eggs

    def lay_egg(self, number):
        self.eggs_laid += number
        return [Egg] * int(number)

    @AdultToBagSet
    def bag(self):
        self.__class__ = Parlad
        self.__init__(self.name)


class Parlad(Worm):
    """Fifth stage 

    Special kind of death resulting in a bag of worms that bursts into dauers after 30 hours.
    Number of dauers generated is based on mass, reduced by an efficiency parameter (.66 by default)

    """
    def __init__(self, name):
        self.stage = 'parlad'
        self.lifespan = self.age
        self.cause_of_death = ('bag')
        self.mass += self.total_egg_mass - (self.eggs_laid * EGGMASS) # Total mass (unlaid eggs + somatic mass)
        self.dauer_potential = int((self.mass * BAG_EFFICIENCY) // STANDARD_LARVA_MASS)
        self.mass_decrement = self.mass / (30 / TIMESTEP) # Will lose this much mass per timestep (converted into dauers) down to 0
        self.note = 'Will burst into {} dauers in 30 hours'.format(self.dauer_potential)
        super(Parlad, self).__init__(name)

    def tax(self):
        """Rather than paying "metabolic tax," going to use this to keep track of mass as it is consumed by matricidal hatching
        """

    def make_checks(self, current_food, prev_food):
        """Parlads check if it's time to burst
        """
        released_dauers = []
        if self.age - self.lifespan >= 30:
            released_dauers.extend([Dauer]*self.dauer_potential)
            self.mass = self.mass - self.dauer_potential * STANDARD_LARVA_MASS
            assert self.mass > 0
            self.die('bag')

        return released_dauers


def create_dead_mass_counter():
    
    mass_counter = 0

    def init_wrapper(func):
        def wraps(self, *args, **kwargs):
            nonlocal mass_counter
            func(self, *args, **kwargs)
            mass_counter += self.mass
        return wraps

    def get_dead_mass():
        pass
    


class Dead(Worm):
    """Final stage

    Catch-all subclass to make sure dead worms don't keep doing the activity of the living.

    Now dead worms are moved into their own object (Dead_worms) after their death is recorded.
    """
    def __init__(self, name, cause_of_death):
        self.stage = 'dead'
        if not hasattr(self, 'cause_of_death'): self.cause_of_death =  cause_of_death
        if not hasattr(self, 'lifespan'): self.lifespan = self.age
        self.note = 'Lifespan: {} days, Cause of death: {}'.format(self.lifespan / 24, self.cause_of_death)
        super(Dead, self).__init__(name)

#%%




Egg.CULL_PERCENT = EGG_CULL_PERCENT
Larva.CULL_PERCENT = LARVA_CULL_PERCENT
Dauer.CULL_PERCENT = DAUER_CULL_PERCENT
Adult.CULL_PERCENT = ADULT_CULL_PERCENT
Parlad.CULL_PERCENT = PARLAD_CULL_PERCENT


import datetime

if args["--database"]:
    with sqlite3.connect(args["--database"]) as conn:
        conn.execute('''CREATE TABLE IF NOT EXISTS Metadata
                    (name TEXT, start_time TEXT, parameter TEXT, status TEXT, error TEXT)''')
        # Inserting values into the Metadata table
        name = 'speedtest'
        start_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        conn.execute('INSERT INTO Metadata (name, start_time, parameter) VALUES (?, ?, ?)', (name, start_time, json.dumps(param)))
        # Committing the transaction
        conn.commit()
        test = Simulation(directory, connection=conn, report_individuals=True)
        test.run()
else:
    test = Simulation(directory)
    test.run()
