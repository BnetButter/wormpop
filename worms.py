import random
import numpy as np
from reporting import CreateCounter, CreateDeathCounter
from database import WormSummary
import math
from genome import Genome
import typing
import functools

from config import load_constants
constants = load_constants()

import simulation_globals


# Simulation time details
SIMULATION_LENGTH = constants['SIMULATION_LENGTH']  # 800 timesteps = 100 days
TIMESTEP = constants['TIMESTEP']  # 3 hr per timestep, 8 timesteps per day

# Initial conditions
STARTING_WORMS = constants['STARTING_WORMS']
STARTING_STAGE = constants['STARTING_STAGE']  # 'egg'
EGGMASS = constants['EGGMASS']  # Nanograms

# Adult constants
MIN_ADULT_MASS = constants['MIN_ADULT_MASS']  # Minimum mass to be an adult (default 800 ng)
MIN_ADULT_AGE = constants['MIN_ADULT_AGE']  # Minimum age in hours to transition to adult
MAX_ADULT_AGE = constants['MAX_ADULT_AGE']  # Max age in hours to transition (larvae past this age die of "arrested development")

# Larva constants
STANDARD_LARVA_MASS = constants['STANDARD_LARVA_MASS']  # ~228 ng
LARVAL_STARVE_PROB = constants['LARVAL_STARVE_PROB']  # Chance to cheat death by starvation, though starvation is probabilistic

# Dauer constants
MIN_DAUER_MASS = constants['MIN_DAUER_MASS']  # ~137 ng
MAX_DAUER_MASS = constants['MAX_DAUER_MASS']  # ~456 ng
DAUER_THRESHOLD = constants['DAUER_THRESHOLD']  # Concentration (mg/mL) that scales probability of dauering = 250,000 ng total food available
DAUER_RATE = constants['DAUER_RATE']  # Number of days at 0 food concentration for a larva to have a 50% chance of dauering/starving
DAUER_EXIT_PROB = constants['DAUER_EXIT_PROB']  # Chance per timestep to exit dauer, based on empirical data

# Bag constants
BAG_THRESHOLD = constants['BAG_THRESHOLD']  # mg/mL (=2500 ng)
BAG_RATE = constants['BAG_RATE']
BAG_EFFICIENCY = constants['BAG_EFFICIENCY']  # Efficiency with which somatic mass of parlads can be converted to dauers

# Food constants
STARTING_FOOD = constants['STARTING_FOOD']  # 10 mg = 1x10^7 ng
FEEDING_AMOUNT = constants['FEEDING_AMOUNT']  # 10 mg added per feeding schedule

# Environment size
FLASK_VOLUME = constants['FLASK_VOLUME'] # default = 5 mL

# Scheduling constants
FEEDING_SCHEDULE = constants['FEEDING_SCHEDULE']  # Frequency of adding food (default = 24 hr)
CULLING_SCHEDULE = constants['CULLING_SCHEDULE']  # Frequency of culling (default = 24 hr)
PERCENT_CULL = constants['PERCENT_CULL']  # Percent of "media" culled at each culling interval

# Metabolic constants
COST_OF_LIVING = constants['COST_OF_LIVING']  # Percent biomass consumed per timestep through metabolism
METABOLIC_EFFICIENCY = constants['METABOLIC_EFFICIENCY']  # Percent food converted to worm or egg mass after consumption

# Culling percentages for each stage
EGG_CULL_PERCENT = constants['EGG_CULL_PERCENT']
LARVA_CULL_PERCENT = constants['LARVA_CULL_PERCENT']
DAUER_CULL_PERCENT = constants['DAUER_CULL_PERCENT']
ADULT_CULL_PERCENT = constants['ADULT_CULL_PERCENT']
PARLAD_CULL_PERCENT = constants['PARLAD_CULL_PERCENT']


L1_ARREST_ENTER_THRESHOLD = constants["L1_ARREST_ENTER_THRESHOLD"]
L1_ARREST_EXIT_THRESHOLD = constants["L1_ARREST_EXIT_THRESHOLD"]

# You can now use these constants in your simulation code
GENOME_VERSION = "0.1"

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

def starve_from_l1_arrest(num_days):
    L = 98.09173354410315
    x0 = 17.484181571141384
    k = -0.47620124366404964

        # Define the reverse sigmoid function
    def reverse_sigmoid(x, L, x0, k):
        return L / (1 + np.exp(-k * (x - x0)))

    return reverse_sigmoid(num_days, L, x0, k)

def egg_curve(x, Y, genome: "Genome"):
    Y = Y/4
    dt = TIMESTEP / 24 # dt = timestep length in days

    eggN = genome.eggN
    eggM = genome.eggM
    eggScale = genome.eggScale

    def f1(x):
        return (eggM * 3.273) * x**(eggN * 1.319) * np.exp(-x / (eggScale * 3.481)) * dt

    def f2(x):
        return (eggM * 4.817) * x**(eggN * 2.201) * np.exp(-x / (eggScale * 2.024)) * dt

    def f3(x):
        return (eggM * 7.326) * x**(eggN * 3.485) * np.exp(-x / (eggScale * 1.111)) * dt

    def f4(x):
        return (eggM * 16.86) * x**(eggN  * 4.157) * np.exp(-x / (eggScale * 0.75)) * dt

    if 0 <= Y < 1/3:
        return (1 - 3*Y) * f1(x) + 3*Y * f2(x)
    elif 1/3 <= Y < 2/3:
        return (2 - 3*Y) * f2(x) + (3*Y - 1) * f3(x)
    elif 2/3 <= Y <= 1:
        return (3 - 3*Y) * f3(x) + (3*Y - 2) * f4(x)
    else:
        return f4(x)

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

    genome: Genome # Just a type hint
    can_dauer = False

    can_arrest = False
    has_arrested = False
    # For Reporting
    _metabolic_efficiency_loss = 0
    _actual_growth = 0
    _eggs_laid_timestep = 0

    _age_hours_entered_larva = 0
    _age_hours_entered_adult = 0
    _age_hours_entered_larva_after_dauer = 0
    _age_hours_entered_parlad = 0
    _age_hours_entered_dauer = 0

    _age_hours_entered_l1 = 0
    _age_hours_entered_larva_after_l1 = 0
    
    _eggs_laid = 0
    
    _summary_table: typing.Optional[WormSummary] = None
    

    def __init__(self, name, genome=None):
        self.name = name
        self._summary_table = WormSummary(self.name) if self._summary_table is None else self._summary_table

        # self.genome = random.sample([NormalAppetite, FatWorm, SkinnyWorm])
        choices = simulation_globals.variants

        #randomly pick a choice if choices isn't decided yet
        if not hasattr(self, "genome"):
            self.genome = genome if genome else random.choice(choices)

    def cull_maybe(self):
        roll = np.random.rand()
        if roll <= self.CULL_PERCENT / 100:
            self.die('culled')
    
    def ageup(self):
        self.age += TIMESTEP

    @die_wrapper
    def die(self, cause_of_death):
        if cause_of_death == "bag":
            self._summary_table.Parlad_span_days = (self.age - self._age_hours_entered_parlad) / 24
        
        self._summary_table.Cause_of_Death = cause_of_death
        self._summary_table.Total_Body_Mass = self.mass
        self._summary_table.Total_Eggs_Laid = getattr(self, "eggs_laid", 0)

        if cause_of_death == "culled":
            simulation_globals.instance._summary.Worms_died_cull += 1
        elif cause_of_death == "starvation":
            simulation_globals.instance._summary.Worms_died_starvation += 1
        elif cause_of_death == "bag":
            simulation_globals.instance._summary.Worms_died_bag += 1
        elif cause_of_death == "old_age":
            simulation_globals.instance._summary.Worms_died_old_age += 1
        elif cause_of_death == "arrested_development":
            simulation_globals.instance._summary.Worms_died_arrested_development += 1



 



        if type(self) == Egg:
            pass
        elif type(self) == Larva:
            if self.has_arrested:
                self._summary_table.LarvaPost_L1_span_days = (self.age - self._age_hours_entered_larva_after_l1) / 24
            else:
                self._summary_table.Larva_span_days = (self.age - self._age_hours_entered_larva) / 24

            if not getattr(self, "has_dauered", False):
                self._summary_table.Larva_span_days = (self.age - self._age_hours_entered_larva) / 24
            else:
                self._summary_table.LarvaPostDauer_span_days = (self.age - self._age_hours_entered_larva_after_dauer) / 24

        elif type(self) == L1Arrest:
            self._summary_table.L1_span_days = (self.age - self._age_hours_entered_l1) / 24
        elif type(self) == Dauer:
            self._summary_table.Dauer_span_days = (self.age - self._age_hours_entered_dauer) / 24
        elif type(self) == Adult:
            self._summary_table.Adult_span_days = (self.age - self._age_hours_entered_adult) / 24
        elif type(self) == Parlad:
            self._summary_table.Parlad_span_days = (self.age - self._age_hours_entered_parlad) / 24
        
        self._summary_table.Life_span_days = self.age / 24

        egg_history = getattr(self, "egg_history", [])
        if egg_history:
            tstart = 0
            tend = 0
            for timestep, egg in egg_history:
                if egg:
                    tstart = timestep
                    break 
            for timestep, egg in reversed(egg_history):
                if egg:
                    tend = timestep
                    break
            self._summary_table.Reproductive_Span = TIMESTEP*(tend - tstart) / 24
        else:
            pass
            
        self._summary_table.Total_Mass = (
            self._summary_table.Total_Body_Mass 
            + self._summary_table.Total_Egg_Mass 
            + self._summary_table.Total_Metabolic_Cost 
            + self._summary_table.Total_Metabolic_Tax 
            + (0 if self._summary_table.Parlad_Mass_Converted is None else self._summary_table.Parlad_Mass_Converted)
        )

        self.__class__ = Dead
        self.__init__(self.name, cause_of_death)


        # Cache the summary table to be committed later
        simulation_globals.instance.bulk_data.append(self._summary_table)

            

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
    def __init__(self, name, *args, **kwargs):
        super().__init__(name, *args, **kwargs)
        self.mass = EGGMASS
        self.stage = 'egg'
        self.age = 0
        self.eggs_laid = 0
        self.egg_age = 0 # Eggs hatch after 15 hours, and larvae are born at age 0
      

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
        self._age_hours_entered_larva = self.age
        self.__class__ = Larva
        self.__init__(self.name)

LarvaToDauerSet, LarvaToDauerGet = CreateCounter()
LarvaToAdultSet, LarvaToAdultGet = CreateCounter()

LarvaToL1ArrestSet, LarvaToL1ArrestGet = CreateCounter()
L1ArrestToLarvaSet, L1ArrestToLarvaGet = CreateCounter()

class Larva(Worm):
    """Second stage

    Larvae eat, grow, test their environment to see if they dauer or starve, and potentially become adults after
    a set time and if in a specific mass range.

    """
    def __init__(self, name, *args, **kwargs):
        self.stage = 'larva'
        self.can_dauer = False
        self.eggs_laid = 0
        # The below "if" statements account for situations like if the simulation is being started with larvae,
        # or if a worm is re-entering larvahood after having been a dauer, but still remembers its larval age.
        if not hasattr(self, 'mass'): self.mass = STANDARD_LARVA_MASS
        if not hasattr(self, 'age'): self.age = 0
        if not hasattr(self, 'larval_age'): self.larval_age = 0
        self.p_awaken = None
        super(Larva, self).__init__(name)

    @LarvaToL1ArrestSet
    def l1_arrest(self):
        self._age_hours_entered_l1 = self.age
        self.__class__ = L1Arrest
        self.__init__(self.name)

    def ageup(self):
        self.age += TIMESTEP
        self.larval_age += TIMESTEP

    def tax(self):
        metabolic_tax = self.mass * self.genome.metabolic_tax
        self._summary_table.Total_Metabolic_Tax += metabolic_tax
        self.mass -= metabolic_tax
        assert self.mass > 0, "tax"

    def get_maintenance(self):
        self.maintenance = self.mass * self.genome.metabolic_tax

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
            
            self.growth_mass = dx if dx > 0 else 0
        else:
            self.growth_mass = 0
        
        
        assert self.growth_mass >= 0

    def eat(self, amount):
        metabolic_efficiency_loss = amount - (amount * METABOLIC_EFFICIENCY)

        simulation_globals.instance._summary.Bacteria_to_worm_ingested_mg += amount/1e6
        simulation_globals.instance._summary.Bacteria_to_worm_metabolic_inefficiency_mg += metabolic_efficiency_loss / 1e6
        simulation_globals.instance._summary.Bacteria_to_worm_somatic_mass_mg = amount * METABOLIC_EFFICIENCY / 1e6
        


        self._summary_table.Total_Food_Consumed += amount
        self._summary_table.Total_Metabolic_Cost += metabolic_efficiency_loss

        self._metabolic_efficiency_loss = metabolic_efficiency_loss
        curr_mass = self.mass
        self.mass += amount * METABOLIC_EFFICIENCY
        self._actual_growth = self.mass - curr_mass

    @staticmethod
    def dauer_pheremone_function(num_worms):
        num_worms = numpy.array(num_worms)

        # Normalize the number of worms
        if len(num_worms) == 0:
            return 0
        else:
            mu = numpy.mean(num_worms)
            sigma = numpy.std(num_worms)
            X_standardized = (num_worms - mu) / sigma
            # A logistic function
            return numpy.mean(1 / (1 + numpy.exp(-X_standardized)))

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

        if self.mass < MIN_DAUER_MASS and not self.has_arrested and not hasattr(self, 'has_dauered'):
            if current_food + prev_food < L1_ARREST_ENTER_THRESHOLD:
                self.l1_arrest()
                return []

        if self.can_dauer:
            if np.random.rand() < self.genome.dauer_probability:
                self.dauer()
                return []

        dauer_multiplier = 1 #self.dauer_pheremone_function(simulation_globals.instance.worm_count)


        # Check dauer/starvation:
        self.p_starve = (1 / DAUER_RATE) * math.exp(-0.5 * (current_food + prev_food) / DAUER_THRESHOLD) * dauer_multiplier
        
        

        roll = np.random.rand() # Random number between 0 and 1


        if roll < self.p_starve:

            if self.can_dauer:
                self.dauer()
                return []

            else: # Unable to dauer -> starve
                newroll = np.random.rand()
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
        self._summary_table.Larva_span_days = (self.age - self._age_hours_entered_larva) / 24
        self._age_hours_entered_dauer = self.age
        self.__class__ = Dauer
        self.__init__(self.name)
    
    @LarvaToAdultSet
    def molt(self):
        """
        Mature to adult
        """
        if not getattr(self, 'has_dauered', False):
            self._summary_table.Larva_span_days = (self.age - self._age_hours_entered_larva) / 24
        else:
            self._summary_table.LarvaPostDauer_span_days = (self.age - self._age_hours_entered_larva_after_dauer) / 24

        self._age_hours_entered_adult = self.age
        self.__class__ = Adult
        self.__init__(self.name)




class L1Arrest(Worm):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stage = 'L1_Arrest'
        self.has_arrested = True
        assert not getattr(self, "has_dauered", False), "L1s cannot have dauered"

        self._summary_table.Larva_span_days = (self.age - self._age_hours_entered_larva) / 24

    def make_checks(self, current_food, prev_food):
        probability = starve_from_l1_arrest(simulation_globals.instance.timestep / 24)

        if np.random.random() > probability:
            self.die("starvation")
            return []
    
        if current_food + prev_food > L1_ARREST_EXIT_THRESHOLD:
            self.exit_arrest()
            return []
 
        return []
        
    @L1ArrestToLarvaSet
    def exit_arrest(self):
        self._summary_table.LarvaPost_L1_span_days = (self.age - self._age_hours_entered_l1) / 24
        self._age_hours_entered_larva_after_l1 = self.age
        self.__class__ = Larva
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
        self.p_starve = None
        super(Dauer, self).__init__(name)

    def make_checks(self, current_food, prev_food):
        """Dauers check if conditions are safe to exit dauer. Currently I'm treating dauers as immortal

        TODO: Add dauer attrition
        """

        self.p_awaken = DAUER_EXIT_PROB * math.sqrt(0.5 * (current_food + prev_food) * 1e6 * FLASK_VOLUME) # Converted back to ng here for convenience
                                                                                                           # Could also just use converted dauer exit probability (3.24e-5 * sqrt(5e6) = 0.0724486)
        roll = np.random.rand()
        if roll < self.p_awaken:
            self.exit_dauer()

        return []
    
    @DauerToLarvaSet
    def exit_dauer(self):
        self._summary_table.Dauer_span_days = (self.age - self._age_hours_entered_dauer) / 24
        self._age_hours_entered_larva_after_dauer = self.age
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
        self.egg_history: List[Tuple[int, int]] = []
        super(Adult, self).__init__(name)

    def ageup(self):
        self.age += TIMESTEP
        self.adult_age += TIMESTEP

    def tax(self):
        tax = self.mass * self.genome.metabolic_tax
        simulation_globals.instance._summary.Bacteria_to_worm_metabolic_tax_mg += tax / 1e6

        self._summary_table.Total_Metabolic_Tax += tax
        self.mass -= tax

    def get_maintenance(self):
        self.maintenance = self.mass * self.genome.metabolic_tax


    def get_growth_mass(self, food_conc):
        """Logistic growth formula:

        Same as in larvae. See above for docstring.

        """

        if food_conc > 0:
            
            dt = TIMESTEP / 24 # dt = timestep length in days

            K = Kr * math.tanh(Ks * food_conc)
            b = bm + (bn / food_conc)

       
            dx = K * self.mass * (1 - (b * self.mass)) * dt
            
            self.growth_mass = dx if dx > 0 else 0
        else:
            self.growth_mass = 0
        
        
        assert self.growth_mass >= 0
        
        # assert self.growth_mass >= 0

    def get_egg_mass(self, food_conc):
        """Mass of eggs (desired to be) produced on a given day. Used in determining appetite.
        Based on empirical measurement

        d_eggs/dt = M * x^n * exp(-x / x0)

        x = adult age (days)

        M, n, and x0 are referenced from the above eggM, eggN, and eggScale lists according to food concentration.

        To get eggs per timestep (rather than eggs per day), the result is multiplied by timestep length / 24

        Returns desired mass to be allocated to egg production, in nanograms


        """

        x = self.adult_age / 24 # Convert adult age from hours to days
        dt = TIMESTEP / 24
        food_avail = food_conc

        eggs = egg_curve(x, food_avail, self.genome)
        self.desired_egg_mass = eggs * EGGMASS

    def eat(self, amount):
        metabolic_efficiency_loss = amount - (amount * METABOLIC_EFFICIENCY)

        simulation_globals.instance._summary.Bacteria_to_worm_ingested_mg += amount/1e6
        simulation_globals.instance._summary.Bacteria_to_worm_metabolic_inefficiency_mg += metabolic_efficiency_loss/1e6
        simulation_globals.instance._summary.Bacteria_to_worm_somatic_mass_mg = amount * METABOLIC_EFFICIENCY / 1e6

        self._summary_table.Total_Food_Consumed += amount
        self._summary_table.Total_Metabolic_Cost += metabolic_efficiency_loss
        
        self._metabolic_efficiency_loss = metabolic_efficiency_loss

        curr_mass = self.mass
        self.mass += amount * METABOLIC_EFFICIENCY
        self.convert_mass()
        self._actual_growth = self.mass - curr_mass

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
        else:
            self.actual_egg_mass = 0
        
        simulation_globals.instance._summary.Bacteria_to_worm_eggs_mg += self.actual_egg_mass / 1e6


        self.mass -= self.actual_egg_mass
        assert self.mass > 0, "mass > actual_egg_mass"
        self.total_egg_mass += self.actual_egg_mass
        self._summary_table.Total_Egg_Mass += self.actual_egg_mass

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
        roll = np.random.rand()
        if roll < self.p_starve:
            self.bag()
            return eggs

        self.p_death = self.genome.life_span * (math.exp(self.age / gompertzTau) - 1) / gompertzA # Probability of dying at given time

        
        roll = np.random.rand()
        if roll < self.p_death:
            self.die('old_age')

        return eggs

    def lay_egg(self, number):
        self.eggs_laid += number

        self.egg_history.append((
            number,
            simulation_globals.instance.timestep
        ))

        Egg_partial = functools.partial(Egg, genome=self.genome)
        self._eggs_laid_timestep = int(number)
        simulation_globals.instance._summary.Worms_born_egg = int(number)
        result = [Egg_partial] * int(number)
        if result:
            simulation_globals.instance.worms_that_laid_eggs.add(self.name)
        return result
    @AdultToBagSet
    def bag(self):
        self._summary_table.Adult_span_days = (self.age - self._age_hours_entered_adult) / 24
        self._age_hours_entered_parlad = self.age
        self.__class__ = Parlad
        self.__init__(self.name)

def CountParladToDauer():
    num_parlads = 0

    def set_parlads(func):
        def wraps(*args, **kwargs):
            nonlocal num_parlads
            val = func(*args, **kwargs)
            num_parlads += len(val)
            return val
        return wraps
    
    def get_parlads():
        nonlocal num_parlads
        val = num_parlads
        num_parlads = 0
        return val, val * STANDARD_LARVA_MASS
    
    return get_parlads, set_parlads


ParladToDauerGet, ParladToDauerSet = CountParladToDauer()
    

    

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
        assert self.mass > 0, f"{self.total_egg_mass - (self.eggs_laid * EGGMASS)}, {self.genome}"

    def tax(self):
        """Rather than paying "metabolic tax," going to use this to keep track of mass as it is consumed by matricidal hatching
        """

    @ParladToDauerSet
    def make_checks(self, current_food, prev_food):
        """Parlads check if it's time to burst
        """
        assert self.mass > 0
        released_dauers = []
        if self.age - self.lifespan >= 30:
            
            simulation_globals.instance._summary.Worms_born_dauer += self.dauer_potential            

            released_dauers.extend([Dauer]*self.dauer_potential)

            new_mass = self.mass - self.dauer_potential * STANDARD_LARVA_MASS
            self._summary_table.Total_Body_Mass -= new_mass
            self._summary_table.Parlad_Mass_Converted = self.dauer_potential * STANDARD_LARVA_MASS
            assert self.mass > 0, f"{self.dauer_potential} dauers released from {self.name} but mass is {self.mass}"
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

# TODO - Change Egg Mass efficiency
# TODO - Change Life Span
# TODO - Parameterize the genome


Egg.CULL_PERCENT = EGG_CULL_PERCENT
Larva.CULL_PERCENT = LARVA_CULL_PERCENT
Dauer.CULL_PERCENT = DAUER_CULL_PERCENT
Adult.CULL_PERCENT = ADULT_CULL_PERCENT
Parlad.CULL_PERCENT = PARLAD_CULL_PERCENT