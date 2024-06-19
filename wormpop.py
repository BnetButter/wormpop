#@Author: Matt Mosley
#2021-10-19
#%%

"""
USAGE: wormpop [--parameters=<string>] [ --database=<string> ] [ --name=<string> ] [--directory=<string>] [ --variants=<string> ] [ --report-individuals ] [ --socket=<string> ]
"""

import pathlib
import math
import numpy
from numpy import random
import csv
import json
import docopt
import collections
import functools
import os
import math
import numpy as np
import websockets
import sys
from websockets.server import WebSocketServerProtocol

from typing import *

from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    Float,
    ForeignKey
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import (
    sessionmaker,
    relationship
)

from sqlalchemy.sql import expression
from sqlalchemy.schema import DefaultClause


Base = declarative_base()

args = docopt.docopt(__doc__)
parameters = args["--parameters"]
database = args["--database"] if args["--database"] is not None else ":memory:"
name = args["--name"] if args["--name"] is not None else "Simulation"
directory = args["--directory"] if args["--directory"] is not None else "Simulation"
report_individuals = args["--report-individuals"]

# Constants:

import json

# Read the JSON file
if not parameters:
    parameters = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "constants.json"
    )

with open(parameters, 'r') as file:
    param = json.load(file)

        
# Simulation time details
SIMULATION_LENGTH = param['SIMULATION_LENGTH']  # 800 timesteps = 100 days
TIMESTEP = param['TIMESTEP']  # 3 hr per timestep, 8 timesteps per day

# Initial conditions
STARTING_WORMS = param['STARTING_WORMS']
STARTING_STAGE = param['STARTING_STAGE']  # 'egg'
EGGMASS = param['EGGMASS']  # Nanograms

# Adult constants
MIN_ADULT_MASS = param['MIN_ADULT_MASS']  # Minimum mass to be an adult (default 800 ng)
MIN_ADULT_AGE = param['MIN_ADULT_AGE']  # Minimum age in hours to transition to adult
MAX_ADULT_AGE = param['MAX_ADULT_AGE']  # Max age in hours to transition (larvae past this age die of "arrested development")

# Larva constants
STANDARD_LARVA_MASS = param['STANDARD_LARVA_MASS']  # ~228 ng
LARVAL_STARVE_PROB = param['LARVAL_STARVE_PROB']  # Chance to cheat death by starvation, though starvation is probabilistic

# Dauer constants
MIN_DAUER_MASS = param['MIN_DAUER_MASS']  # ~137 ng
MAX_DAUER_MASS = param['MAX_DAUER_MASS']  # ~456 ng


#DAUER_THRESHOLD = param['DAUER_THRESHOLD']  # Concentration (mg/mL) that scales probability of dauering = 250,000 ng total food availablea
#DAUER_RATE = param['DAUER_RATE']  # Number of days at 0 food concentration for a larva to have a 50% chance of dauering/starving

DAUER_EXIT_PROB = param['DAUER_EXIT_PROB']  # Chance per timestep to exit dauer, based on empirical data

DOUBLE_FEED_INTERVAL_DAYS = param["DOUBLE_FEED_INTERVAL_DAYS"]

# Bag constants
BAG_THRESHOLD = param['BAG_THRESHOLD']  # mg/mL (=2500 ng)
BAG_RATE = param['BAG_RATE']
BAG_EFFICIENCY = param['BAG_EFFICIENCY']  # Efficiency with which somatic mass of parlads can be converted to dauers

# Food constants
STARTING_FOOD = param['STARTING_FOOD']  # 10 mg = 1x10^7 ng
FEEDING_AMOUNT = param['FEEDING_AMOUNT']  # 10 mg added per feeding schedule

# Environment size
FLASK_VOLUME = param['FLASK_VOLUME'] # default = 5 mL

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



WORLD_MAX_X = param["WORLD_MAX_X"]
WORLD_MAX_Y = param["WORLD_MAX_Y"]

def normal():
    normal_distribution = np.random.normal()
    # Clip the distribution to the range [-1, 1]
    clipped_distribution = np.clip(normal_distribution, -1, 1)
    return clipped_distribution

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

def get_column_default(column):
    if column.default is None:
        return None
    if isinstance(column.default, DefaultClause):
        if isinstance(column.default.arg, expression.Function):
            return None
        return column.default.arg
    return column.default.arg


def distance(a, b):
    return math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

def move_towards_food(agent_pos, food_pos, speed):
    # Extract coordinates
    x, y = agent_pos
    x_f, y_f = food_pos
    
    # Calculate the Euclidean distance
    distance = math.sqrt((x_f - x)**2 + (y_f - y)**2)
    
    # Calculate the unit vector components
    unit_vector_x = (x_f - x) / distance
    unit_vector_y = (y_f - y) / distance
    
    # Calculate the movement vector components
    move_x = speed * unit_vector_x
    move_y = speed * unit_vector_y
    
    # Update the agent's position

    return move_x, move_y


class Food:
    """
    Represents a food lawn
    """

    def __init__(self, id, pos, amount):
        self.id = id
        self.pos = pos
        self.amount = amount
        self.my_worms = {}
        self.summed_appetite = 0
        self.food_history = [self.food_concentration]
    
    @property
    def food_concentration(self):
        return self.amount / 1e6 / FLASK_VOLUME

    @property
    def radius(self):
        return math.sqrt(self.amount / FLASK_VOLUME / math.pi) / 2
    
class FoodSQL(Base):
    __tablename__ = "Food"
    id = Column(Integer, primary_key=True)

    food_id = Column(Integer)
    timestep = Column(Integer)
    x = Column(Integer)
    y = Column(Integer)
    amount = Column(Float)
    radius = Column(Integer)
    summed_appetite = Column(Float)

    @classmethod
    def get_schema(cls):
        schema = {}
        for column in cls.__table__.columns:
            column_type = str(column.type)
            if column_type.startswith("VARCHAR") or column_type.startswith("STRING"):
                column_type = "string"
            elif column_type.startswith("FLOAT"):
                column_type = "float"
            elif column_type.startswith("INTEGER"):
                column_type = "integer"
           
            default_value = get_column_default(column)

            schema[column.name] = {
                "type": column_type,
                "default": default_value
            }

        return schema


class Genome(Base):
    __tablename__ = "Genome"
    variant = Column(String, primary_key=True)

    appetite: float = Column(Float, default=1)
    life_span: float = Column(Float, default=1)
    metabolic_tax: float = Column(Float, default=0.035)

    eggN: float = Column(Float, default=1)
    eggM: float = Column(Float, default=1)
    eggScale: float = Column(Float, default=1)

    # # Number of days at 0 food concentration for a larva to have a 50% chance of dauering/starving
    dauer_rate = Column(Float, default=2)

    # Concentration (mg/mL) that scales probability of dauering = 250,000 ng total food availablea
    dauer_threshold = Column(Float, default=0.05)
    @classmethod
    def get_schema(cls):
        schema = {}
        for column in cls.__table__.columns:
            column_type = str(column.type)
            if column_type.startswith("VARCHAR") or column_type.startswith("STRING"):
                column_type = "string"
            elif column_type.startswith("FLOAT"):
                column_type = "float"
            elif column_type.startswith("INTEGER"):
                column_type = "integer"
           
            default_value = get_column_default(column)

            schema[column.name] = {
                "type": column_type,
                "default": default_value
            }

        return schema


class WormTimestep(Base):
    __tablename__ = "worms"
    id = Column(Integer, primary_key=True, autoincrement=True)
    Worm_Name = Column(String)
    Timestep = Column(Integer)
    Age_hours = Column(Float)
    Stage = Column(String)
    Mass = Column(Float)
    Egg_Mass = Column(Float)
    Eggs_Laid = Column(Integer)
    Available_Food = Column(Float)
    Total_Appetite = Column(Float)
    Desired_Growth = Column(Float)
    Desired_Eggs = Column(Float)
    Actual_Egg_Investment = Column(Float)
    Metabolic_Cost = Column(Float)
    Amount_Eaten = Column(Float)
    Chance_of_Starvation = Column(Float)
    Chance_of_Dauer_Awakening = Column(Float)
    Chance_of_Death = Column(Float)
    Notes = Column(String)
    Variant = Column(String)
    PosX = Column(Float)
    PosY = Column(Float)
    Colony = Column(Integer)
    

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
    instance: "Simulation" = None
    variants = []

    location_history = []

    clients: Set[WebSocketServerProtocol] = set()

    food_location: List[Food] = [ Food(id, [random.randint(0, WORLD_MAX_X), random.randint(0, WORLD_MAX_Y)], STARTING_FOOD/10) for id in range(10) ]

    @classmethod
    def get_food_location(cls):
        return [
            {
                "pos": food.pos,
                "amount": food.amount,
                "radius": food.radius,
            } for food in cls.food_location
        ]

    @classmethod
    def nearest_food(cls, pos):
        return min([f for f in cls.food_location if f.amount > 0], key=lambda f: distance(f.pos, pos))

    @classmethod
    def update_food_lawn(cls, worms: "Worms"):
        for food in cls.food_location:
            food.my_worms = {}
        
        for worm in worms:
            nearest_food = cls.nearest_food(worm.pos)
            nearest_food.my_worms[worm.name] = worm
        

    def __init__(self, output_location, number_worms=STARTING_WORMS, starting_stage=STARTING_STAGE, starting_food=STARTING_FOOD, length=SIMULATION_LENGTH, report_individuals=False, connection=None):
        self.worms: List[Worm] = Worms()
        self.worms.initialize_worms(number_worms, starting_stage)
        self.dead: List[Worm] = Dead_worms()
 
        self.path = pathlib.Path(output_location)
        self.timestep = 0
        self.time = 0
        self.length = length
        self.report_individuals = report_individuals
        self.connection = connection
        self.bulk_data = []
        self.variants = []
        self.worm_count = [] # list of number of worms in each timestep

        self.food_bulk_data = []
    
    @property
    def food(self):
        return sum(food.amount for food in self.food_location)

    @property
    def food_concentration(self):
        return self.food / 1e6 / FLASK_VOLUME

    
    async def serve(self, socket: WebSocketServerProtocol, path: str):
        print("connected")
        try:
            self.clients.add(socket)
            #await socket.send(json.dumps(self.location_history))
            async for _ in socket:
                pass # clients only listen
        finally:
            self.clients.remove(socket)

    @classmethod
    def load_variants(cls, data: dict, session):
        for d in data["variants"]:
            assert "variant" in d, "Must name the variant"
            G = Genome(**d)
            session.add(G)
            cls.variants.append(G)

        session.commit()


    async def iterate_once(self):
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

        # Halve the feeding schedule
        global WORLD_MAX_X
        global WORLD_MAX_Y

        global FEEDING_SCHEDULE
        if self.time / 24 % DOUBLE_FEED_INTERVAL_DAYS == 0:
            FEEDING_SCHEDULE *= 2
            WORLD_MAX_X = int(WORLD_MAX_X * 1.3)
            WORLD_MAX_Y = int(WORLD_MAX_Y * 1.3)


        # Cull/add bacteria, if applicable:
        if self.time % CULLING_SCHEDULE == 0: self.cull(PERCENT_CULL)

        
        

        # Add food and randomize location
        if self.time % FEEDING_SCHEDULE == 0:
            for food in self.food_location:
                food.pos = [random.randint(0, WORLD_MAX_X), random.randint(0, WORLD_MAX_Y)]
                food.amount += FEEDING_AMOUNT / len(self.food_location)
        self.worms.move()
        
        self.update_food_lawn(self.worms)


        # Calculate appetite
        

        self.worms.compute_appetite() # For simplicity, worms only detect environment once at the start of each time step

        for food_lawn in self.food_location:
            food_lawn.food_history.append(food_lawn.food_concentration)
        # Feed worms, grow worms
        self.worms.eat()
        
        # Metabolic upkeep
        self.worms.tax()
        
        # Run Checks
        self.worms.make_checks() # Using food concentration detected before feeding so you're only starving if you didn't get enough to eat

        self.worm_count.append(len(self.worms))
            
        
        new_location = [w.get_location_history() for w in self.worms]


        self.location_history.append(new_location)

        # Report outcome of timestep
        self.report()

        await self.broadcast_new_location(new_location, Simulation.get_food_location())

        await asyncio.sleep(1/120)


    async def broadcast_new_location(self, new_location, food_location):
        for client in list(self.clients):
            try:
                await client.send(json.dumps({
                    "worms": new_location,
                    "food": food_location
                }))
            except:
                print("Error sending message to client", file=sys.stderr)

    def cull(self, percent):
        """Periodic culling
        Removes a set percentage of the "media," e.g. 10% of all worms and food to simulate prediation.

        TODO stage specific culling
        
        """

        pct_cull = percent / 100

        self.worms.cull(pct_cull)
        for food_lawn in Simulation.food_location:
            food_lawn.amount *= (1 - pct_cull)


    
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
            'p_awaken', 'p_death', 'note', 'variant'
        ]

        if self.report_individuals:

            for w in self.worms:
                w.note = 'Born at timestep {}'.format(self.timestep) if not hasattr(w, 'note') else w.note

                for a in attributes:
                    if not hasattr(w, a):
                        setattr(w, a, None)

                reportlist = [getattr(w, a) for a in attributes]

                wormts = WormTimestep(
                    Worm_Name=w.name,
                    Timestep=self.timestep,
                    Age_hours=w.age,
                    Stage=w.stage,
                    Mass=w.mass,
                    Egg_Mass=w.current_egg_progress,
                    Eggs_Laid=w.eggs_laid,
                    Available_Food=w.sensed_food,
                    Total_Appetite=w.appetite,
                    Desired_Growth=w.growth_mass,
                    Desired_Eggs=w.desired_egg_mass,
                    Actual_Egg_Investment=w.actual_egg_mass,
                    Metabolic_Cost=w.maintenance,
                    Amount_Eaten=w.portion,
                    Chance_of_Starvation=w.p_starve,
                    Chance_of_Dauer_Awakening=w.p_awaken, 
                    Chance_of_Death=w.p_death,
                    Notes=w.note,
                    Variant=w.genome.variant,
                    PosX=w.pos[0],
                    PosY=w.pos[1],
                    Colony=w.colony
                )

                self.bulk_data.append(wormts)

                w.note = ''


            for food in Simulation.food_location:
                foodts = FoodSQL(
                    timestep=self.timestep,
                    x=food.pos[0],
                    y=food.pos[1],
                    amount=food.amount,
                    radius=food.radius,
                    summed_appetite=food.summed_appetite
                )

                self.food_bulk_data.append(foodts)

            if self.timestep % 10 == 0:
                self.connection.bulk_save_objects(self.bulk_data)
                self.connection.bulk_save_objects(self.food_bulk_data)
                self.food_bulk_data = []
                self.bulk_data = []
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
        
        reportlist = [self.timestep, self.time, self.time / 24, self.food_concentration * 1e6 * FLASK_VOLUME, self.food_concentration]

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
        parlad_to_dauer, parlad_to_dauer_mass = ParladToDauerGet()



        if header:
            with open(self.stage_transition, "w") as fp:
                writer = csv.writer(fp, delimiter="\t")
                writer.writerow([
                    "Timestep",
                    "egg_to_larva", "egg_to_larva_mass", 
                    "larva_to_adult", "larva_to_adult_mass", 
                    "larva_to_dauer","larva_to_dauer_mass",
                    "adult_to_bag","adult_to_bag_mass",
                    "dauer_to_larva", "darva_to_larva_mass", 
                    "adult_laid_egg", "adult_laid_egg_mass",
                    "parlad_to_dauer", "parlad_to_dauer_mass",
                
                ])
            
            with open(self.death_transition, "w") as fp:
                
                writer = csv.writer(fp, delimiter="\t")
                fields = ["Timestep"]
                for _class, value in die_ind.items():
                    for cause_of_death in value.keys():
                        for metric in [ "ind", "mass" ]:
                            fields.append(f"{_class}-{cause_of_death}-{metric}")
                writer.writerow(fields)
        
            with open(self.variant_count, "w") as fp:
                writer = csv.writer(fp, delimiter="\t")
                fields = ["Timestep"] + [ variant.variant for variant in Simulation.variants ]
                writer.writerow(fields)

                
        with open(self.stage_transition, "a+") as fp:
            writer = csv.writer(fp, delimiter="\t")
            writer.writerow([self.timestep, 
                    egg_to_larva, egg_to_larva_mass,
                    larva_to_adult, larva_to_adult_mass,
                    larva_to_dauer, larva_to_dauer_mass,
                    adult_to_bag, adult_to_bag_mass,
                    dauer_to_larva, dauer_to_larva_mass,
                    eggs_laid, eggs_laid * EGGMASS,
                    parlad_to_dauer, parlad_to_dauer_mass
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
        

        counter = collections.defaultdict(int)
        for w in self.worms:
            counter[w.genome.variant] += 1
        

        with open(self.variant_count, "a+") as fp:
            writer = csv.writer(fp, delimiter="\t")
            data = [self.timestep] + [ counter[variant.variant] for variant in Simulation.variants ]
            writer.writerow(data)


    async def run(self):
        """Run function
        
        Run simulation for set number of timesteps (three hour increments)
        Standard length of 100 days means running for 800 timesteps

        """
        self.path.mkdir(exist_ok=True)
        self.summary_path = self.path / 'summary.tsv'
        self.death_transition = self.path / 'death_transitions.tsv'
        self.stage_transition = self.path / 'stage_transitions.tsv'
        self.variant_count = self.path / "variant_count.tsv"
    
        if self.report_individuals:
            self.individual_path = self.path / 'indivduals'
            self.individual_path.mkdir(exist_ok=True)
        
        with open(self.path / 'parameters.json', "w") as fp:
            json.dump(param, fp, indent=4)

        self.report(header=True) # Initial conditions/header for output file

        for i in range(1, self.length):
            await self.iterate_once()
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

            pos = [random.randint(0, WORLD_MAX_X), random.randint(0, WORLD_MAX_Y)]
            genome = random.choice(Simulation.variants)

            self.append(stagedict[starting_stage](name, pos, genome, i))
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
    def compute_appetite(self):
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

        for food_lawn in Simulation.food_location:
            appetite = 0
            for name, worm in food_lawn.my_worms.items():
                worm.sensed_food = food_lawn.food_concentration
                worm.get_growth_mass(food_lawn.food_concentration)
                worm.get_egg_mass(food_lawn.food_concentration)
                worm.get_maintenance()
                worm.appetite = (worm.growth_mass + worm.desired_egg_mass + worm.maintenance) / METABOLIC_EFFICIENCY
                appetite += worm.appetite


            food_lawn.summed_appetite = appetite
        


    def eat(self):
        """Worms eat as much as they can based on their growth requirements and appetites of other worms.

        Confused about how portion is handled in the previous model, since portion is calculated and then a second restriction:
        (portion*appetite) / (portion + appetite) appears to be applied. I think this is to keep worms from consuming all the 
        available food. If we know empirically that worms grow at a certain rate in a certain concentration, then I think it makes
        the most sense to assume they eat at least that much bacteria, though.
        
        Returns total amount consumed
        """



        for food_lawn in Simulation.food_location:
            mass = food_lawn.amount
            if mass > food_lawn.summed_appetite or food_lawn.summed_appetite == 0:
                for name, worm in food_lawn.my_worms.items():
                    worm.portion = worm.appetite
                    if distance(worm.pos, food_lawn.pos) < food_lawn.radius:
                        if worm.portion > food_lawn.amount:
                            worm.portion = food_lawn.amount
                        worm.eat(worm.portion)
                        food_lawn.amount -= worm.portion
                    assert food_lawn.amount >= 0, "Food mass cannot be negative"
            else:
                for name, worm in food_lawn.my_worms.items():
                    worm.portion = (worm.appetite / food_lawn.summed_appetite) * mass
                    if distance(worm.pos, food_lawn.pos) < food_lawn.radius:
                        if worm.portion > food_lawn.amount:
                            worm.portion = food_lawn.amount
                        worm.eat(worm.portion)
                        food_lawn.amount -= worm.portion

                    assert food_lawn.amount >= 0, "Food mass cannot be negative"

    def make_checks(self):
        """Runs checks applicable to each worm

        Since this is the only way for new worms to enter the simulation, each check function returns an empty list if there are no new 
        worms, or a list of class objects of the appropriate worm sublcass (Egg or Dauer, for instance). The new arrivals are then appended
        to the Worms object.
        """
        

        for food_lawn in Simulation.food_location:
            
            current_food = food_lawn.food_history[-1]
            prev_food = food_lawn.food_history[-2]

            new_arrivals = [w.make_checks(current_food, prev_food) for _, w in food_lawn.my_worms.items()]
       
            flat_list = [item for sublist in new_arrivals for item in sublist]
            self += [w('worm_' + str(self.total_worm_number + i + 1)) for i, w in enumerate(flat_list)]
            self.total_worm_number += len(flat_list)
    

    def move(self):
        for w in self:
            w.move()



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
    

    genome: Genome # Just a type hint

    def __init__(self, name, start_pos, genome, colony):
        self.name = name
        self._eggs_laid = 0
        self.pos = start_pos
        self.genome = genome
        self.colony = colony
        self.speed = 0

        # self.genome = random.sample([NormalAppetite, FatWorm, SkinnyWorm])



    def cull_maybe(self):
        roll = random.rand()
        if roll <= self.CULL_PERCENT / 100:
            self.die('culled')
    
    def ageup(self):
        self.age += TIMESTEP

    def move(self):
        pass

    def _move(self, x: int, y: int):
        # check within world bound
        if self.pos[0] + x < 0 or self.pos[0] + x >= WORLD_MAX_X:
            x *= -1
        if self.pos[1] + y < 0 or self.pos[1] + y >= WORLD_MAX_Y:
            y *= -1

        self.pos[0] += x
        self.pos[1] += y
    
    def get_location_history(self) -> dict:
        return {
            "name": self.name,
            "pos": self.pos,
            "stage": type(self).__name__,
            "stage_short": type(self).__name__[0].upper(), 
        }

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
    def __init__(self, name, start_pos, genome, colony, *args, **kwargs):
        super().__init__(name, start_pos, genome, colony, *args, **kwargs)
        self.mass = EGGMASS
        self.stage = 'egg'
        self.age = 0
        self.eggs_laid = 0
        self.egg_age = 0 # Eggs hatch after 15 hours, and larvae are born at age 0
       

    def move(self):
        """
        Eggs don't move
        """
        pass 

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
        self.__init__(self.name, self.pos, self.genome, self.colony)

LarvaToDauerSet, LarvaToDauerGet = CreateCounter()
LarvaToAdultSet, LarvaToAdultGet = CreateCounter()

class Larva(Worm):
    """Second stage

    Larvae eat, grow, test their environment to see if they dauer or starve, and potentially become adults after
    a set time and if in a specific mass range.

    """
    def __init__(self, name, pos, genome, colony, *args, **kwargs):
        self.stage = 'larva'
        self.can_dauer = False
        self.eggs_laid = 0
        # The below "if" statements account for situations like if the simulation is being started with larvae,
        # or if a worm is re-entering larvahood after having been a dauer, but still remembers its larval age.
        if not hasattr(self, 'mass'): self.mass = STANDARD_LARVA_MASS
        if not hasattr(self, 'age'): self.age = 0
        if not hasattr(self, 'larval_age'): self.larval_age = 0
        self.p_awaken = None
        super(Larva, self).__init__(name, pos, genome, colony)
    
    def move(self):
        speed = 5
        
        closest_food = Simulation.nearest_food(self.pos)
        
        dy, dx = 0, 0
        # if inside food lawn, don't move towards center
        if distance(self.pos, closest_food.pos) > closest_food.radius:
            dx, dy = move_towards_food(self.pos, closest_food.pos, speed)

        dx += speed * normal()
        dy += speed * normal()

        self._move(dx, dy)



    def ageup(self):
        self.age += TIMESTEP
        self.larval_age += TIMESTEP

    def tax(self):
        self.mass -= self.mass * self.genome.metabolic_tax
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
        self.mass += amount * METABOLIC_EFFICIENCY * self.genome.appetite


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


        dauer_multiplier = 1 #self.dauer_pheremone_function(Simulation.instance.worm_count)

        # Check dauer/starvation:
        self.p_starve = (1 / self.genome.dauer_rate) * math.exp(-0.5 * (current_food + prev_food) / self.genome.dauer_threshold) * dauer_multiplier
        
        

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
        self.__init__(self.name, self.pos, self.genome, self.colony)
    
    @LarvaToAdultSet
    def molt(self):
        """Mature to adult
        """
        self.__class__ = Adult
        self.__init__(self.name, self.pos, self.genome, self.colony)


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
    def __init__(self, name, start_pos, genome, colony):
        self.stage = 'dauer'
        self.has_dauered = True
        self.eggs_laid = 0
        if not hasattr(self, 'mass'): self.mass = STANDARD_LARVA_MASS
        if not hasattr(self, 'age'): self.age = 15 # assuming 30 hours from parlad bagginng -> 15 hours for eggs to hatch, 15 hours for dauers to develop
        self.p_starve = None
        
        super(Dauer, self).__init__(name, start_pos, genome, colony)

    
    def move(self):
        """
        
        Dauers don't move"""
        pass
        
  

    def make_checks(self, current_food, prev_food):
        """Dauers check if conditions are safe to exit dauer. Currently I'm treating dauers as immortal

        TODO: Add dauer attrition
        """

        self.p_awaken = DAUER_EXIT_PROB * math.sqrt(0.5 * (current_food + prev_food) * 1e6 * FLASK_VOLUME) # Converted back to ng here for convenience
                                                                                                           # Could also just use converted dauer exit probability (3.24e-5 * sqrt(5e6) = 0.0724486)
        roll = random.rand()
        if roll < self.p_awaken:
            self.exit_dauer()

        return []
    
    @DauerToLarvaSet
    def exit_dauer(self):
        self.__class__ = Larva
        self.__init__(self.name, self.pos, self.genome, self.colony)

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

    def __init__(self, name, pos, genome, colony):
        self.stage = 'adult'
        self.adult_age = 0
        self.total_egg_mass = 0
        self.min_somatic_mass = (self.mass + MIN_ADULT_MASS) / 2
        self.note = 'Min adult mass set to {}'.format(self.min_somatic_mass)
        self.bag_rate = BAG_RATE
        self.bag_threshold = BAG_THRESHOLD
        super(Adult, self).__init__(name, pos, genome, colony)

    def ageup(self):
        self.age += TIMESTEP
        self.adult_age += TIMESTEP

    def tax(self):
        self.mass -= self.mass * self.genome.metabolic_tax

    def get_maintenance(self):
        self.maintenance = self.mass * self.genome.metabolic_tax
    

    def move(self):
        speed = 10

        dy, dx = 0, 0
        closest_food = Simulation.nearest_food(self.pos)

        # if inside food lawn, don't move towards center
        if distance(self.pos, closest_food.pos) > closest_food.radius:
            dx, dy = move_towards_food(self.pos, closest_food.pos, speed)

        dx += speed * normal()
        dy += speed * normal()

        self._move(dx, dy)


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
        assert self.mass > 0, "mass > actual_egg_mass"
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

        self.p_death = self.genome.life_span * (math.exp(self.age / gompertzTau) - 1) / gompertzA # Probability of dying at given time

        
        roll = random.rand()
        if roll < self.p_death:
            self.die('old_age')

        return eggs

    def lay_egg(self, number):
        self.eggs_laid += number

        Egg_partial = functools.partial(Egg, genome=self.genome, start_pos=list(self.pos), colony=self.colony)

        return [Egg_partial] * int(number)

    @AdultToBagSet
    def bag(self):
        self.__class__ = Parlad
        self.__init__(self.name, self.pos, self.genome, self.colony)


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
    def __init__(self, name, pos, genome, colony):
        self.stage = 'parlad'
        self.lifespan = self.age
        self.cause_of_death = ('bag')
        self.mass += self.total_egg_mass - (self.eggs_laid * EGGMASS) # Total mass (unlaid eggs + somatic mass)
        self.dauer_potential = int((self.mass * BAG_EFFICIENCY) // STANDARD_LARVA_MASS)
        self.mass_decrement = self.mass / (30 / TIMESTEP) # Will lose this much mass per timestep (converted into dauers) down to 0
        self.note = 'Will burst into {} dauers in 30 hours'.format(self.dauer_potential)
        super(Parlad, self).__init__(name, pos, genome, colony)
        assert self.mass > 0, f"{self.total_egg_mass - (self.eggs_laid * EGGMASS)}, {self.genome}"

    def move(self):
        # Do Parlads move?
        speed = 0
        dx, dy = normal() * speed, normal() * speed

        self._move(dx, dy)


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

            dauer_init = functools.partial(Dauer, genome=self.genome, start_pos=list(self.pos), colony=self.colony)

            released_dauers.extend([dauer_init]*self.dauer_potential)
            self.mass = self.mass - self.dauer_potential * STANDARD_LARVA_MASS
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
        super(Dead, self).__init__(name, self.pos, self.genome, self.colony)

#%%

# TODO - Change Egg Mass efficiency
# TODO - Change Life Span
# TODO - Parameterize the genome


Egg.CULL_PERCENT = EGG_CULL_PERCENT
Larva.CULL_PERCENT = LARVA_CULL_PERCENT
Dauer.CULL_PERCENT = DAUER_CULL_PERCENT
Adult.CULL_PERCENT = ADULT_CULL_PERCENT
Parlad.CULL_PERCENT = PARLAD_CULL_PERCENT


def worm_world_animation_data(database_file):
    import collections
    # Create connection to database
    engine = create_engine(f"sqlite:///{database_file}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)

    
    
    with Session() as session:
        # read worm data
        worm_data = collections.defaultdict(list)
        worms: Iterable[WormTimestep] = session.query(WormTimestep).all()
        for w in worms:
            worm_data[w.Timestep].append({
                "timestep": w.Timestep,
                "stage": w.Stage,
                "name": w.Worm_Name,
                "x": w.PosX,
                "y": w.PosY,
                "colony": w.Colony,
                "variant": w.Variant,
            })
        
        food_data = collections.defaultdict(list)
        foods: Iterable[FoodSQL] = session.query(FoodSQL).all()
        for f in foods:
            food_data[f.timestep].append({
                "food_id": f.food_id,
                "amount": f.amount,
                "radius": f.radius,
                "x": f.x,
                "y": f.y,
            })
        
        return {
            "worm": worm_data,
            "food": food_data,
        }

            
    


async def main():
    if args["--variants"]:
        variants_file = args["--variants"]
        with open(variants_file) as fp:
            variants_data = json.load(fp)
    

    if args["--database"]:
        # Establish connection
        engine =  create_engine(f"sqlite:///{directory}/{args['--database']}")

    else:
        engine = create_engine(f"sqlite://")

    # Create all the table
    Base.metadata.create_all(engine)

    # Open the session

    

    Session = sessionmaker(bind=engine)
    with Session() as session:
        Simulation.load_variants(variants_data, session)
        Simulation.instance = test = Simulation(directory, connection=session, report_individuals=report_individuals)
    
        if args["--socket"]:
            port = int(args["--socket"])
            print("Server started on port", port)
            async with websockets.serve(test.serve, "localhost", port) as server:
                await test.run()
        else:
            await test.run()

if __name__ == "__main__":
    import asyncio

    loop = asyncio.get_event_loop()
    main_task = loop.create_task(main())

    def cleanup(task: asyncio.Task):
        # retrieve the exception
        
        if task.exception():
            stack = task.get_stack()
            for frame in stack:
                filename = frame.f_code.co_filename
                lineno = frame.f_lineno
                print(f"Error: {task.exception()} at {filename}:{lineno}")
            exit()

    main_task.add_done_callback(cleanup)

    loop.run_forever()
