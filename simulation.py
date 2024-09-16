# simulation.py
import pathlib
import numpy as np
from worms import Worm, Dead
from genome import Genome
from database import SimulationSummary
from reporting import CreateCounter, CreateDeathCounter

from config import load_constants
constants = load_constants()

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
    dynamic_table = None

    def __init__(self, output_location, number_worms=STARTING_WORMS, starting_stage=STARTING_STAGE, starting_food=STARTING_FOOD, length=SIMULATION_LENGTH, report_individuals=False, connection=None, engine=None):
        self.worms = Worm()
        self.worms.initialize_worms(number_worms, starting_stage)
        self.dead = Dead()
        self.food = starting_food
        self.food_concentration = self.food / 1e6 / FLASK_VOLUME  # convert to mg / mL
        self.food_history = [self.food_concentration] # Used to keep track of how much food each worm has seen
        self.path = pathlib.Path(output_location)
        self.timestep = 0
        self.time = 0
        self.length = length
        self.report_individuals = report_individuals
        self.connection = connection
        self.bulk_data = []
        self.bulk_death_data = [] # death data handled differently
        self.variants = []
        self.worm_count = [] # list of number of worms in each timestep
        self.engine = engine
        self._summary = SimulationSummary()
        self.worms_that_laid_eggs = set()

    @classmethod
    def load_variants(cls, data: dict, session):
        for d in data["variants"]:
            assert "variant" in d, "Must name the variant"
            
            # Check if the variant already exists
            existing_variant = session.query(Genome).filter_by(variant=d["variant"]).first()
            
            # Add the new record
            G = Genome(**d)
            session.add(G)
            cls.variants.append(G)
        
        print(data)
        session.commit()


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
        if self.time % CULLING_SCHEDULE == 0: 
            self.cull(PERCENT_CULL)

        if self.time % FEEDING_SCHEDULE == 0:
            self._summary.Bacteria_In_mg += FEEDING_AMOUNT
            self.food += FEEDING_AMOUNT

        # Calculate appetite
        self.food_concentration = self.food / 1e6 / FLASK_VOLUME # Convert from nanograms to mg/mL
        self.food_history.append(self.food_concentration)
        self.worms.compute_appetite(self.food_concentration) # For simplicity, worms only detect environment once at the start of each time step

        # Feed worms, grow worms
        amount_consumed = self.worms.eat(self.food)
        self.food -= amount_consumed

        # Metabolic upkeep
        self.worms.tax()
        
        # Run Checks
        self.worms.make_checks(self.food_history) # Using food concentration detected before feeding so you're only starving if you didn't get enough to eat

        self.worm_count.append(len(self.worms))

        # Report outcome of timestep
        self.report()


    def cull(self, percent):
        """Periodic culling
        Removes a set percentage of the "media," e.g. 10% of all worms and food to simulate prediation.

        TODO stage specific culling
        
        """

        pct_cull = percent / 100

        self.worms.cull(pct_cull)
        culled = self.food * pct_cull
        Simulation.instance._summary.Bacteria_Culled_mg += culled
        Simulation.instance._summary.Bacteria_Culled_Percent = Simulation.instance._summary.Bacteria_Culled_mg / Simulation.instance._summary.Bacteria_In_mg
        self.food -= culled



    
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
            w:Worm
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
                    Metabolic_Efficiency_Loss=w._metabolic_efficiency_loss,
                    Age_days = w.age / 24,
                    Actual_Growth = w._actual_growth,
                    Able_To_Dauer = w.can_dauer,
                    Eggs_Laid_Timestep = w._eggs_laid_timestep

                )

                self.bulk_data.append(wormts)

                w.note = ''

            if self.timestep % 10 == 0:
                self.connection.bulk_save_objects(self.bulk_data)
                self.bulk_data = []
                self.connection.commit()

        self.dead.extend([w for w in self.worms if w.stage == 'dead'])
        self.worms[:] = [w for w in self.worms if w.stage != 'dead']

        # Group reporting:
        
        if header:
            with open(self.summary_path, 'w+') as file:
                file.write('\t'.join(['Timestep','Time (hours)','Time (days)', 'Food Mass (ng)', 'Food Conc (mg/mL)', 'Number Worms', 'Number Eggs','Number Larvae', 'Number Dauer',
                'Number Adults', 'Number Parlads',"Number L1 Arrest", 'Number Dead','Total Worm Mass (ng)','Egg Mass','Larva Mass','Dauer Mass','Adult Mass','Parlad Mass','Dead Mass','Eggs Laid',
                'Died of old age', 'Died of starvation','Died of bagging','Died of predation','Died of arrested development', ])+'\n')

        stages = ['egg','larva','dauer','adult','parlad', 'L1_Arrest']
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

        larva_to_l1, larva_to_l1_mass = LarvaToL1ArrestGet()
        l1_to_larva, l1_to_larva_mass = L1ArrestToLarvaGet()


        transition = StageTransition(
            timestep=self.timestep,
            egg_to_larva=egg_to_larva, egg_to_larva_mass=egg_to_larva_mass,
            larva_to_adult=larva_to_adult, larva_to_adult_mass=larva_to_adult_mass,
            larva_to_dauer=larva_to_dauer, larva_to_dauer_mass=larva_to_dauer_mass,
            adult_to_bag=adult_to_bag, adult_to_bag_mass=adult_to_bag_mass,
            dauer_to_larva=dauer_to_larva, dauer_to_larva_mass=dauer_to_larva_mass,
            adult_laid_egg=eggs_laid, adult_laid_egg_mass=eggs_laid*EGGMASS,
            parlad_to_dauer=parlad_to_dauer, parlad_to_dauer_mass=parlad_to_dauer_mass,
            larva_to_l1arrest=larva_to_l1, larva_to_l1arrest_mass=larva_to_l1_mass,
            l1arrest_to_larva=l1_to_larva, l1arrest_to_larva_mass=l1_to_larva_mass
        )

        self.bulk_data.append(transition)


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
                    "larva_to_l1arrest", "larva_to_l1arrest_mass",
                    "l1arrest_to_larva", "l1arrest_to_larva_mass",
                
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
            
            metadata = MetaData(bind=self.engine)

            self.dynamic_table = create_dynamic_table(metadata, 'dynamic_stage_transition', die_ind)
            metadata.create_all()

        
        
        death_data = create_entry_data(self.timestep, die_ind, death_metrics)

        self.bulk_death_data.append(death_data)

        if self.timestep % 10 == 0 and self.dynamic_table is not None:
            self.connection.execute(
                self.dynamic_table.insert(),
                self.bulk_death_data
            )
            self.bulk_death_data = []

       



                
        with open(self.stage_transition, "a+") as fp:
            writer = csv.writer(fp, delimiter="\t")
            writer.writerow([self.timestep, 
                    egg_to_larva, egg_to_larva_mass,
                    larva_to_adult, larva_to_adult_mass,
                    larva_to_dauer, larva_to_dauer_mass,
                    adult_to_bag, adult_to_bag_mass,
                    dauer_to_larva, dauer_to_larva_mass,
                    eggs_laid, eggs_laid * EGGMASS,
                    parlad_to_dauer, parlad_to_dauer_mass,
                    larva_to_l1, larva_to_l1_mass,
                    l1_to_larva, l1_to_larva_mass
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
    

    def on_simulation_end(self):
        self._summary.Num_Timestep = self.timestep
        self._summary.Bacteria_remaining_mg = self.food
        alive = self._summary.Worms_alive_at_last_timestep = len(self.worms)
        dead = self._summary.Worms_dead_at_last_timestep = len(self.dead)

        if dead:
            self._summary.Worms_died_cull_percent = self._summary.Worms_died_cull / dead
            self._summary.Worms_died_bag_percent = self._summary.Worms_died_bag / dead
            self._summary.Worms_died_starvation_percent = self._summary.Worms_died_starvation / dead
            self._summary.Worms_died_old_age_percent = self._summary.Worms_died_old_age / dead
            self._summary.Worms_died_arrested_development_percent = self._summary.Worms_died_arrested_development / dead
        
        total_born = self._summary.Worms_born_egg + self._summary.Worms_born_dauer
        if total_born:
            self._summary.Worms_born_egg_percent = self._summary.Worms_born_egg / total_born
            self._summary.Worms_born_dauer_percent = self._summary.Worms_born_dauer / total_born
        
        self._summary.Worms_Laid_Eggs_no = len(self.worms_that_laid_eggs)
        if alive + dead:
            self._summary.Worms_Laid_Eggs_percent = self._summary.Worms_Laid_Eggs_no / (alive + dead)

        self.connection.add(self._summary)
        self.connection.commit()

        self.connection.bulk_save_objects(self.bulk_data)
        self.connection.commit()


    def run(self):
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
            self.iterate_once()
            if self.timestep % 10 == 0:
                print('{} Timesteps, Food = {} mg/mL, {} Worms Alive, {} Worms Dead'.format(self.timestep, round(self.food_concentration,2), len(self.worms), len(self.dead)))
