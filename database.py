from sqlalchemy import Column, Integer, String, Float, DateTime
from sqlalchemy.ext.declarative import declarative_base
import os

Base = declarative_base()

class WormTimestep(Base):
    __tablename__ = "worms"
    id = Column(Integer, primary_key=True, autoincrement=True)
    Worm_Name = Column(String)
    Timestep = Column(Integer)
    Available_Food = Column(Float)
    Age_hours = Column(Float)
    Age_days = Column(Float)
    Stage = Column(String)
    Mass = Column(Float)
    Total_Appetite = Column(Float)
    Amount_Eaten = Column(Float)
    Metabolic_Efficiency_Loss = Column(Float) # NEW FIELD
    Desired_Growth = Column(Float)
    Actual_Growth = Column(Float) # NEW FIELD
    Desired_Eggs = Column(Float)
    Actual_Egg_Investment = Column(Float)
    Egg_Mass = Column(Float)
    Eggs_Laid = Column(Integer)
    Eggs_Laid_Timestep = Column(Integer) # NEW FIELD
    Metabolic_Cost = Column(Float)
    Chance_of_Starvation = Column(Float)
    Chance_of_Dauer_Awakening = Column(Float)
    Chance_of_Death = Column(Float)
    Able_To_Dauer = Column(Integer) # NEW FIELD
    Notes = Column(String)
    Variant = Column(String)

class WormSummary(Base):
    __tablename__ = "worm_summary"
    Worm_Name = Column(String, primary_key=True)
    Larva_span_days = Column(Float)

    L1_span_days = Column(Float)
    LarvaPost_L1_span_days = Column(Float)

    Dauer_span_days = Column(Float)
    LarvaPostDauer_span_days = Column(Float)
    Adult_span_days = Column(Float)
    Parlad_span_days = Column(Float)
    Life_span_days = Column(Float)

    Total_Food_Consumed = Column(Float)
    Total_Eggs_Laid = Column(Integer)
    Total_Mass = Column(Float)
    Total_Body_Mass = Column(Float)
    Total_Egg_Mass = Column(Float, default=0)

    Total_Metabolic_Cost = Column(Float)
    Total_Metabolic_Tax = Column(Float)
    Parlad_Mass_Converted = Column(Float)
    Reproductive_Span = Column(Float)
    Cause_of_Death = Column(String)

    def __init__(self, Worm_Name):
        super().__init__(Worm_Name=Worm_Name)
        self.Total_Food_Consumed = 0
        self.Total_Eggs_Laid = 0
        self.Total_Egg_Mass = 0
        self.Total_Body_Mass = 0
        self.Total_Metabolic_Cost = 0
        self.Total_Metabolic_Tax = 0
        self.Parlad_Mass_Converted = 0


class SimulationSummary(Base):
    __tablename__ = "simulation_summary"

    id = Column(Integer, primary_key=True, autoincrement=True)
    Timestamp = Column(DateTime)
    Commit = Column(String)
    RepoURL = Column(String)
    Branch = Column(String)
    Experimentor = Column(String, default=os.environ.get("USER", 'default_user'))

    Num_Timestep = Column(Integer)

    # 1. Bacteria all
    Bacteria_In_mg = Column(Float)
    Bacteria_Culled_mg = Column(Float) 
    Bacteria_Culled_Percent = Column(Float)

    Bacteria_to_worm_ingested_mg = Column(Float)
    Bacteria_to_worm_ingested_percent = Column(Float) # TODO
    Bacteria_to_worm_somatic_mass_mg = Column(Float)
    Bacteria_to_worm_somatic_mass_percent = Column(Float) # TODO
    
    Bacteria_to_worm_eggs_mg = Column(Float)
    Bacteria_to_worm_eggs_percent = Column(Float) # TODO
    Bacteria_to_worm_metabolic_tax_mg = Column(Float)
    Bacteria_to_worm_metabolic_tax_percent = Column(Float) #TODO

    Bacteria_to_worm_metabolic_inefficiency_mg = Column(Float)
    Bacteria_to_worm_metabolic_inefficiency_percent = Column(Float) # TODO
    
    Bacteria_remaining_mg = Column(Float)
    Bacteria_average = Column(Float) # TODO
    Bacteria_max = Column(Float) # TODO
    Bacteria_min = Column(Float) # TODO

    Worms_born_dauer = Column(Integer)
    Worms_born_dauer_percent = Column(Float)
    Worms_born_egg = Column(Integer)
    Worms_born_egg_percent = Column(Float)

    
    Worms_died_cull = Column(Integer)
    Worms_died_cull_percent = Column(Float)
    Worms_died_starvation = Column(Integer)
    Worms_died_starvation_percent = Column(Float)
    Worms_died_bag = Column(Integer)
    Worms_died_bag_percent = Column(Integer)
    Worms_died_old_age = Column(Integer)
    Worms_died_old_age_percent = Column(Float)
    Worms_died_arrested_development = Column(Integer)
    Worms_died_arrested_development_percent = Column(Float)


    Worms_alive_at_last_timestep = Column(Integer)
    Worms_dead_at_last_timestep = Column(Integer)


    Worms_Laid_Eggs_no = Column(Integer)
    Worms_Laid_Eggs_percent = Column(Integer) # TODO
    Worms_average_laid_eggs = Column(Integer) # TODO
    Worms_average_repro_span = Column(Float) # TODO
    Worms_average_repro_span_exclude_culled = Column(Float) #TODO
    Worms_average_repro_span_excude_culled_starve = Column(Float) #TODO

    Average_dauer_no = Column(Float) # TODO
    Average_Dauer_percent = Column(Float) # TODO
    Average_Dauer_from_Parlads = Column(Float) # TODO
    Average_Dauer_from_Parlads_percent = Column(Float) #TODO

    Average_Dauer_from_Larvae = Column(Float) #TODO
    Average_Dauer_from_Larvae_percent = Column(Float) #TODO

    Timestep_with_highest_population = Column(Integer) #TODO
    Maximum_population_number = Column(Integer) # TODO
    Lowest_population_number_after_max = Column(Integer) #TODO
    Lowest_population_number_after_max_timestep = Column(Integer) #TODO

    Average_number_of_worms = Column(Float) # TODO
    Average_number_of_adults = Column(Float) # TODO
    Average_number_of_parlads = Column(Float) # TODO
    Average_number_of_larva = Column(Float) # TODO
    Average_number_of_eggs = Column(Float) # TODO
    Average_number_of_dauer = Column(Float) # TODO

    def __init__(self):
        self.Timestamp = datetime.datetime.now()
        self.Commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode('utf-8')
        self.RepoURL = subprocess.check_output(["git", "config", "--get", "remote.origin.url"]).strip().decode('utf-8')
        self.Branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).strip().decode('utf-8')
        
        self.Num_Timestep = Column(Integer)

        self.Bacteria_In_mg = STARTING_FOOD
        self.Bacteria_Culled_mg = 0
        self.Bacteria_Culled_Percent = 0

        self.Bacteria_to_worm_ingested_mg = 0
        self.Bacteria_to_worm_ingested_percent = 0
        self.Bacteria_to_worm_somatic_mass_mg = 0
        self.Bacteria_to_worm_somatic_mass_percent = 0
        
        self.Bacteria_to_worm_eggs_mg = 0
        self.Bacteria_to_worm_eggs_percent = 0
        self.Bacteria_to_worm_metabolic_tax_mg = 0
        self.Bacteria_to_worm_metabolic_tax_percent = 0

        self.Bacteria_to_worm_metabolic_inefficiency_mg = 0
        self.Bacteria_to_worm_metabolic_inefficiency_percent = 0
        
        self.Bacteria_remining_mg = 0
        self.Bacteria_average = 0
        self.Bacteria_max = 0
        self.Bacteria_min = 0

        self.Worms_born_dauer = 0
        self.Worms_born_dauer_percent = 0
        self.Worms_born_egg = 0
        self.Worms_born_egg_percent = 0
        
        self.Worms_died_cull = 0
        self.Worms_died_cull_percent = 0
        self.Worms_died_starvation = 0
        self.Worms_died_starvation_percent = 0
        self.Worms_died_bag = 0
        self.Worms_died_bag_percent = 0
        self.Worms_died_old_age = 0
        self.Worms_died_old_age_percent = 0
        self.Worms_died_arrested_development = 0
        self.Worms_died_arrested_development_percent = 0
        self.Worms_alive_at_last_timestep = 0
        self.Worms_dead_at_last_timestep = 0

        self.Worms_Laid_Eggs_no = 0
        self.Worms_Laid_Eggs_percent = 0
        self.Worms_average_laid_eggs = 0
        self.Worms_average_repro_span = 0
        self.Worms_average_repro_span_exclude_culled = 0
        self.Worms_average_repro_span_excude_culled_starve = 0

        self.Average_dauer_no = 0
        self.Average_Dauer_percent = 0
        self.Average_Dauer_from_Parlads = 0
        self.Average_Dauer_from_Parlads_percent = 0

        self.Average_Dauer_from_Larvae = 0
        self.Average_Dauer_from_Larvae_percent = 0

        self.Timestep_with_highest_population = 0
        self.Maximum_population_number = 0
        self.Lowest_population_number_after_max = 0
        self.Lowest_population_number_after_max_timestep = 0

        self.Average_number_of_worms = 0
        self.Average_number_of_adults = 0
        self.Average_number_of_parlads = 0
        self.Average_number_of_larva = 0
        self.Average_number_of_eggs = 0
        self.Average_number_of_dauer = 0


class StageTransition(Base):
    __tablename__ = 'stage_transition'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    timestep = Column(Integer, nullable=False)
    egg_to_larva = Column(Integer, nullable=False)
    egg_to_larva_mass = Column(Float, nullable=False)
    larva_to_adult = Column(Integer, nullable=False)
    larva_to_adult_mass = Column(Float, nullable=False)
    larva_to_dauer = Column(Integer, nullable=False)
    larva_to_dauer_mass = Column(Float, nullable=False)
    adult_to_bag = Column(Integer, nullable=False)
    adult_to_bag_mass = Column(Float, nullable=False)
    dauer_to_larva = Column(Integer, nullable=False)
    dauer_to_larva_mass = Column(Float, nullable=False)
    adult_laid_egg = Column(Integer, nullable=False)
    adult_laid_egg_mass = Column(Float, nullable=False)
    parlad_to_dauer = Column(Integer, nullable=False)
    parlad_to_dauer_mass = Column(Float, nullable=False)
    larva_to_l1arrest = Column(Integer, nullable=False)
    larva_to_l1arrest_mass = Column(Float, nullable=False)
    l1arrest_to_larva = Column(Integer, nullable=False)
    l1arrest_to_larva_mass = Column(Float, nullable=False)

def create_dynamic_table(metadata, table_name, die_ind):
    columns = [
        Column('id', Integer, primary_key=True, autoincrement=True),
        Column('timestep', Integer, nullable=False)
    ]
    
    for _class, causes in die_ind.items():
        for cause_of_death in causes.keys():
            columns.append(Column(f"{_class}_{cause_of_death}_ind", Integer, nullable=False))
            columns.append(Column(f"{_class}_{cause_of_death}_mass", Float, nullable=False))
    
    dynamic_table = Table(table_name, metadata, *columns)
    return dynamic_table