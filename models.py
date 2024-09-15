from sqlalchemy import (
    Column, Integer, String, Float, DateTime, create_engine, ForeignKey
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
import datetime
import os

Base = declarative_base()

# SQLAlchemy Models

class Genome(Base):
    __tablename__ = "genome"
    variant = Column(String, primary_key=True)
    
    appetite = Column(Float, default=1)
    life_span = Column(Float, default=1)
    metabolic_tax = Column(Float, default=0.035)

    eggN = Column(Float, default=1)
    eggM = Column(Float, default=1)
    eggScale = Column(Float, default=1)
    
    dauer_probability = Column(Float, default=0)

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
           
            default_value = column.default.arg if column.default else None

            schema[column.name] = {
                "type": column_type,
                "default": default_value
            }

        return schema

class WormTimestep(Base):
    __tablename__ = "worm_timestep"
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
    Metabolic_Efficiency_Loss = Column(Float)
    Desired_Growth = Column(Float)
    Actual_Growth = Column(Float)
    Desired_Eggs = Column(Float)
    Actual_Egg_Investment = Column(Float)
    Egg_Mass = Column(Float)
    Eggs_Laid = Column(Integer)
    Eggs_Laid_Timestep = Column(Integer)
    Metabolic_Cost = Column(Float)
    Chance_of_Starvation = Column(Float)
    Chance_of_Dauer_Awakening = Column(Float)
    Chance_of_Death = Column(Float)
    Able_To_Dauer = Column(Integer)
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
    Timestamp = Column(DateTime, default=datetime.datetime.now)
    Commit = Column(String)
    RepoURL = Column(String)
    Branch = Column(String)
    Experimentor = Column(String, default=os.environ.get("USER", 'default_user'))
    Num_Timestep = Column(Integer)
    
    # Bacteria
    Bacteria_In_mg = Column(Float)
    Bacteria_Culled_mg = Column(Float)
    Bacteria_Culled_Percent = Column(Float)
    Bacteria_to_worm_ingested_mg = Column(Float)
    Bacteria_to_worm_ingested_percent = Column(Float)
    Bacteria_to_worm_somatic_mass_mg = Column(Float)
    Bacteria_to_worm_somatic_mass_percent = Column(Float)
    Bacteria_to_worm_eggs_mg = Column(Float)
    Bacteria_to_worm_eggs_percent = Column(Float)
    Bacteria_to_worm_metabolic_tax_mg = Column(Float)
    Bacteria_to_worm_metabolic_tax_percent = Column(Float)
    Bacteria_to_worm_metabolic_inefficiency_mg = Column(Float)
    Bacteria_to_worm_metabolic_inefficiency_percent = Column(Float)
    Bacteria_remaining_mg = Column(Float)
    
    # Worms and their life stages
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
    Worms_Laid_Eggs_percent = Column(Integer)
    
    # Population metrics
    Timestep_with_highest_population = Column(Integer)
    Maximum_population_number = Column(Integer)
    Lowest_population_number_after_max = Column(Integer)
    Lowest_population_number_after_max_timestep = Column(Integer)
    Average_number_of_worms = Column(Float)
    Average_number_of_adults = Column(Float)
    Average_number_of_parlads = Column(Float)
    Average_number_of_larva = Column(Float)
    Average_number_of_eggs = Column(Float)
    Average_number_of_dauer = Column(Float)


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

