# genome.py
from sqlalchemy import Column, String, Float
from sqlalchemy.ext.declarative import declarative_base
from utils import get_column_default
from database import Base
from simulation_globals import *
from sqlalchemy.orm import Session
#Base = declarative_base()

class Genome(Base):
    __tablename__ = "Genome"
    variant = Column(String, primary_key=True)

    appetite: float = Column(Float, default=1)
    life_span: float = Column(Float, default=1)
    metabolic_tax: float = Column(Float, default=0.035)

    eggN: float = Column(Float, default=1)
    eggM: float = Column(Float, default=1)
    eggScale: float = Column(Float, default=1)
    
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
           
            default_value = get_column_default(column)

            schema[column.name] = {
                "type": column_type,
                "default": default_value
            }

        return schema

def load_variants(data: dict, session: Session):
    """
    Loads variants from a given dictionary and inserts them into the database if they don't already exist.

    Args:
        data (dict): Dictionary containing the variant data.
        session (Session): SQLAlchemy session used to interact with the database.

    Raises:
        AssertionError: If a variant does not have a "variant" key.
    """
    for d in data["variants"]:
        assert "variant" in d, "Must name the variant"
        
        # Check if the variant already exists
        existing_variant = session.query(Genome).filter_by(variant=d["variant"]).first()
        
        if existing_variant:
            print(f"Variant {d['variant']} already exists: {existing_variant}")
        else:
            # Add the new record
            G = Genome(**d)
            session.add(G)
            variants.append(G)
    
    session.commit()  # Commit changes to the database