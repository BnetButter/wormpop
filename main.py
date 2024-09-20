#!/usr/bin/env python3

"""
USAGE: simulation [--parameters=<string>] [ --database=<string> ] [ --name=<string> ] [--directory=<string>] [ --variants=<string> ] [ --socket=<port> ] [ --report-individuals ]
"""

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import os, sys
import json
import docopt
from pathlib import Path
from config import set_config_path
import simulation_globals

args = docopt.docopt(__doc__)
#set path of config file so it can be loaded in other modules
parameters = args["--parameters"]
set_config_path(parameters)

from simulation import Simulation
from database import Base
from genome import load_variants

def main():
    args = docopt.docopt(__doc__)
    parameters = args["--parameters"]
    set_config_path(parameters)

    database = args["--database"] or ":memory:"
    directory = args["--directory"] or "Simulation"

    if args["--variants"]:
        variants_file = args["--variants"]
        with open(variants_file) as fp:
            variants_data = json.load(fp)

    # Establish connection
    directory = args["--directory"] if args["--directory"] else "Simulation"
    database_file = Path(directory) / args['--database']
    if database_file.is_file():
        print(f"Deleting {database_file.stem}")
        database_file.unlink()
    engine = create_engine(f"sqlite:///{directory}/{args['--database']}")

    Base.metadata.create_all(engine)

    Session = sessionmaker(bind=engine)
    # Create the session and run the simulation
    with Session() as session:
        #loaded in as a list of genomes of each variants
        #print(simulation_globals.variants)
        load_variants(variants_data, session)
        #print(simulation_globals.variants)
        #sys.exit()
        simulation = Simulation(
            directory, 
            connection=session, 
            report_individuals=args["--report-individuals"], 
            engine=engine
        )    
        try:
            simulation.run()
        finally:
            simulation.on_simulation_end()

if __name__ == "__main__":
    main()
