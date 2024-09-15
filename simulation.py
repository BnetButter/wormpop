from sqlalchemy import create_engine, sessionmaker
from models import WormTimestep, WormSummary, Genome, SimulationSummary, StageTransition, Worms, Dead_worms, Worm, Egg, Larva, Adult, Dauer, Parlad
import os
import json
import docopt
from constants import SIMULATION_LENGTH, TIMESTEP, STARTING_WORMS, STARTING_STAGE, STARTING_FOOD, FLASK_VOLUME, FEEDING_SCHEDULE, CULLING_SCHEDULE, PERCENT_CULL

def main():
    args = docopt.docopt(__doc__)
    parameters = args["--parameters"]
    database = args["--database"] or ":memory:"
    directory = args["--directory"] or "Simulation"

    # Establish connection
    engine = create_engine(f"sqlite:///{directory}/{args['--database']}")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)

    # Create the session and run the simulation
    with Session() as session:
        Simulation.load_variants(variants_data, session)
        sim = Simulation(directory, session, args["--report-individuals"], engine)
        try:
            sim.run()
        finally:
            sim.on_simulation_end()

if __name__ == "__main__":
    main()
