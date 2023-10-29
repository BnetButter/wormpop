Windows:

    Press the Windows key.

    Type "cmd" or "PowerShell" and press Enter to open the Command Prompt or PowerShell.

    Use the ‘cd’ command to navigate to the directory where your wormpop.py script is located. For example, if your script is in C:\Example\Directory\Path, you can use the following command in Command Prompt:

            cd C:\Example\Directory\Path

    You can use a text editor like Notepad or a code editor like Visual Studio Code to open wormpop.py. For example, to open it with Notepad, you can run the following command in Command Prompt:

            notepad wormpop.py  

            (Remember: you cannot run a simulation on notepad)

    To run the simulation on a python file, you need to run the following command in Command prompt:

            python3 wormpop.py

    Check and uncomment the code:

        Use your text editor to locate the last line of code and check if the lines "test = Simulation ('speed test')" and "test.run()" are commented out (prefixed with a #). If they are, simply delete the # to uncomment them.

        You can modify constants, all UPPERCASE variable below “# constants”,  in the code as needed but remember to note the original values and change them back when you're done running your simulations.

        Save your changes by pressing “Ctrl + S”

    End the simulation:

        Press “Ctrl + C” to stop the simulation

----------------------------------------------------------------------------------------------------------
Linux:

    Open the Linux terminal

    Navigate to the directory that contains wormpop.py

        cd Example/Directory/Path/)

    If you wish to go back a directory

        cd ..

    Open wormpop.py

        nano wormpop.py

    Go to the last lines of code and check to see if "test = Simulation('speedtest')" and "test.run()" are commented out. (A "#" before them signifies that they are commented out.)

        If so, delete the "#" to uncomment them so that the simulation can run.

    Modify and constants at the top of the code to test certain changes to constants in the simulation environment.

        The constants are at the top of the code and are usually in all capital letters.

        Make sure to note the original values of the constants before changing them so that they can be changed back later.

    Exit the code and save any changes.

        Ctrl + X | When asked "Save modified buffer" hit Y for yes or N for no

    Run the simulation

        python3 wormpop.py

    If you wish to cancel the simulation early

        Ctrl + Z
---------------------------------------------------------------------------------------------------------
How to use the constants.json commandline argument for the simulation:

    python3 wormpop.py --parameters=constants.json

    OR

    python3 wormpop.py < constants.json

        All simulation parameters are set to the defaults in constants.json

    Full commandline argument

        python3 wormpop.py [--parameters=<string>] [--database=<string>] [--name=<string>] [--directory=<string>]

            --database specifies the name of the SQLITE database to use. By default, this is :memory:, so it writes to a in memory sqlite3 database

            --name allows you to store the name of the simulation run. By default the value is "Simulation". The name is stored in the sqlite3 database. This is used mostly for metadata purposes

            --directory specifies the output of the TSV file
