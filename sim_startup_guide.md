First Time Program and Python Setup on Windows:

    Open Windows terminal

    Type "python" to open up the microsoft store and click the "Get" button on the Python page to install Python for terminal usage

    Go to the directory that contains the OURE simulation (Refer to the Windows section of the startup guide if you do not know the commands for navigating directories inside of terminal)

    Type "python3 wormpop.py --parameters=constants.json" to attempt to run the program

    You should then be told by the terminal that you lack some of the modules used in the simulation

    Use the command "pip3 install NameOfTheModule" to install the module needed through terminal. Repeat this process until the terminal stops telling you that you are missing a module

    Once all of the modules it asks for are installed using the above method, the program should be able to run on your profile

-------------------------------------------------------------------------------------------------------
Windows:

    Press the Windows key.

    Type "cmd" or "PowerShell" and press Enter to open the Command Prompt or PowerShell.

    Use the ‘cd’ command to navigate to the directory where your wormpop.py script is located. For example, if your script is in C:\Example\Directory\Path, you can use the following commands in Command Prompt:

        To move through all directory paths at once:

            cd Example\Directory\Path

        To move through one directory at a time:

            cd Example

            cd Directory

            cd Path

       (You can also look inside the directory you are currently in to see what is inside it and what the next directory path you have to navigate to is called with the "ls" or "dir" command)

       If you wish to go back a directory:

           cd ..

    You can use a text editor like Notepad or a code editor like Visual Studio Code to open wormpop.py. For example, to open it with Notepad, you can run the following command in Command Prompt:

            notepad wormpop.py  

            (Remember: you cannot run a simulation on notepad)

    To run the simulation on a python file, you need to run the following command in Command prompt:

            python3 wormpop.py 

            or

            py wormpop.py

    Check and uncomment the code:

        Use your text editor to locate the last line of code and check if the lines "test = Simulation ('speed test')" and "test.run()" are commented out (prefixed with a #). If they are, simply delete the # to uncomment them.

        Please refer to the "How to use the constants.json commandline argument for the simulation" section of the guide for simulation constant modification.

            Make sure to note the original values of the constants before changing them so that they can be changed back later.

        Save your changes by pressing “Ctrl + S”

    End the simulation:

        Press “Ctrl + C” to stop the simulation

----------------------------------------------------------------------------------------------------------
Linux:

    Open the Linux terminal

    To move through all directory paths at once: 

        cd Example\Directory\Path

    To move through one directory at a time:

        cd Example

        cd Directory

        cd Path

    (You can also look inside the directory you are currently in to see what is inside it and what the next directory path you have to navigate to is called with the "ls" command)

    If you wish to go back a directory

        cd ..

    Open wormpop.py

        nano wormpop.py

    Go to the last lines of code and check to see if "test = Simulation('speedtest')" and "test.run()" are commented out. (A "#" before them signifies that they are commented out.)

        If so, delete the "#" to uncomment them so that the simulation can run.

    Please refer to the "How to use the constants.json commandline argument for the simulation" section of the guide for simulation constant modification.

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

    To modify constant values in the constants.json file:

        Windows:

            notepad constants.json

            Modify any constant values you wish to have ran differently in the simulation

            Ctrl + S to save the modified code

        Linux:

            nano constants.json

            Modify any constant values you wish to have ran differently in the simulation

            Close the modified file with Ctrl + X | When asked "Save modified buffer" hit Y for yes or N for no

