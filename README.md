# SVV-FD Repository for TA's

Welcome TA's, these are the instructions on how to install and use the repository of the Flight Dynamics assignment of SVV.

## Installation

These are the installation instructions for the SVV-FD repository. It is recommended to make use of PyCharm as your IDE to use this repository. Please see [Using PyCharm](#using-pycharm) and [Using Git](#using-git) for further instructions.

The word `terminal` is used to reference both the Unix/MacOS terminal and the Windows command prompt. Alternatively, if using PyCharm as the IDE, it refers to the terminal of the project repository.

Make sure you have an installation of `Python 3.8` (download the latest version [here](https://www.python.org/downloads/)). 
Other versions might work, but are not tested.

### PyCharm

1. Install PyCharm Professional [here](https://www.jetbrains.com/pycharm/download/)
2. Run PyCharm Community
3. Click on `Get from VCS`
4. On the next screen select all defaults should be selected already, but if not:  
  `Repository URL` on the left  
  `Version control: Git` on the right  
   copy and paste you SSH Git URL `https://gitlab.tudelft.nl/SVV21/test/tas-flight-dynamics/fd-code-2021.git` into the `URL` field:
5. Go to `File` > `Settings`
6. Go to `File` > `Settings` > `Project: FD-Code-2021`. Select `Python Interpreter` on the right and click the gear icon in the top right > `Show all` > `Add (+)` > `Virtualenv Environment` > `OK`.
   This creates a virtual environment: an isolated workspace for the repository in which you can install dependencies locally, without cluttering your system.
7. Open a terminal window from within PyCharm.  
8. Enter `pip install -r requirements.txt`, hit `[enter]` and all the required dependencies will be installed automatically.


### Git

1. Install Git [here](https://git-scm.com/downloads) and run Git Bash.
2. Type `ssh-keygen -t rsa -b 4096 -C "your_email@student.tudelft.nl"`.
3. Hit `[enter]` three times.
4. Copy your public key that is stored in `~/.ssh/id_rsa.pub`. This file is located in your user folder (for example `catC:\users\%username%\.ssh\id_rsa.pub`). You can open this file in any text editor. It is just a string. Copy everything inside it. If you can't find the file using your Windows Explorer, just type `cat ~/.ssh/id_rsa.pub` in the terminal.
5. Add your SSH key to your GitLab account by clicking the button `Add SSH key` that appears on the main page of the TU Delft GitLab website.
6. Now paste the file contents from your public key file `C:\users\%username%\.ssh\id_rsa.pub` into the SSH key field, Give it a title and click below on `Add SSH Key`.
7. This tells GitLab that this is your computer.

## Using the repository

The repository includes all the scripts and directories necessary for the TA's to grade the FD assignments. A brief description of the contents of each directory is shown below. Make sure to check all individual files and do not hesitate to correct bugs if you spot any. Please communicate to the supervisors if you wish to make any further improvements to the repository and always make proper use of Git (see the section [Using Git](#using-git)).

### Repository Architecture
* The `data` directory contains all the data sources that are necessary for the assignment. This includes the flight test measurements file `matlab.mat`, `fuelarms.txt` and `fuelmasses.txt` for the centre of gravity calculations, and the reference file of the post flight datasheet.
* The `mass_balance_report` directory contains the scripts and data necessary for the calculation of the centre of gravity.
* The `model` directory contains the main scripts of the repository: 
  * The scripts `symmetric.py` and `asymmetric.py` contain the specification of the numerical models, methods for verification purposes, and the improvements made to the simulation to validate results. All relevant plots are produced prior to the parameter tuning.
  * The script `Cit_par.py` contains all the relevant stability derivatives and coefficients as posted on Brightspace, with an addition of a few methods to parse and update the values from the flight test data.
  * The script `stationary.py` contains all the methods to get the aerodynamic and stability coefficients of the aircraft during the flight test from the static measurements. All relevant plots are produced with the option of adding the theoretical trends.
  * The script `eigen_methods.py` contains all the methods necessary to obtain all the eigenmotion parameters in the numerical models.
* The directory `plots` contains all the resulting plots from running the scripts `symmetric.py`, `asymmetric.py`, and `stationary.py`.
* The directory `thrust` contains all the necessary scripts and files to produce the `thrust.dat` files using the `thrust.exe` as provided on Brightspace for the thrust reduction calculations.
### Using Git

Make sure to use branches when modifying the repository to avoid making incorrect changes on `main`. or more information, check the [help section](https://gitlab.tudelft.nl/help) and the [PyCharm help section](https://www.jetbrains.com/help/pycharm/manage-branches.html).

Here are few examples of the commands that can be used in the `terminal`. 
```
git status                          (Check status of a file in repository)
git add name_file                   (Add file to Git)
git add .                           (Add all the files to Git)
git commit -m ‘commit message’      (Make a commit in local branch)
git push origin name_branch         (Pushes branch to origin)
git branch                          (See branches)
git branch name_branch              (Create new branch)
git checkout name_branch            (Switch to another branch)
git merge name_branch               (Merge branch to main while on main)
git branch -d name_branch           (Delete branch locally)
```

## Credits
Developers:
* Daniel Martini
* Andrei Badea
* Andrej Nikolajevic
* Anique Altena
* Kipras Paliušis
* Midas Gossye

Editor:
* Fernando Corte Vargas

Under supervision of the staff of the course AE3212-II - Simulation, Verification & Validation:
* Dr.ir. Wouter van der Wal
* Dr.ir. Julien van Campen
* Dr.ir. Alexander in ’t Veld

## License

See license [here](license).
