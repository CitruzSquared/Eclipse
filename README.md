# Eclipse
Eclipse calculator programs

# In common for both:
All times, unless otherwise mentioned, is given in the unit of days.
Day 0 represents a specific northern vernal equinox. (That must occur at midnight)

# "eclipse_date_predictor.py"
"eclipse_date_predictor.py" is a python file.

**How to use:**
Open it and change the variable values to the values of your world, then run it.
*Make sure to follow the comments*
Saros mode will give dates for every eclipse in the specific saros cycle of the member eclipse entered.

Solar eclipses are given in the form: "[A, B, C, D, E] NM: X".
A, B, C, D, E represent approximations of the times of P1, U1, Greatest Eclipse, U2, P2 respectively.
"-" means that particular contact did not happen.
X is the time of New Moon. (conjunction of the Moon and Sun in Ecliptic Longitude)

Lunar eclipses are given in the form: "[A, B, C, D, E, F, G] FM: X, D: Y, M: Z".
A, B, C, D, E, F, G represent approximations of the times of P1, U1, U2, Greatest Eclipse, U3, U4, P4 respectively.
"-" means that particular contact did not happen.
X is the time of Full Moon. (opposition of the Moon and Sun in Ecliptic Longitude)
D is the duration of the entire eclipse in days.
M is the magnitude of eclipse.

# "solar_eclipse_visualizer"
Download the folder.
"solar_eclipse_visualizer.pyde" is a processing 3 file in python mode.
https://py.processing.org/tutorials/gettingstarted/ is the tutorial for getting it. I recommend processing 3, not processing 4.
*The folder and the .pyde file must have the same name.*
"sketch.properties" is needed for the .pyde file to run.

**How to use:**
Copy all the variables from "eclipse_date_predictor.py".
One extra variable is required, earth_hours_in_day, and it should be set according to ur world.
The options showGrid and longitudeCorrect can be turned on or off.
The resolution of the image can be changed.
*Read the comments*
To generate a map for a specific eclipse, select a NM value generated from "eclipse_date_predictor.py", and put it as the value as the variable "time" on line 2.
Then you can run the program, and it will automatically save a file called "X.png" where X is the time value entered into the folder.

The labeled diagram for the map generated can be found here: 
![image](https://user-images.githubusercontent.com/23460281/210059411-89fce457-80ae-47ab-a3c7-c409d2699da0.png)
If the path of totality is shown in dark red instead of green, it means the eclipse is annular, and the Moon is not big enough to cover the sun.

The program breaks sometimes and I do not know why.
