#----------------------------------------------------------------------------#
#                         CognitiveSpaceWork                                 #
#----------------------------------------------------------------------------#


Author: Paul Hughes

Readme version: 1.1

#---

Purpose:

#---


To show my ability in python, not just to show my ability to deal with orbital dynamics, but also to show
that I can quickly produce thoughtful applications, and utilize common industry packages such as pandas.I also 
want to show that I can produce readable, and maintainable code, althouth this code isn't the cleanest I admit

#---

Quick Start Guide:

#---

If the dependencies are installed, run the main.py program in your favorite IDE, or in command line.
The dependency folders were created by PDM so if you use an IDE, I believe the path has to be specified in the configuration,
alternatively if you use pdm you have to use pdm install

Once the program is begun, a gui produced using the built in python library tkinter should pop up. There you can specify the input
keplerian elements to be simulated. Since keplerian elements are defined in a number of ways, I decided to at least make the anomaly term 
interchangable. Select what variable you'd like and press simulate. Alternatively, I found some orbital elements for the ISS, so I thought it 
would be fun to have those. I believe the program will catch any sort of common invalid input. 

Then, you can select some perturbations. Given my exeriences I've coded orb dynamics programs over and over. I decided to add the 
accelerations due to the obliquity of the Earth, since they are among the largest for a LEO or GEO sat. The math comes from Vallado's Astrodynamics book.
I would have used a more accurate set of elements for propagation, except the problem statement explicitly asked for propagation in cartesian. For the same 
reason, I coded a simple runge kutta 4th order ode solver. Anything more would have been a little overkill for cartesian elements.

Once you hit submit, the app will then produce the track given the inputs. When it is finished the first windows close, and an 
output dialogue box will pop up. You can pick any two variables to plot against eachother, and then have the plot print or save in 
root folder of the program. It will plot the first variable from the dropdown as the x axis, with the second for the y.
Alternatively, you can use the 3-D model which allows you to produce a 3d plot of the track, with a heat map for speed. 
You can also add earth to the plot because I thought that was a fun addition, although earth is not accurate to the shape. There also is a save csv button
so that you can save all the data from the track, which should help if you want to do any checking.I would have customized labels or done more to make the graphs 
look pretty but this is already overkill.

#---

Additional Info:

#---

If you'd like to change the physical constants, or the rate at which data is recorded, go to the change constants file in utils/loading_constants.

#---

TroubleShooting/Caution

#---

If you input a character to the inputs that are expecting, the dialogue box that pops up restarts the kernel and the program, so just be aware of that.
Additionally the tkinter library I used is not especially thread safe so try to keep only one instance of the program running at a time. If the initial 
run config window pops up as multiple windows, simply restart. This is a bug of the library I'm using.

Additionally, my personal pc ended up being upgraded to windows 11, and the pdm library was acting strange in powershell. If installing the dependencies from 
pdm does not work, attempt to use anaconda or spyder to run the file.

#---

Future Work

#---

If I continue developing this project I would do some things different. For one I would probably seperate the functions from the config files into seperate folders,
and work on the readability. Some of that problem came from spending only a small amount of time on this. I also would add the options to use different orbital 
element sents like equinoctial elements, additional perturbation terms, and an interactive version of the 3D graph export. I also probably would work on making the gui more stable first, and then start adding some testing. I realize this was only supposed to be a homework problem that should only take a couple hours, but I 
got carried away as I wanted to impress, and go a little above and beyond to show my worth.
