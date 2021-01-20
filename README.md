# CMOR Mass Index
This script allows the user to select single-station radar observations of a meteor population, filter them, and generate a fit to the amplitude distribution in order to produce a population mass index. 
The user must provide the parameters of the desired population (date/time, radiant, velocity, certain thresholds, etc) and must configure how the functions are called depending on the goal. Examples are provided in the file, for comparing a shower to non-shower activity, and for observing the mass index over time.

Note: example data has been removed.

### Requirements:
- WMPL tools: The package and instructions can be found [here](https://github.com/wmpg/WesternMeteorPyLib)

- Multinest: Refer to [Pokorný and Brown (2016)](https://www.aanda.org/articles/aa/abs/2016/08/aa28134-16/aa28134-16.html).
The link in that article currently throws a 404 error. The source code is given in the zip file. This must be extracted to the Multinest directory and built according to its own instructions, then overwrite the default `priors.in` file with the one provided. Multinest has additional required FORTRAN packages in order to compile. All I can say to those is to follow the error messages.
- This code has only been tested on Linux and I cannot guarantee (nor do I expect) that it would work on Windows. In particular, Multinest almost certainly will not work. Mac functionality is also unknown.