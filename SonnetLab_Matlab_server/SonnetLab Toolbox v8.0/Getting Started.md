# Getting Started with SonnetLab

Starting with SonnetLab is simple. If you haven't installed it yet, though, **none of this will work**.
To install it, just go to `Set Path` in the `ENVIRONMENT` section in the `HOME` tab at the top of the matlab window. Then click *Add with Subfolders...* and select the *Scripts* folder in your release of SonnetLab.  
You can also look in the *Documentation*, *Tutorials*, or *Examples* folder if you want more information on any SonnetLab functionality, or need more detail on the installation process.   


The first thing you'll want to do is open a Sonnet Project.
You can open a project that already exists like this:
```
aProject = SonnetProject('Welcome.son')
```
Or open a new project like this:
```
aNewProject = SonnetProject()
```

If you want to save either project to a new file, you can do so like this:
```
aProject.saveAs('NewFile.son')
```
and if your project already has a file associated with it (because you used `saveAs` or opened it from a file like with `SonnetProject('Welcome.son')`)
you can save any changes like this
```
aProject.save()
```

### Making Changes
One of the simplest things to do is to change the frequency sweep settings of a project.
For example, you can add an ABS frequency sweep from 1 to 3 GHz (or whatever the project's units are, GHz are the default. This document assumes GHz for all frequency units)  
like this:
```
aProject.addAbsFrequencySweep(1, 3)
```
A sweep can also be created on its own with `CreateSweep`, which is recommended for projects using sonnet version 17 and up.
In fact, this method will not work for projects versions of Sonnet before 17.
This method is recommended because it is then easier to keep track of and change the sweep, and easier to add/remove it from sweep sets  
A sweep can be created like so:
```
aSweep = CreateSweep('ESWEEP', 'Y', 'N', 1, 20, 100) %Create an exponential sweep, enabled for simulation, not adaptive, with 100 analasys frequencies, going from 1GHz to 20GHz
```
The options for the type (case insensitive), here `ESWEEP`, are `SWEEP`, `ABS_ENTRY`, `ABS` (identical to `ABS_ENTRY`), `DC`, `ESWEEP`, `LSWEEP`, `SIMPLE`, and `STEP`  
The new sweep can now be added to the project like this:
```
aIndex = aProject.addSweep(aSweep) % the value of aIndex is the index of this sweep in the current sweep set
```
The new sweep can be added to a particular sweep set (by index) like this:
```
aIndex = aProject.addSweepToSet(1, aSweep) % add aSweep to sweep set at index 1
% the value of aIndex is the index of this sweep in the sweep set at index 1
```
New sweep sets can be created like this:
```
aSweepSetIndex = aProject.addSweepSet()
aSweepSetIndex2 = aProject.addSweepSet('N') %Add a new sweep that is not enabled
% aSweepSetIndex is the index of the newly created sweep set
```
The handle to those sweep sets can be retrieved like this:
```
aSweepSet = aProject.getSweepSet(aSweepSetIndex)
```
And new sweeps to the sweep set
```
aSweepIndex = aSweepSet.addSweep('LIST', 'Y', 'N', [10, 20, 30]) % a new list sweep, enabled, not adaptive, with frequencies of 10, 20, and 30 GHz
aSweepIndex2 = aSweepSet.addSweep(CreateSweep('LIST', 'Y', 'N', [40, 50, 60])) % a new list sweep, enabled, not adaptive, with frequencies of 40, 50, and 60 GHz
aSweepIndex3 = aSweepSet.addSweep(CreateDefaultSweep('LIST', [70, 80, 90])) % a new list sweep, (enabled, not adaptive by default), with frequencies of 70, 80, and 90 GHz
% The new sweeps are now at indices aSweepIndex, aSweepIndex2, and aSweepIndex3 inside the set aSweepSet
```
And to get the sweeps back from the sweep set
```
aListSweep = aProject.getSweepFromSet(aSweepSetIndex, aSweepIndex)
```
And now the sweep can be changed and still impact the project, because it is a handle
```
aListSweep.setList([1, 2, 3]) % Set the list of frequencies to 1, 2, and 3 GHz
```
If you now save the project, using `aProject.save()`, the list of frequencies will be 1, 2, and 3 GHz


### Viewing the project

The project can be prieviewed in several ways. 
Using `aProject.draw3d()` or `aProject.drawCircuit()` will show a 3D preview of the circuit in matlab. `aProject.draw2d(layer_index)` or `aProject.drawLayer(layer_index)` will
display a top-down 2d prieview of the circuit layer `layer_index`. Finally, a text preview of all polygons is possible with `aProject.displayPolygons()`, displaying, for example:
```
#     ID      Centroid        Mean Point    Size       Type       Level      Metal Type
-----------------------------------------------------------------------------------------
1     222  	(180  ,  125)	(181  ,  123)	150  	          	    0       Lossless
2     223  	(90   ,   85)	(108  ,   84)	190  	          	    0       Lossless
3     225  	(270  ,  165)	(254  ,  166)	170  	          	    0       Lossless
```

### Units

`aProject.changeFrequencyUnit('MHZ')` changes the frequency unit to megahertz  
In all the following functions `unit_string` is case insensitive (i.e. attempting to use *millihertz* for frequency will give you *megahertz*)
`aProject.changeFrequencyUnit(unit_string)` where `unit_string` is one of `Hz`, `kHz`, `MHz`, `GHz`, `THz`, `PHz`  
`aProject.changeInductanceUnit(unit_string)` where `unit_string` is one of `H`, `mH`, `uH`, `nH`, `pH`, `fH`  
`aProject.changeLengthUnit(unit_string)` where `unit_string` is one of `MIL`, `um`, `mm`, `cm`, `in`, `m`  
`aProject.changeAngleUnit(unit_string)` the only supported unit is `DEG` 
`aProject.changeResistanceUnit(unit_string)` where `unit_string` is one of `OH`, `kOH`, `MOH`  
`aProject.changeCapacitanceUnit(unit_string)` where `unit_string` is one of `F`, `mF`, `uF`, `nF`, `pF`, `fF`  


### Geometry and Netlist projects

You can check whether a project is a geometry or netlist project by using `aProject.isNetlistProject()` or `aProject.isGeometryProject()`  
Projects creates with `aProject = SonnetProject()` are by default geometry projects
```
aBool = aProject.isNetlistProject()
aBool = aProject.isGeometryProject()
```

### Simulating and Working with Sonnet on your PC

It is possible to run simulations and open projects in Sonnet from Matlab. Currently, this is only tested on Windows.
Simulating the project is as simple as running `aProject.simulate()`. It is also possible to pass options in a  string to `simulate()`, 
these are explained more in depth in the documentation comments in 'SonnetProject.m'. Estimate the memory usage of the simulation with `aProject.estimateMemoryUsage()`.
Provide a specific Sonnet path with `aProject.estimateMemoryUsage('C:\Program Files\sonnet.12.52')`. Delete the simulation data for a project with `aProject.cleanProject()`.
Current density calculations can be enabled with `aProject.enableCurrentCalculations()` and disabled with `aProject.disableCurrentCalculations`. (Be aware that current calculations can add significant simulation time to a project).


Open a project in the Sonnet GUI with `aProject.openInSonnet()`. The project will then be opened in Sonnet, and execution is blocked until Sonnet is closed.  
A non blocking call can be made with `aProject.openInSonnet(false)` (setting the `isWaitForGuiToClose` parameter to false), and opened with a particular version of sonnet, like 12.52, with `aProject.openInSonnet(false,'C:\Program Files\sonnet.12.52')`.

View response data in Sonnet with `aProject.viewResponseData()`. A path can again be specified with: `aProject.viewResponseData('C:\Program Files\sonnet.12.52')`.  
Similarly, view currents with `aProject.viewCurrents()`, and open with a particular version with `aProject.viewCurrents('C:\Program Files\sonnet.12.52')`.

The file contents corresponding to the project (what would be written using `aProject.save()`) can be retrieved with `aProject.stringSignature()`
