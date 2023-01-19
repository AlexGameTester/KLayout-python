Design files were generated with Klayout-python. Klayout-python is a Python library developed by AQS lab (MIPT) for 
automating the design of superconducting quantum circuits. It uses the KLayout[1] layout design program API. 
KLayout itself allows to write and execute arbitrary python code through embedded interpreter, and primary aimed at 
designing 2D structures for various application domains.
KLayout-python design code is executed within the "KLayout" infrastrucre and specializes on designing microwave and 
superconducting qubit planar designs, including drawing patterns, simulation and domain specific design rules checkers.
There are other frameworks for quantum structures designing based upon KLayout [2-4] that are maintained by the scientific community.

[1] - https://www.klayout.de/
[2] - https://iqm-finland.github.io/KQCircuits/ - also aimed at superconducting qubit's designing
and considered superior to Klayout-python.
[3] - https://klayoutmatthias.github.io/xsection/ - simulates deposition process through any of crossections of the design.
[4] - https://github.com/siepic/SiEPIC-Tools/blob/master/README.md - semiconductor photonics designing


To successfully use this library you should:
1. download this library from github/shamil777/KLayout-python
2. open KLayout, press 'F5' to open the Macro Editor.
3. In the leftmost panel choose the "Python" tab.
4. Right-click to the leftmost panel's empty region and choose 'Add Location'. 
Enter path to the 'KLayout-python' directory.
5. Create KLAYOUT_PYTHONPATH environment variable
(this variable is used by KLayout python interpreter). 
This step depends on your OS so google 'how to' do this step.
7. In KLayout's macro editor: Locate and launch KLAYOUT_PYTHONPATH.py. Assign output to the KLAYOUT_PYTHONPATH environment variable.
8. Restart KLayout. 
9. Open macro editor and try to launch some example from "Klayout-python/Examples" folder.
10. Congrats! Have fun. If not -> leave an issue.

For developers:
1. Docstring format - numpy
2. Read commit_template.txt
3. If you created a new primitive, please attach its basic usage
example into examples folder. Furthermore, you may schematics with
geometry parameters specified as dimensions (see ./classLib/schematics/ for example).
4. If you meet any troubles in contributing - you may contact me directly via
"kadyrmetov@phystech.edu"