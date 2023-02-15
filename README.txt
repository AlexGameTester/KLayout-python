Things to keep in mind (and I wish a had more time to fix them): 
- no backward compatibility is provided. Hence, a lot of `Examples` content that has its history will not work from scratch.
- no test are covering `Examples`. Therefore, there is no up-to-date info about examples that won't work.
- Previous statements are true for any design that has not been updated for more than 2 weeks.
- Installation process is not automated. There are 9 step instruction below to actually run this package.

To successfully use this library you should:
0. Install KLayout v.0.26.8 (or later, but no guarantee) 
1. download this library from github/shamil777/KLayout-python
2. open KLayout then press 'F5' to open macros editor.
3. In left panel choose the "python" tab.
4. Right-click to empty space in this tab -> 'Add Location'. Enter path to the 'KLayout-python' directory.
5. Create KLAYOUT_PYTHONPATH environment variable (this variable is used by KLayout python interpreter in the same way PYTHONPATH is utilized by system-wide python interpreter).
	This step depends on your OS so google 'how to' do this step.
7. In KLayout's macro editor: Locate and launch "KLAYOUT_PYTHONPATH.py". Assign output to the KLAYOUT_PYTHONPATH environment variable.
8. Restart KLayout.
9. Open macro editor and try to launch some example from "Klayout-python/Examples" folder.
10. Congrats! Have fun. If not -> leave an issue.
