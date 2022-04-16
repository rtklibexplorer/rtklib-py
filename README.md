## **rtklib-py : A python subset of RTKLIB for PPK solutions**


### **What is rtklib-py?:**

Rtklib-py is a free and open source subset of RTKLIB written in Python. It was originally based on CSSRlib but has been rewritten to align very closely with the demo5 version of RTKLIB.  It currently supports only PPK solutions using GPS, GLonass, and Galileo constellations.  It supports most but not all options available in PPK solutions in RTKLIB and produces solutions very similar but not identical to the RTKLIB solutions.  Debug trace is supported and trace level 3 messages are aligned to closely match the demo5 RTKLIB messages.  Funtction names, comments, and variables are also closely aligned to the demo5 RTKLIB code.

### **How to use:** 

Run_ppk.py is the top level script to run a PPK solution.  The rtklib-py package includes two sample data sets, one is a u-blox F9P rover mounted on the roof of a car, the second is a data set from the 2021 Google Smartphone Decimeter Challenge.  Run_ppk.py is configured to generate a solution for the u-blox data set but includes commented out lines to run the GSDC example.  Config parameters are in config_f9p.py and config_phone.py and align closely to the config parameters in RTKLIB.  

### **Purpose:** 

This code is not meant to be a replacement for RTKLIB but as a tool to either:
  1)  Use as a "map" to explore the inner details of how RTKLIB works
  2)  Use as a development environment to experiment with enhancements or adjustments to the RTKLIB algorithms.  Due to the close alignment between the two packages, these can then be fairly easily ported to the C/C++ version of RTKLIB.  In particular it is hoped that it will be useful to competitors in the Google Smartphone Decimeter Challenge.
  3)  Pieces of this code can be cut and paste into more custom solutions
