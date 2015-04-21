*******************************************************
PyChemkin - A Python interface for working with CHEMKIN
*******************************************************

Introduction
============

PyChemkin provides an easy interface for running CHEMKIN automatically using
Python scripts.  It writes input files to the CHEMKIN client, retrieves output,
and also parses data.  

Currently PyChemkin has capabilities for running homogeneous batch reactors, plug flow reactors,
and jet-stirred reactors.  It can calculate ignition delays automatically based on
output profiles and also run brute force ignition delay sensitivities.  

Usage
=====

Write the full path of your chemkin installation folder into the ``chemkin_path`` file.  i.e.::

/home/reaction/chemkin15131_linuxx8664

You can run any python script in the ``examples`` folder to test PyChemkin.

License
=======

Copyright (c) 2015 by Connie W. Gao (connie.w.gao@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the 'Software'),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

Installation
============

On any linux system, running PyChemkin simply requires adding the PyChemkin directory to the
PYTHONPATH like so::

PYTHONPATH=PyChemkin_directory:$PYTHONPATH

Then one can import the chemkin module inside any standard python script::

import chemkin