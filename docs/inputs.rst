.. This is a comment. Note how any initial comments are moved by
   transforms to after the document title, subtitle, and docinfo.

.. demo.rst from: http://docutils.sourceforge.net/docs/user/rst/demo.txt

.. |EXAMPLE| image:: _images/temp.png
   :width: 1em

**********************
Inputs
**********************



NJOY input files are able to call multiple modules in sequence, which can be used together to perform complicated manipulations of the nuclear data. Broadly speaking, an NJOY input file will be structured as follows.

.. literalinclude:: exampleInputs/genericInput
  :language: html

Here, the first module name is specified, and the inputs are divided into "cards" (i.e. lines) of inputs. Each module has its own number of possible card, and each card has its own number of inputs. 

Modules in NJOY are rarely called in isolation, as they require the output of some other module as input. LEAPR, however, can be called alone, as it does not require any auxiliary input files (just the input file is sufficient to call LEAPR). THERMR requires that an input PENDF be provided, which can originate from earlier NJOY modules. 

.. Simple :math:`\mbox{H}` in :math:`\mbox{H}_2\mbox{O}`

Simple
========================================================
This simple example is for instructional purposes only. Note that parentheses following the forward slash will indicate the card number.

..
  COMMENT: .. contents:: Table of Contents
   :emphasize-lines: 1


-------------------------------------------------------------------------------


**Card 1**

.. literalinclude:: exampleInputs/simple
   :language: html
   :lines: 1-4

First, we specify the module (THERMR). 

The first card is used to specify input and output files. Here, thermal scattering input (which is the output from LEAPR) is tape24, and the PENDF (which could come from RECONR and BROADR, for example) is in tape23. THERMR will write its output file to tape25.
  

.. glossary::
    nendf 
        Thermal scattering input, which is the output from LEAPR. nendf can be set to zero if no LEAPR run is performed prior to this calculation, but doing so does limit THERMR's capabilities. In the directory that NJOY is run, nendf should appear as "tape24"

    nin
        Input PENDF file, which can be obtained from (for example) RECONR and BROADR. This should be available in the NJOY bin directory as "tape23".

    nout
        Desired output file to which THERMR will write the final PENDF. Here, it's "tape25"




.. -------------------------------------------------------------------------------

.. .. literalinclude:: exampleInputs/simple_H_H2O
..    :language: html
..    :lineno-start: 0
..    :lines: 1-3


-------------------------------------------------------------------------------

**Card 2**

.. literalinclude:: exampleInputs/simple
   :language: html
   :lines: 1-13

.. glossary::
    matde 
        Material ID for the material on the ``nendf`` tape. This will be equal to the ``mat`` value defined in the corresponding LEAPR run (LEAPR card 4, value 1).

    matdp 
        Material ID for the material on the ``nin`` tape, and will match the material identification used in preceding modules. For example, if tape23 was created using the RECONR module, then ``matdp`` will equal RECONR's ``mat`` specification from (RECONR card 3, value 1).

    nbin
        Number of equiprobable angle bins to be used for inelastic scattering distributions.
        


