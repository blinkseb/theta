.. _overview:

********
Overview
********

``theta`` can be used in different ways, or using different *layers*:
 #. The *plugin layer*: extending the theta functionality with plugins -- such as adding new template morphing implementation, or adding new statistical methods
 #. The *cfg-file layer*: the main theta program written in C++ operates on text-based configuration files which specify the statistical model and the statistical methods to run. The configuration file format is relatively simple, but the configuration tends to becomes very large for more complex models.
 #. The *theta-auto layer*: `theta-auto` is a set of python scripts which automatize the task of creating the configuration files for the *cfg-file layer* and also collect and display the result.

In the order the layers were listed, the power and flexibility in general decreases, but the usability increases. 
For most users, *theta-auto* is best suited, and it is recommended that new users start with *theta-auto*.

Documentation overview
-----------------------

There are several suources of documentation for ``theta``, listed with decreasing relevance:
 * The documentation on these pages cover *theta-auto*, which should be enough for most users.
 * A general, physics-focused introduction is available as pdf `here <http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/theta.pdf>`_. It includes some documentation of the cfg-file layer, but does not mention *theta-auto* yet.
 * For a documentation of the C++ code (including the cfg-file format), please refer to the `doxygen documentation <http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/testing/html/index.html>`_.
 * The source code is also available online via `trac <https://ekptrac.physik.uni-karlsruhe.de/trac/theta/>`_.

Note that the source code contains a lot more documentation than this document at many places already.




