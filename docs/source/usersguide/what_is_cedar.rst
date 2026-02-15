What is Cedar?
==============

Cedar: Computational Environment for Dynamics of Advanced Reactors

Cedar is an open-source multiphysics modeling code used to investigate the
transients of nuclear reactor systems. It was initially created with
applications to space nuclear reactors in mind; however, the modeling
methodology can be applied to any reactor system. It contains many physics
models that can be coupled together to capture the tightly-coupled multiphysics
that drive system performance in reactor systems.

Motivation
----------

There are a plethora of modeling tools available for nuclear reactor systems and
each has its pros and cons. Some are limited in the scope of physics included.
Some are limited in terms of ease-of-use due to steep learning curves. Others
are difficult to obtain due to licensing requirements.

Cedar aims to address these issues by housing the majority of the models in a
light-weight, intuitive, extensible, open-source framework that includes many
examples and tutorials to help users in the early stages.

Fidelity
--------

Cedar is considered a bi-fidelity framework. In this context, bi-fidelity refers
to high-fidelity models for some areas (such as heat transfer and neutron
transport) and low-fidelity models for other areas (such as gas flow and
turbomachinery). However, it is anticipated that users can add alternative
models to increase fidelity as they see fit.

Can I Trust Cedar's Results?
----------------------------

Cedar includes benchmarks for all of its models. A benchmark includes some
reference data which is compared against the outputs of Cedar. This reference
data could be an analytical solution, outputs of another tool (model-to-model
benchmark), or experimental data. All benchmarks are accompanied by a report
that includes the problem statement, expected outputs, and error metrics.


