Constant Material
=================

``cedar.materials.ConstantMaterial()``

Unlike the other material properties objects which represent a specific
material, the ``ConstantMaterial`` is used to represent any material with
constant values provided by the user. In the data shown below, the object is
instantiated with all values equal to 1.

Density
-------

1 [kg/m3]

Thermal Conductivity
--------------------

Valid Temperature Range:

.. math::
   0 \le T \le 10000

.. figure:: img/constantmaterial_k.png
   :scale: 100 %
   :align: center

   Thermal Conductivity

Specific Heat Capacity
----------------------

Valid Temperature Range:

.. math::
   0 \le T \le 10000

.. figure:: img/constantmaterial_cp.png
   :scale: 100 %
   :align: center

   Specific Heat Capacity