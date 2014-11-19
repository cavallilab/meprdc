Please inspect the almost scripts comments for details on how they work:
// this is a comment

A few practial notes:
Currently only one-bond RDCs are supported. Bondlengths are assumed to be fixed to idealized values.
Make sure all RDCs have are on the same scale. That is, if HA-CA RDCs and H-N RDCs have been measured in the same alignment media one of them should be scaled by the appropriate Dmax-ratio so their overall scale match. Be sure to account for differences in experimental uncertainty also, if applicable.

pmerdc - constraint class specifics:
set_dmax sets the initial value of the product of Dmax and the degree of alignment s. This can be kept constant if set_autoscale is set to false. A good guess is usually the maximally observed experimental scaled by 110% - going to larger values will also work but will take longer to converge - if set_autoscale is true.



