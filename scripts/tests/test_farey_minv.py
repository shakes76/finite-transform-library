'''
Test for the pure independent modules for number theory
'''
import _libpath #add custom libs
import finitetransform.numbertheory as nt
import finitetransform.farey as farey

x = 5
m = 641

##########
# test the multiplicative inverse mod m
inv = nt.minverse(x, m)
print inv

identity = (inv*x)%m
print identity

##########
# test the Farey finite vector function
vector = farey.farey(1, 2)
m = farey.toFinite(vector, m)
print vector, "->", m
