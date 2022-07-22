# I am a python code

print("Welcome to Seismology computational laboratory introduction")
print("Code: A simple calculation to find the finite difference coefficients for fourth order approximation of the second derivative of a function")

# Initiate some packages we will eventually need:
import numpy as np
#import matplotlib.pyplot as plt


m_list = [[1,1,1,1,1], [2,1,0,-1,-2], [4,1,0,1,4], [8,1,0,-1,-8], [16,1,0,1,16]]

A = np.array(m_list)

inv_A = np.linalg.inv(A)

dx = 1.
s = np.array([0,0,2./(dx)**2,0,0])
w = np.linalg.inv(A).dot(s)

print("Matrix A:")
print(A)

print("Inverse of A")
print(inv_A)

print("Expected result (y=s):")
print(s)

print("Finding the coefficients(x=w):")
print(w)
