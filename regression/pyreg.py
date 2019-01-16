#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model

import argparse

parser = argparse.ArgumentParser(description='Process some command-line options')

parser.add_argument("-f", "--file", help="File stem to use when opening input file.")
parser.add_argument("-o", "--output", help="File stem to use when opening output file.")
parser.add_argument("-i", "--index", help="Sweep index.")

args = parser.parse_args()

# Read in training data
data = np.genfromtxt("regression/{file}.csv".format(file=args.file), dtype=float, delimiter=',')


# Define arrays of training data
y =  data[:,-1]
x = data[:,0]
X = x.reshape(-1,1)

# Generate test data
x_plot = np.linspace(224.17, 224.57, 500)
X_test = x_plot[:,np.newaxis]

# Plot data points
fig, (plt1, plt2) = plt.subplots(2, 1)

plt1.scatter(x,y, color='black', label="Fit Points")
plt1.set_title("Background vs Frequency")

plt1.set(ylabel="Background Level", xlabel="Frequency (MHz)")
plt2.set(ylabel="Background Level", xlabel="Frequency (MHz)")

colors = ["red", "blue", "orange", "gold", "teal"]

for count, deg in enumerate([3]):
    poly = PolynomialFeatures(degree=deg) # changed
    x_ = poly.fit_transform(X)

    clf = linear_model.LinearRegression()
    result = clf.fit(x_,y)
    coefs = clf.coef_
    intercept = clf.intercept_
    R2 = clf.score(x_, y, sample_weight=None)

    X_trans = poly.fit_transform(X_test)
    Y_predict = clf.predict(X_trans)
    y_predict = clf.predict(x_)

    coefs = np.insert(coefs, 0, intercept)
    coefs = np.delete(coefs, 1)
#    print(coefs)

with open("regression/{output}_coefficients.csv".format(output=args.output),"wrb") as outfile:
        np.savetxt(outfile, coefs.reshape(1, 4), delimiter=",", fmt="%.20f")
outfile.close()

#plt1.plot(X_test, Y_predict, color=colors[count], linewidth=2, label="n=%d, " % deg + "R^2 = %.5f" % R2)

#print("Predictors ({0}): ".format(deg), coefs)
#print("Intercept ({0}): ".format(deg), intercept)
#print("R^2: ({0})".format(R2))

#difference = np.subtract(y, y_predict)
#plt2.scatter(x, difference, color=colors[count], label="n=%d, " % deg)

#plt1.legend(loc='lower left') 
#plt.show()

