#! /usr/bin/env python


def pascal(n):

    line = [1]
    for k in range(n):
        line.append(line[k] * (n-k) / (k+1))

    return line

length = 4
non_unique = length ** 2

print pascal(length)

choices = ['A', 'B']

for i in range(non_unique):
    sequence = ''
