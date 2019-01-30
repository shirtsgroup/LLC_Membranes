#!/usr/bin/env python

import sqlite3 as sql
import names

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()
penalty = 0.25
query = "SELECT name, sigma, alpha, hurst, python_MSD, python_MSD_CI_lower, python_MSD_CI_upper from msd WHERE penalty = %.2f ORDER BY python_MSD DESC" % penalty
output = crsr.execute(query).fetchall()

table = ""
table += r"\begin{table}[h]" + "\n" + r"  \centering" + "\n" + r"  \begin{tabular}{cccc}" + "\n" + r"  \toprule" + "\n"
table += r"  System & $\sigma$ ($nm$) & $\alpha$ & $H$ \\" + "\n" + r"  \midrule" + "\n"
#table += r"  System & $\sigma$ ($nm$) & $\alpha$ & $H$ & Simulated MSD ($nm^2$)\\" + "\n" + r"  \midrule" + "\n"

#for i in output:
#	table += r"  %s & %.2f & %.2f & %.2f & %.2f [%.2f, %.2f] \\" % (names.res_to_name[i[0]], i[1], i[2], i[3], i[4], i[5], i[6]) + "\n"

for i in output:
	table += r"  %s & %.2f & %.2f & %.2f \\" % (names.res_to_name[i[0]], i[1], i[2], i[3]) + "\n"


table += r"  \bottomrule" + "\n" + r"  \end{tabular}"
print(table)
